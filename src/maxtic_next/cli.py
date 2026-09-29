"""命令行入口（argparse 封装，参数名等价原版 ``MaxTiC.py``）。

原版裸 ``sys.argv`` 参数映射：

* ``species_tree`` / ``constraints``：位置参数（``constraints`` 在 ``--from-ale`` 模式下
  接受一个或多个 ``.uml_rec`` 文件）；
* ``ls=LOCAL_SEARCH`` → ``--local-search`` / ``--ls``；
* ``t=TEMPERATURE`` → ``--temperature`` / ``--t``；
* ``r=RANDOMIZATION_TYPE`` → ``--random-type`` / ``--r``；
* ``d=MIN_DISTANCE`` → ``--min-transfer-distance`` / ``--d``（按 phylogenetic distance 过滤）；
* ``ts=THRESHOLD_CONSTRAINTS`` → ``--threshold-constraints`` / ``--ts``；
* ``rd=RANDOM_TREES`` → ``--random-trees`` / ``--rd``；
* 额外新增 ``--seed``（默认 42，保证无参可复现）。

可选插入点（默认不影响核心排序结果）：
* ``--from-ale``：把 ``.uml_rec`` 文件列表经 ALE 适配器转为约束；
* ``--from {ale,ranger,eccetera,artra,alerax,auto}`` / ``--from-auto``：多软件兼容，
  把 constraints 视为对应上游工具输出（或 auto 自动检测），经注册表转为统一约束；
* ``--target-clade <node>`` / ``--target-clade-ancestor-map``：目标类群剪裁；
* ``--dry-run``：正式排序前静态预检；退出码取自结构化 ``report["ok"]``
  ，合法输入返回 0；
* ``--no-html``：关闭默认开启的交互式 HTML 报告；
* ``--ale-min-support`` / ``--ale-min-family-size`` / ``--ale-cache-dir`` / ``--ale-source``：
  ALE 适配器参数。``--ale-source`` 默认 **trf**（ALE 官方 MaxTiC
  集成口径）；
* ``--min-endpoint-hit-rate`` / ``--quiet-adapters``：适配器可观测性开关
  （端点命中率告警阈值 / 抑制 stderr 例行摘要）。

取值域校验：所有数值/枚举选项都在 **argparse 层** 用 ``type`` /
``choices`` 卡住取值域，非法输入退出码 2 并给出可操作的提示；旧的 ``type=int`` /
``type=float`` 会让 ``--r 3``、``--ts 5``、``--ls -5`` 之类拼错静默通过。
"""

import argparse
import math
import sys
from typing import List, Optional

from maxtic_next import __version__ as PACKAGE_VERSION
from maxtic_next import api
from maxtic_next.config import (
    DEFAULT_MCMC_ITERS,
    DEFAULT_MIN_ENDPOINT_HIT_RATE,
    DEFAULT_MIN_TRANSFER_DIST,
    DEFAULT_NEAR_OPTIMAL_TOP_K,
    DEFAULT_OUTPUT_STYLE,
    DEFAULT_RANDOM_TREES,
    DEFAULT_RANDOM_TYPE,
    DEFAULT_SEED,
    DEFAULT_TEMPERATURE,
    DEFAULT_THRESHOLD_CONSTRAINTS,
    DEFAULT_TIME_FOR_SEARCH,
    MCMC_TEMPERATURE_AUTO,
    RANDOM_TYPE_CHOICES,
    THRESHOLD_CONSTRAINTS_MAX,
    THRESHOLD_CONSTRAINTS_MIN,
)


def _nonneg_int(name: str):
    """argparse type：非负整数。"""

    def _check(raw: str) -> int:
        try:
            val = int(raw)
        except (TypeError, ValueError):
            raise argparse.ArgumentTypeError(f"{name} 需为整数，收到 {raw!r}") from None
        if val < 0:
            raise argparse.ArgumentTypeError(f"{name} 必须 >= 0，收到 {val}")
        return val

    return _check


def _pos_int(name: str):
    """argparse type：正整数。"""

    def _check(raw: str) -> int:
        try:
            val = int(raw)
        except (TypeError, ValueError):
            raise argparse.ArgumentTypeError(f"{name} 需为整数，收到 {raw!r}") from None
        if val <= 0:
            raise argparse.ArgumentTypeError(f"{name} 必须 > 0，收到 {val}")
        return val

    return _check


def _nonneg_float(name: str):
    """argparse type：非负有限浮点数。"""

    def _check(raw: str) -> float:
        try:
            val = float(raw)
        except (TypeError, ValueError):
            raise argparse.ArgumentTypeError(f"{name} 需为数值，收到 {raw!r}") from None
        if math.isnan(val) or math.isinf(val):
            raise argparse.ArgumentTypeError(f"{name} 必须是有限数值，收到 {raw!r}")
        if val < 0:
            raise argparse.ArgumentTypeError(f"{name} 必须 >= 0，收到 {val}")
        return val

    return _check


def _pos_float(name: str):
    """argparse type：正有限浮点数。"""

    def _check(raw: str) -> float:
        try:
            val = float(raw)
        except (TypeError, ValueError):
            raise argparse.ArgumentTypeError(f"{name} 需为数值，收到 {raw!r}") from None
        if math.isnan(val) or math.isinf(val):
            raise argparse.ArgumentTypeError(f"{name} 必须是有限数值，收到 {raw!r}")
        if val <= 0:
            raise argparse.ArgumentTypeError(f"{name} 必须 > 0，收到 {val}")
        return val

    return _check


def _unit_interval(name: str):
    """argparse type：``[0, 1]`` 闭区间内的比例。"""

    def _check(raw: str) -> float:
        try:
            val = float(raw)
        except (TypeError, ValueError):
            raise argparse.ArgumentTypeError(f"{name} 需为数值，收到 {raw!r}") from None
        if not (THRESHOLD_CONSTRAINTS_MIN <= val <= THRESHOLD_CONSTRAINTS_MAX):
            raise argparse.ArgumentTypeError(
                f"{name} 是比例，取值域 [{THRESHOLD_CONSTRAINTS_MIN}, "
                f"{THRESHOLD_CONSTRAINTS_MAX}]，收到 {val}"
            )
        return val

    return _check


def _mcmc_temperature(raw: str) -> float:
    """argparse type：``auto``（默认，哨兵 0.0）或正温度。"""
    if raw.strip().lower() in ("auto", ""):
        return MCMC_TEMPERATURE_AUTO
    try:
        val = float(raw)
    except (TypeError, ValueError):
        raise argparse.ArgumentTypeError(
            f"--mcmc-temperature 需为 'auto' 或正数，收到 {raw!r}"
        ) from None
    if math.isnan(val) or math.isinf(val) or val <= 0:
        raise argparse.ArgumentTypeError(f"--mcmc-temperature 需为 'auto' 或正数，收到 {raw!r}")
    return val


def build_parser() -> argparse.ArgumentParser:
    """构建参数解析器。"""
    parser = argparse.ArgumentParser(
        prog="MaxTiC-Next",
        description=(
            "MaxTiC-Next —— MaxTiC 的 Python 3 重写：基于最大时间一致性（MTC）"
            "对物种树内部节点排序，以兼容水平基因转移（HGT）推导的加权时间约束。"
        ),
    )
    parser.add_argument("--version", action="version", version=f"MaxTiC-Next {PACKAGE_VERSION}")
    parser.add_argument(
        "species_tree",
        help="物种树文件路径（Newick，内部节点标签在 bootstrap 字段；"
        ".gz / gzip 流与单成员 .tar.gz 自动解压）",
    )
    parser.add_argument(
        "constraints",
        nargs="+",
        help="约束文件路径（空格或逗号双格式；.gz / gzip 流与单成员"
        " .tar.gz 自动解压，多成员归档会给出 tar 解包指引）；"
        "使用 --from-ale 时为一个或多个 .uml_rec 文件",
    )
    parser.add_argument(
        "--seed",
        type=_nonneg_int("随机种子"),
        default=DEFAULT_SEED,
        help=f"随机种子，驱动 mix 平局与局部搜索（默认 {DEFAULT_SEED}；须 >= 0）",
    )
    parser.add_argument(
        "--local-search",
        "--ls",
        dest="local_search",
        type=_nonneg_float("局部搜索时长"),
        default=DEFAULT_TIME_FOR_SEARCH,
        help=f"局部搜索时长（秒），0 表示关闭（默认 {DEFAULT_TIME_FOR_SEARCH}；须 >= 0）",
    )
    parser.add_argument(
        "--temperature",
        "--t",
        dest="temperature",
        type=_pos_float("Metropolis 温度"),
        default=DEFAULT_TEMPERATURE,
        help=f"Metropolis 温度（默认 {DEFAULT_TEMPERATURE}；须 > 0）",
    )
    parser.add_argument(
        "--random-type",
        "--r",
        dest="random_type",
        type=int,
        choices=list(RANDOM_TYPE_CHOICES),
        default=DEFAULT_RANDOM_TYPE,
        help="随机化类型 0/1/2（默认 0；其它取值会被拒绝——"
        "旧实现把 3/-1 等静默当作 0，等于用未随机化的数据做对照）",
    )
    parser.add_argument(
        "--min-transfer-distance",
        "--d",
        dest="min_transfer_distance",
        type=_nonneg_float("最小转移距离"),
        default=DEFAULT_MIN_TRANSFER_DIST,
        help="最小转移距离阈值，按 phylogenetic distance 列过滤"
        f"（默认 {DEFAULT_MIN_TRANSFER_DIST}；须 >= 0。"
        "注意：输入无距离列时该选项被忽略）",
    )
    parser.add_argument(
        "--threshold-constraints",
        "--ts",
        dest="threshold_constraints",
        type=_unit_interval("约束权重阈值比例"),
        default=DEFAULT_THRESHOLD_CONSTRAINTS,
        help=f"约束权重阈值比例，取值域 [{THRESHOLD_CONSTRAINTS_MIN}, "
        f"{THRESHOLD_CONSTRAINTS_MAX}]"
        f"（默认 {DEFAULT_THRESHOLD_CONSTRAINTS}）",
    )
    parser.add_argument(
        "--random-trees",
        "--rd",
        dest="random_trees",
        type=_nonneg_int("随机树采样数量"),
        default=DEFAULT_RANDOM_TREES,
        help="随机树采样数量，0 表示关闭（默认 0；须 >= 0）",
    )

    # ---- ALE 适配器----
    parser.add_argument(
        "--from-ale",
        action="store_true",
        help="将 constraints 位置参数视为 ALE .uml_rec 文件列表，经适配器转为约束",
    )
    parser.add_argument(
        "--ale-min-support",
        type=_unit_interval("最小支持度"),
        default=0.05,
        help="ALE 适配器：单基因家族内最小支持度阈值（默认 0.05）",
    )
    parser.add_argument(
        "--ale-min-family-size",
        type=_nonneg_int("最小基因家族规模"),
        default=5,
        help="ALE 适配器：最小基因家族规模（调和树叶子数，默认 5）",
    )
    parser.add_argument(
        "--ale-cache-dir", default=None, help="ALE 适配器：基因家族文件级缓存目录（断点续传）"
    )
    parser.add_argument(
        "--ale-source",
        choices=["rec", "trf"],
        default="trf",
        help="ALE 适配器：约束来源。trf=**默认**，按转移事件生成"
        "（ALE 官方 MaxTiC 集成 constraints_from_transfers 的"
        "口径：donor 的父节点须早于 receptor 的子节点）；"
        "rec=按调和事件生成（与 ALE 的 .uml_rec 调和记录对应）。"
        "两者产出的约束集规模与权重量级**不同**，切换即换掉整个输入；"
        "实际口径见 stderr 摘要与 run_metadata['adapter_diagnostics']",
    )
    parser.add_argument(
        "--min-endpoint-hit-rate",
        type=_unit_interval("端点命中率阈值"),
        default=DEFAULT_MIN_ENDPOINT_HIT_RATE,
        help="适配器端点命中率阈值（默认 "
        f"{DEFAULT_MIN_ENDPOINT_HIT_RATE}）：上游输出的转移事件中"
        "能落进物种树标签集的比例低于该值时在 stderr 警告，"
        "为 0 时直接报错（不再静默产出空约束集）。"
        "设为 0 关闭告警（0 命中仍报错）",
    )
    parser.add_argument(
        "--quiet-adapters",
        dest="quiet_adapters",
        action="store_true",
        help="不打印适配器在 stderr 上的例行口径声明与诊断摘要"
        "（诊断量仍写入 Result.run_metadata，错误照旧抛出）",
    )
    parser.add_argument(
        "--ale-parallel",
        choices=["process", "thread"],
        default="process",
        help="ALE 适配器：并行解析模式 process=进程池（默认，绕过 GIL）/ thread=线程池",
    )

    # ---- 多软件兼容框架：统一 --from / --from-auto（ranger/eccetera/artra/alerax）----
    parser.add_argument(
        "--from",
        dest="from_tool",
        choices=["ale", "ranger", "eccetera", "artra", "alerax", "auto"],
        default=None,
        help="上游工具适配器：把 constraints 视为该工具输出文件/目录，"
        "经注册表分发转为统一约束；auto=按内容/扩展名自动检测"
        "（ALE 适配器参数 --ale-* 对所有调和适配器通用）",
    )
    parser.add_argument(
        "--from-auto",
        dest="from_auto",
        action="store_true",
        help="等价 --from auto：自动检测上游工具输出格式并分发",
    )
    parser.add_argument(
        "--artra-transfer-kind",
        choices=["all", "replacing", "additive"],
        default="all",
        help="ARTra 适配器：统计的转移类别（默认 all=加性+替换）",
    )

    # ---- 目标类群剪裁----
    parser.add_argument(
        "--target-clade", default=None, help="目标类群根节点标签；默认仅保留两端点均在类群内的约束"
    )
    parser.add_argument(
        "--target-clade-ancestor-map",
        action="store_true",
        help="开启'外部祖先映射到根'的保守剪裁规则（需配合 --target-clade；"
        "默认关闭，避免静默引入虚假偏序）",
    )

    # ---- 预检----
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="仅做静态预检（物种树格式/二叉性/标签唯一性/约束格式/"
        "权重合法性/端点合法性/d 列语义），不运行正式排序。"
        "退出码：0=通过，1=存在 error 级问题",
    )

    # ---- HTML 报告----
    parser.add_argument("--no-html", action="store_true", help="关闭默认开启的交互式 HTML 报告")

    # ---- 线性扩展上的 Metropolis-Hastings 采样（默认关闭）----
    parser.add_argument(
        "--mcmc",
        action="store_true",
        help="开启在物种树拓扑相容的线性扩展上的 Metropolis-Hastings 采样"
        "（默认关闭，不影响主流程结果）。注意：这是一条**尚未通过收敛"
        "诊断验证**的 MH 链，其样本不得当作后验样本使用",
    )
    parser.add_argument(
        "--mcmc-iters",
        type=_pos_int("MCMC 迭代次数"),
        default=DEFAULT_MCMC_ITERS,
        help=f"MCMC 链步数（默认 {DEFAULT_MCMC_ITERS}；必须 > 0）",
    )
    parser.add_argument(
        "--mcmc-temperature",
        type=_mcmc_temperature,
        default=MCMC_TEMPERATURE_AUTO,
        help="MCMC Boltzmann 温度：'auto'（默认）表示 "
        "max(总权重,1)/100，有效值会打印并记入 run_metadata；"
        "也可给出正数以固定温度运行",
    )
    parser.add_argument(
        "--mcmc-burn-in",
        type=_nonneg_int("MCMC burn-in"),
        default=None,
        help="MCMC burn-in 步数（默认取 --mcmc-iters 的 50%%），该段样本只推进链、不进入统计",
    )
    parser.add_argument(
        "--mcmc-thin",
        type=_pos_int("MCMC thin"),
        default=None,
        help="MCMC 抽稀间隔（默认 max(1, iters//100)），用于降低相邻样本的自相关",
    )

    # ---- 增量 scoring（默认开启）与检查点续跑（默认关闭）----
    # 增量与全量在 4 组输入 × 5 个种子 ×
    # {stdout, ranked newick, best_order, values, 三个 TSV} 上**逐字节一致**，
    # 实测提速 3.1-8.2 倍；`--no-incremental` 保留原版全量口径作为逃生门。
    parser.add_argument(
        "--incremental",
        dest="incremental",
        action="store_true",
        default=True,
        help="增量 scoring（**默认已开启**）：局部搜索每次区间旋转仅算 "
        "增量 O(b-a)，而非全量 O(|E|)。此开关保留用于显式表达",
    )
    parser.add_argument(
        "--no-incremental",
        dest="incremental",
        action="store_false",
        help="退回原版的全量重算口径（逐次比较更保守，但慢 3-8 倍）",
    )
    parser.add_argument(
        "--checkpoint",
        default=None,
        help="检查点文件路径：搜索中断后可从检查点续跑。默认关闭。仅在 --local-search > 0 时生效",
    )
    parser.add_argument(
        "--checkpoint-interval",
        type=_pos_float("检查点保存间隔"),
        default=60.0,
        help="检查点自动保存间隔（秒），默认 60 秒。仅 --checkpoint 生效时使用",
    )
    parser.add_argument(
        "--local-search-max-iters",
        type=_nonneg_int("局部搜索迭代上限"),
        default=0,
        help="局部搜索迭代次数上限（0=只受 --ls 时长约束）。"
        "设为正数可让固定 --seed 的局部搜索路径跨机器可复现",
    )
    parser.add_argument(
        "--near-optimal-top-k",
        dest="near_optimal_top_k",
        type=_pos_int("近优解收集容量"),
        default=DEFAULT_NEAR_OPTIMAL_TOP_K,
        help="近优解（near-optimal orders）收集器容量：局部搜索期间"
        f"为「稳健性/敏感性摘要」保留多少个**去重**近优排序，"
        f"也是摘要能看到多少个近优序的硬上限"
        f"（默认 {DEFAULT_NEAR_OPTIMAL_TOP_K}；须 >= 1）。"
        "仅在 --ls > 0 时有意义；调大它可让敏感性摘要考察更多"
        "近优拓扑（内存与 summary 计算量随 K·n^2 增长）。生效值"
        "记入 run_metadata['params'] 与 HTML 报告",
    )

    # ---- 两阶段流程的"约束生成"阶段输出----
    parser.add_argument(
        "-o",
        "--constraints-out",
        default=None,
        help="仅 --from <tool> / --from-ale 模式：把转换得到的统一"
        "约束集写出到该路径（5 列 ALE 逗号格式，可被文本解析器"
        "读回）后提前返回，不进入正式排序；供 Snakemake/Nextflow"
        " 两阶段流程使用。文本约束模式下无效（输入已是统一格式）",
    )

    # ---- 输出文件命名风格 ----
    parser.add_argument(
        "--output-style",
        choices=["short", "legacy"],
        default=DEFAULT_OUTPUT_STYLE,
        help="三输出文件命名风格：short=简洁（默认，如 "
        "<前缀>.mt.informative.tsv）/ legacy=与原版逐字节等价的长名",
    )

    # ---- 输出前缀 ----
    parser.add_argument(
        "-p",
        "--output-prefix",
        dest="output_prefix",
        default=None,
        help="输出文件前缀（默认取约束文件路径）；设置后三输出文件与 HTML 报告均写到该前缀路径下",
    )

    # ---- 覆盖保护----
    parser.add_argument(
        "-f",
        "--force",
        action="store_true",
        help="允许覆盖已存在的输出文件（默认拒绝覆盖并以退出码 3 中止，"
        "以免静默破坏上一次运行的结果）",
    )

    return parser


def main(argv: Optional[List[str]] = None) -> None:
    """CLI 主函数。

    退出码约定：

    * ``0`` —— 成功（含 ``--dry-run`` 预检通过）；
    * ``1`` —— ``--dry-run`` 预检存在 error 级问题（依据结构化 ``report["ok"]``）；
    * ``2`` —— argparse 取值域错误（由 argparse 自身返回）；
    * ``3`` —— 输出文件已存在且未加 ``--force``；
    * ``4`` —— 运行时异常（输入不满足算法前置条件等）。
    """
    args = build_parser().parse_args(argv)
    try:
        result = api.rank(
            species_tree_path=args.species_tree,
            constraints_path=args.constraints,
            seed=args.seed,
            local_search=args.local_search,
            temperature=args.temperature,
            random_type=args.random_type,
            min_transfer_distance=args.min_transfer_distance,
            threshold_constraints=args.threshold_constraints,
            random_trees=args.random_trees,
            output_prefix=args.output_prefix,
            from_ale=args.from_ale,
            from_tool=("auto" if args.from_auto else args.from_tool),
            artra_transfer_kind=args.artra_transfer_kind,
            target_clade=args.target_clade,
            ancestor_map=args.target_clade_ancestor_map,
            dry_run=args.dry_run,
            html_report=not args.no_html,
            ale_min_support=args.ale_min_support,
            ale_min_family_size=args.ale_min_family_size,
            ale_cache_dir=args.ale_cache_dir,
            ale_source=args.ale_source,
            ale_parallel=args.ale_parallel,
            adapter_min_endpoint_hit_rate=args.min_endpoint_hit_rate,
            adapter_quiet=args.quiet_adapters,
            mcmc=args.mcmc,
            mcmc_iters=args.mcmc_iters,
            mcmc_temperature=args.mcmc_temperature,
            mcmc_burn_in=args.mcmc_burn_in,
            mcmc_thin=args.mcmc_thin,
            incremental=args.incremental,
            checkpoint_path=args.checkpoint,
            checkpoint_interval=args.checkpoint_interval,
            constraints_out=args.constraints_out,
            output_style=args.output_style,
            force=args.force,
            local_search_max_iterations=args.local_search_max_iters,
            top_k=args.near_optimal_top_k,
        )
    except FileExistsError as exc:
        sys.stderr.write(f"{exc}\n")
        sys.exit(3)
    except (ValueError, KeyError, RecursionError) as exc:
        sys.stderr.write(f"错误：{exc}\n")
        sys.exit(4)

    # 预检（dry-run）退出码：只取结构化 ok 字段。
    # 旧实现用 `"错误" in report_str` 子串匹配，而结论行恒含"错误 N"，
    # 导致完全合法的输入也返回退出码 1。
    if args.dry_run:
        meta = getattr(result, "run_metadata", {}) or {}
        if meta.get("dry_run"):
            if not meta.get("dry_run_ok", False):
                sys.exit(1)


if __name__ == "__main__":
    main()
