"""Python API：``rank`` / ``build_constraints``。

本模块是 ``rank`` 的统一编排入口，也是各可选插入点（约束转换 / 目标类群剪裁 /
预检 / HTML 报告）的挂载处。所有可选特性**默认关闭，关闭时不影响核心排序结果**：

* ``from_ale`` / ``from_tool``：把上游工具输出经适配器转为约束；
* ``target_clade`` / ``ancestor_map``：目标类群剪裁；
* ``dry_run``：正式排序前的静态预检；预检结论同时以**结构化字段**
  ``Result.run_metadata["dry_run_ok"]`` 回传（CLI 退出码不得靠子串匹配）；
* ``html_report``：交互式 HTML 报告；
* ``force``：允许覆盖既有产物（，CLI ``--force``）。

本轮（第二遍）新增的三个可选项，**默认值都保持历史行为**：

* ``top_k``（CLI ``--near-optimal-top-k``）：近优解收集器容量 —— 稳健性/敏感性摘要
  能看过多少个去重近优排序（，默认 50）；
* ``stop_check``：协作式取消谓词，在局部搜索 / MCMC 的**迭代边界**真正打断计算
  （Studio 的「取消」按钮）；
* 压缩输入：``.gz`` 单文件与上游官方 ``.tgz`` 示例（如 ALE 的
  ``examples/reconciliations.tgz``）都可直接作为输入，无需先手工解压。
"""

import os
import sys
from typing import Callable, Dict, List, Optional, Tuple, Union

from maxtic_next.config import (
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
    DEFAULT_MCMC_ITERS,
)
from maxtic_next.constraints.constraint import ConstraintSet
from maxtic_next.ranking.ranker import Ranker, Result
from maxtic_next.random_ import RandomWrapper

#: 适配器端点命中率告警阈值的"未指定"哨兵：实际默认值由适配器层
#: (``adapters._diagnostics.DEFAULT_MIN_ENDPOINT_HIT_RATE``) 提供，延迟导入以保持
#: ``import maxtic_next.api`` 轻量（适配器及其可选依赖只在真正用到时加载）。
ADAPTER_HIT_RATE_UNSPECIFIED = None


def build_constraints(constraints_path: str) -> ConstraintSet:
    """读取并解析约束文件，返回 ``ConstraintSet``。

    注：按距离过滤由 ``ConstraintSet.filter_by_distance`` 在排序阶段执行，
    此处仅做解析。
    """
    from maxtic_next.io.parsing import read_constraints_file

    return read_constraints_file(constraints_path)


def _expand_archive_inputs(
    cons_list: List[str], tool: Optional[str]
) -> Tuple[List[str], Dict[str, List[str]]]:
    """把输入里的 **tar 归档**（上游官方 ``.tgz`` / ``.tar.gz`` 示例）展开为成员输入。

    ：MaxTiC 的 ``examples/reconciliations.tgz`` 里是 1000 个
    ``*.uml_rec``，旧实现要求用户先手工 ``tar -zxvf`` 且**不给任何提示**。这里只处理
    "一个归档里有多份上游输出"；**单文件** gzip（``.gz`` / gzip 流）由读取层
    （``io.compression.open_text``）透明解压，不经过本函数。

    Returns:
        ``(展开后的输入列表, {归档路径: 展开得到的输入列表})``；第二个返回值供 stderr
        提示与 ``run_metadata["archives_expanded"]`` 留痕（成员活在临时目录里）。
    """
    from maxtic_next.io.compression import expand_archive_inputs

    return expand_archive_inputs(cons_list, tool)


def _archive_input_notes(archive_expansions: Dict[str, List[str]]) -> List[str]:
    """为每个被展开的归档生成一条 stderr 提示（含临时解包位置，便于复核）。"""
    notes: List[str] = []
    for archive, inputs in archive_expansions.items():
        where = inputs[0] if len(inputs) == 1 else os.path.dirname(inputs[0])
        notes.append(
            f"NOTE: 输入 {archive} 是 tar 归档，已自动展开为 {len(inputs)} 项输入"
            f"（解包于 {where}）。三输出文件与 HTML 报告的前缀取第一项输入；如需为"
            "每个成员分别留存产物，请显式指定 -p/--output-prefix 或分别运行。"
            "单文件 gzip（.gz）无需解包，读取层已透明解压；含多个成员的归档若被当作"
            "**一个**文本输入读取，会给出含 tar 解包命令的提示而非 traceback。"
        )
    return notes


def rank(
    species_tree_path: str,
    constraints_path: Union[str, List[str]],
    seed: int = DEFAULT_SEED,
    local_search: float = DEFAULT_TIME_FOR_SEARCH,
    temperature: float = DEFAULT_TEMPERATURE,
    random_type: int = DEFAULT_RANDOM_TYPE,
    min_transfer_distance: float = DEFAULT_MIN_TRANSFER_DIST,
    threshold_constraints: float = DEFAULT_THRESHOLD_CONSTRAINTS,
    random_trees: int = DEFAULT_RANDOM_TREES,
    output_prefix: Optional[str] = None,
    print_summary: bool = True,
    *,
    from_ale: bool = False,
    from_tool: Optional[str] = None,
    target_clade: Optional[str] = None,
    ancestor_map: bool = False,
    dry_run: bool = False,
    html_report: bool = True,
    ale_min_support: float = 0.05,
    ale_min_family_size: int = 5,
    ale_cache_dir: Optional[str] = None,
    ale_source: str = "trf",
    ale_parallel: str = "process",
    artra_transfer_kind: str = "all",
    adapter_min_endpoint_hit_rate: Optional[float] = ADAPTER_HIT_RATE_UNSPECIFIED,
    adapter_quiet: bool = False,
    mcmc: bool = False,
    mcmc_iters: int = DEFAULT_MCMC_ITERS,
    mcmc_temperature: float = MCMC_TEMPERATURE_AUTO,
    mcmc_burn_in: Optional[int] = None,
    mcmc_thin: Optional[int] = None,
    incremental: bool = True,
    checkpoint_path: Optional[str] = None,
    checkpoint_interval: float = 60.0,
    constraints_out: Optional[str] = None,
    output_style: str = DEFAULT_OUTPUT_STYLE,
    force: bool = False,
    local_search_max_iterations: int = 0,
    top_k: int = DEFAULT_NEAR_OPTIMAL_TOP_K,
    stop_check: Optional[Callable[[], bool]] = None,
) -> Result:
    """运行 MaxTiC 排序流程，返回 ``Result``。

    Args:
        species_tree_path: 物种树（Newick）路径。
        constraints_path: 约束文件路径（文本双格式），或上游工具输出文件/目录列表。
        seed: 随机种子（驱动 mix 平局与局部搜索），默认 42。
        local_search: 局部搜索时长（秒），0 表示关闭。
        temperature: Metropolis 温度（``local_search > 0`` 时必须 > 0）。
        random_type: 随机化类型，只接受 ``0`` / ``1`` / ``2``。
        min_transfer_distance: 最小转移距离阈值（按 phylogenetic distance 列过滤）。
        threshold_constraints: 约束权重阈值比例，取值域 ``[0, 1]``。
        random_trees: 随机树采样数量，0 表示关闭。
        output_prefix: 输出文件前缀，默认使用**第一个**约束文件路径（多文件时会在
            ``run_metadata["constraint_files_merged"]`` 与 stdout 中显式记录）。
        print_summary: 是否打印 stdout 摘要。
        from_ale: 是否使用 ALE 适配器模式（``constraints_path`` 视为 ``.uml_rec`` 列表）。
            向后兼容开关，等价 ``from_tool="ale"``。
        from_tool: 上游工具适配器名（``"ale"`` / ``"ranger"`` / ``"eccetera"`` /
            ``"artra"`` / ``"alerax"`` / ``"auto"``）。设置后 ``constraints_path`` 视为
            该工具输出文件/目录列表，经注册表分发到对应适配器转为统一约束；
            ``"auto"`` 时按内容/扩展名自动检测。优先级高于 ``from_ale``。
        target_clade: 目标类群根节点标签；``None`` 表示不剪裁。
        ancestor_map: 是否开启"外部祖先映射到根"的保守剪裁规则（需配合 ``target_clade``）。
        dry_run: 仅做静态预检并打印报告，不运行正式排序。预检结论同时写入
            ``Result.run_metadata`` 的 ``dry_run_ok`` / ``dry_run_report``，供调用方
            以**结构化字段**判定退出码。适配器模式下 ``dry_run`` 不再只是
            "文件可读"检查：会真正解析上游输出并报告端点命中率 / 零约束问题。
        html_report: 是否生成交互式 HTML 报告（默认 True；``--no-html`` 对应 False）。
        ale_min_support: ALE 适配器单家族内最小支持度阈值。
        ale_min_family_size: ALE 适配器最小基因家族规模。
        ale_cache_dir: ALE 适配器文件级缓存目录（断点续传）。
        ale_source: ALE 约束来源 ``"trf"``（**默认**）/ ``"rec"``。默认取 ``trf``，
            即 ALE 官方 MaxTiC 集成（``maxtic/constraints_from_reconciliations.py``
            的 ``constraints_from_transfers``：供体的父节点须早于受体的子节点）所
            规定的口径；``rec``（按调和事件）仍可选，但两者产出的约束集
            **规模与权重量级都不同**，切换即换掉整个输入。
        adapter_min_endpoint_hit_rate: 适配器端点命中率告警阈值：
            低于该比例时向 stderr 告警，命中率为 0 时直接报错（不再产出"静默的
            空约束集"）。``None`` 表示沿用适配器层默认（0.5）。
        adapter_quiet: 抑制适配器在 stderr 上的例行摘要与告警。
            **只影响打印**：诊断量仍完整写入
            ``Result.run_metadata["adapter_diagnostics"]``；程序性错误（0 命中、
            上游 XML 解析失败）依然抛出。
        ale_parallel: ALE 并行解析模式 ``"process"``（默认）/ ``"thread"``。
        artra_transfer_kind: ARTra 转移类别。
        mcmc: 是否开启线性扩展上的 Metropolis–Hastings 采样（默认 False）。
        mcmc_iters: MH 链步数（必须 > 0）。
        mcmc_temperature: Boltzmann 温度；``0.0``（默认）表示 ``auto``，
            有效值 ``max(total_weight, 1) / 100``，会被打印并记入
            ``run_metadata["params"]``。
        mcmc_burn_in: burn-in 步数；``None`` 表示取 ``mcmc_iters`` 的一半。
        mcmc_thin: 抽稀间隔；``None`` 表示 ``max(1, mcmc_iters // 100)``。
        incremental: 是否启用增量 scoring。**默认开启**：
            增量与全量在逐字节等价性扫描（4 组输入 × 5 个种子 × 7 个产物字段全一致）
            下等价，由 tests/test_output_equivalence.py 的成对用例与 tests/test_delivery_guards.py 看守；
            传 False（CLI ``--no-incremental``）可退回全量重算口径。
        checkpoint_path: 检查点文件路径；``None``（默认）禁用。
        checkpoint_interval: 检查点自动保存间隔（秒），默认 60 秒。
        constraints_out: 仅适配器模式有效。写出统一约束集后提前返回；
            ``dry_run=True`` 时同样会写出该文件（不再静默忽略 ``-o``）。
        output_style: 三输出文件命名风格，``"short"``（默认）或 ``"legacy"``。
        force: 允许覆盖已存在的产物文件（对应 CLI ``--force``）。
        local_search_max_iterations: 局部搜索迭代上限（0 = 仅受时长约束）。
            用于把"固定 seed 可复现"的承诺扩展到 ``--ls > 0`` 路径。
        top_k: **近优解收集容量**：局部搜索期间为"稳健性/敏感性摘要"
            保留多少个**去重**近优排序；必须为 ``>= 1`` 的整数，默认 ``50``
            （``config.DEFAULT_NEAR_OPTIMAL_TOP_K``，与历史行为一致）。它同时是摘要能
            看到多少个近优序的**硬上限**，因此 ``local_search > 0`` 时才有效果，
            生效值记入 ``run_metadata["params"]`` 与
            ``run_metadata["near_optimal_top_k"]``。对应 CLI
            ``--near-optimal-top-k``。
        stop_check: **协作式取消回调**（``() -> bool``）：在局部搜索与
            MCMC 链的**迭代边界**打断计算，返回"截至取消点"的结果而不是抛异常。
            Studio（独立仓库）的「取消」按钮经其 ``cancellation.StopToken`` 注入该参数。
            默认 ``None``（不检查）；是否真的打断见 ``run_metadata["cancelled"]``。

    Returns:
        包含排序结果与输出文件路径的 ``Result``（dry_run 时含 ``dry_run_report``）。
    """
    from maxtic_next.io.parsing import read_constraints_file, read_newick_file

    # ---- 参数即校验：dry-run 路径也必须拒绝非法取值 ----
    Ranker._validate_params(
        min_transfer_distance=min_transfer_distance,
        threshold_constraints=threshold_constraints,
        temperature=temperature,
        random_type=random_type,
        time_for_search=local_search,
        random_trees=random_trees,
        mcmc_iters=mcmc_iters,
        mcmc_temperature=mcmc_temperature,
        mcmc_burn_in=mcmc_burn_in,
        mcmc_thin=mcmc_thin,
        local_search_max_iterations=local_search_max_iterations,
        near_optimal_top_k=top_k,
    )

    # 规范化约束路径为列表
    if isinstance(constraints_path, (list, tuple)):
        cons_list = [str(p) for p in constraints_path]
    else:
        cons_list = [str(constraints_path)]

    # ---- ：tar 归档输入（上游官方示例 .tgz）自动展开为成员文件 ----
    tree = read_newick_file(species_tree_path)

    # 适配器名（优先级：--from > --from-ale）
    tool = from_tool if from_tool else ("ale" if from_ale else None)

    # 归档展开需要知道工具名（AleRax 的输出树归档要按目录约定取文件），故排在
    # ``tool`` 判定之后、预检与解析之前：此后所有分支看到的都是"真正的输入文件"。
    cons_list, archive_expansions = _expand_archive_inputs(cons_list, tool)
    archive_notes = _archive_input_notes(archive_expansions)
    for archive_note in archive_notes:
        # 走 stderr：stdout 仍是可机读的原版口径摘要（与多文件前缀提示同一通道）
        sys.stderr.write(archive_note + "\n")
        sys.stderr.flush()
    if tool is None and constraints_out is not None:
        # -o/--constraints-out 仅在适配器模式下有意义（写"转换后的统一约束"）；
        # 文本模式下静默忽略会误导用户，显式报错。
        raise ValueError(
            "constraints_out/-o 仅在适配器模式（--from <tool> / --from-ale）下有效；"
            "文本约束模式请直接使用输入文件（已是统一格式）。"
        )

    # ---- 预检（可选插入点，早于正式排序）----
    adapter_diagnostics: Dict = {}
    if dry_run:
        from maxtic_next.dry_run import dry_run_check, format_dry_run

        report = dry_run_check(
            tree,
            cons_list,
            from_ale=(tool is not None),
            min_transfer_distance=min_transfer_distance,
            target_clade=target_clade,
        )
        # ：适配器模式的预检不再只做"文件可读"检查——真正解析上游输出，
        # 并把"零约束 / 端点命中率过低"如实报为问题（旧实现对垃圾输入判"通过"）。
        if tool is not None:
            issues = report["issues"]
            try:
                cset_dry = _convert_via_registry(
                    tool,
                    tree,
                    cons_list,
                    constraints_out=constraints_out,
                    ale_min_support=ale_min_support,
                    ale_min_family_size=ale_min_family_size,
                    ale_cache_dir=ale_cache_dir,
                    ale_parallel=ale_parallel,
                    ale_source=ale_source,
                    artra_transfer_kind=artra_transfer_kind,
                    min_endpoint_hit_rate=adapter_min_endpoint_hit_rate,
                    quiet=adapter_quiet,
                )
                resolved = cset_dry[0]
                cset = cset_dry[1]
                adapter_diagnostics = _adapter_diagnostics_of(cset)
            except Exception as exc:  # noqa: BLE001 - 预检须把异常转为报告项
                issues.append(
                    {
                        "severity": "error",
                        "message": f"上游输出（格式 {tool}）解析失败：{type(exc).__name__}: {exc}",
                    }
                )
                exc_diagnostics = getattr(exc, "diagnostics", None)
                if exc_diagnostics:
                    adapter_diagnostics = _as_dict(exc_diagnostics)
            else:
                valid_labels = set(str(x) for x in tree.internal_node_labels()) | set(
                    tree.get_leaves_names()
                )
                hit = sum(
                    1
                    for c in cset.constraints
                    if c.donor in valid_labels and c.receptor in valid_labels
                )
                if not cset.constraints:
                    issues.append(
                        {
                            "severity": "error",
                            "message": f"上游输出（格式 {resolved}）未产出任何约束："
                            "正式排序将在 0 总权重上产出一个与数据无关的"
                            "「平凡最优」排序（ 的实际后果）。"
                            "请检查工具输出格式、命名（类群名 vs 数字内部节点 ID）"
                            "与 --ale-min-support / --ale-min-family-size 阈值。",
                        }
                    )
                elif hit == 0:
                    issues.append(
                        {
                            "severity": "error",
                            "message": f"上游输出（格式 {resolved}）的 {len(cset.constraints)} 条约束"
                            "端点全部不在物种树中：正式排序将退化为无信息约束。"
                            "请核对物种树叶名与上游工具所用命名。",
                        }
                    )
                elif hit < len(cset.constraints):
                    issues.append(
                        {
                            "severity": "warning",
                            "message": f"上游输出（格式 {resolved}）的 {len(cset.constraints)} 条约束中"
                            f"仅 {hit} 条端点命中物种树标签"
                            f"（命中率 {hit * 100.0 / len(cset.constraints):.1f}%）。",
                        }
                    )
            # 计数与 ``ok`` 必须在**两条**路径上重算：若只在正常路径重算，
            # "适配器抛异常"那条路径追加了 error 项却不会更新 ``n_errors``/``ok``，
            # 报告就会一边打印 ``[error]``、一边判定"通过"并以退出码 0 结束。
            report["n_errors"] = sum(1 for it in issues if it["severity"] == "error")
            report["n_warnings"] = sum(1 for it in issues if it["severity"] == "warning")
            report["ok"] = report["n_errors"] == 0
        report_str = format_dry_run(report)
        if print_summary:
            print(report_str)
        return Result(
            constraint_file=cons_list[0] if cons_list else species_tree_path,
            # ：退出码取自这个结构化字段，而非对 report_str 的子串匹配
            run_metadata={
                "dry_run": True,
                "dry_run_ok": bool(report["ok"]),
                "dry_run_n_errors": int(report["n_errors"]),
                "dry_run_n_warnings": int(report["n_warnings"]),
                # ：适配器诊断随预检一并回传
                "adapter_diagnostics": adapter_diagnostics,
                # ：归档输入展开为哪些成员（可追溯）
                "archives_expanded": dict(archive_expansions),
                "archive_notes": list(archive_notes),
            },
            dry_run_report=report_str,
        )

    # ---- 构建约束集：上游工具适配器（注册表分发）或文本解析 ----
    constraint_files_merged = list(cons_list)
    multi_file_note = ""
    if tool is not None:
        resolved, cset = _convert_via_registry(
            tool,
            tree,
            cons_list,
            constraints_out=constraints_out,
            ale_min_support=ale_min_support,
            ale_min_family_size=ale_min_family_size,
            ale_cache_dir=ale_cache_dir,
            ale_parallel=ale_parallel,
            ale_source=ale_source,
            artra_transfer_kind=artra_transfer_kind,
            min_endpoint_hit_rate=adapter_min_endpoint_hit_rate,
            quiet=adapter_quiet,
        )
        adapter_diagnostics = _adapter_diagnostics_of(cset)
        adapter_diagnostics.setdefault("from_tool", resolved)
        # 两阶段流程的"约束生成"阶段：仅写出约束文件即返回，不进入正式排序
        if constraints_out is not None:
            if print_summary:
                print(f"[约束生成] {resolved} 约束已写出至：{constraints_out}")
            return Result(
                constraint_file=cons_list[0] if cons_list else species_tree_path,
                run_metadata={
                    "from_tool": resolved,
                    "constraints_out": constraints_out,
                    "adapter_diagnostics": adapter_diagnostics,
                    "archives_expanded": dict(archive_expansions),
                    "archive_notes": list(archive_notes),
                },
            )
        constraint_file = cons_list[0]
        if len(cons_list) > 1 and output_prefix is None:
            multi_file_note = _multi_file_prefix_note(cons_list)
    else:
        # 文本双格式：支持一个或多个约束文件（CLI synopsis "CONSTRAINTS [...]" 与
        # 手册 03 的约定），全部解析后合并为单一 ConstraintSet（边键聚合在
        # EdgeBuilder 中进行，与文件顺序无关）。
        if not cons_list:
            raise ValueError("未提供约束文件。")
        csets = [read_constraints_file(p) for p in cons_list]
        merged = ConstraintSet()
        for cs in csets:
            for c in cs:
                merged.add(c)
        cset = merged
        constraint_file = cons_list[0]
        if len(cons_list) > 1 and output_prefix is None:
            # ：产物前缀只体现第一个文件，必须显式告知而非静默
            multi_file_note = _multi_file_prefix_note(cons_list)
    if multi_file_note:
        # 走 stderr：stdout 只承载"与原版逐字节对齐"的摘要，命名口径这类
        # 运行说明不得污染可被机器解析的 stdout。
        sys.stderr.write(multi_file_note + "\n")
        sys.stderr.flush()

    # ---- 目标类群剪裁（可选插入点，约束解析后、建边前）----
    if target_clade is not None:
        from maxtic_next.prune import prune_constraints

        cset = prune_constraints(tree, cset, target_clade, ancestor_map=ancestor_map)

    rng = RandomWrapper(seed)
    ranker = Ranker(
        tree=tree,
        cset=cset,
        rng=rng,
        min_transfer_distance=min_transfer_distance,
        threshold_constraints=threshold_constraints,
        temperature=temperature,
        random_type=random_type,
        time_for_search=local_search,
        random_trees=random_trees,
        constraint_file=constraint_file,
        output_prefix=output_prefix,
        mcmc=mcmc,
        mcmc_iters=mcmc_iters,
        mcmc_temperature=mcmc_temperature,
        mcmc_burn_in=mcmc_burn_in,
        mcmc_thin=mcmc_thin,
        incremental=incremental,
        checkpoint_path=checkpoint_path,
        checkpoint_interval=checkpoint_interval,
        output_style=output_style,
        force=force,
        constraint_files_merged=constraint_files_merged,
        local_search_max_iterations=local_search_max_iterations,
    )
    result = ranker.run(
        print_summary=print_summary,
        html_report=html_report,
        seed=seed,
        top_k=top_k,
        stop_check=stop_check,
    )
    if multi_file_note:
        # ：stderr 已声明，同时结构化留痕，供 API 调用方自检
        result.run_metadata["multi_file_prefix_note"] = multi_file_note
    if archive_expansions:
        # ：归档输入被展开成哪些成员（成员活在临时目录里）必须可追溯
        result.run_metadata["archives_expanded"] = dict(archive_expansions)
        result.run_metadata["archive_notes"] = list(archive_notes)
    return result


def _convert_via_registry(
    tool: str,
    tree,
    cons_list: List[str],
    *,
    constraints_out: Optional[str],
    ale_min_support: float,
    ale_min_family_size: int,
    ale_cache_dir: Optional[str],
    ale_parallel: str,
    ale_source: str,
    artra_transfer_kind: str,
    min_endpoint_hit_rate: Optional[float] = None,
    quiet: bool = False,
):
    """经注册表把上游输出转为统一约束集，返回 ``(解析到的工具名, ConstraintSet)``。

    ：``auto`` 模式下只做**一次**格式嗅探（旧实现先 ``detect_formats`` 再
    ``convert_auto``，同一批文件被嗅探两遍，无扩展名的大目录代价翻倍）。

    ：``auto`` 遇到**混合格式**时不再只报错 —— 按各自工具的适配器分组
    解析后经 ``merge_constraint_sets`` 合并，同时把"供体端点层级不同 / 权重口径不可
    通约"这一事实打到 stderr 并记入 ``diagnostics``（``convention_conflict`` /
    ``semantics_conflict``），由调用方随 ``Result`` 上报，不再静默相加。
    """
    from maxtic_next.constraints.adapters import registry

    conv_kwargs = dict(
        min_support=ale_min_support,
        min_family_size=ale_min_family_size,
        cache_dir=ale_cache_dir,
        parallel=ale_parallel,
        output_path=constraints_out,
        ale_source=ale_source,
        artra_transfer_kind=artra_transfer_kind,
        min_endpoint_hit_rate=min_endpoint_hit_rate,
        quiet=quiet,
    )
    if tool == "auto":
        try:
            resolved = registry.detect_formats(cons_list)
        except registry.MixedFormatError as exc:
            # 嗅探结果已在异常里（分组），不再二次读取文件
            sets = [
                registry.convert(t, tree, paths, **conv_kwargs)
                for t, paths in sorted(exc.groups.items())
            ]
            merged = registry.merge_constraint_sets(sets, quiet=bool(quiet))
            return "auto+" + "+".join(sorted(exc.groups)), merged
        if resolved is None:
            raise ValueError(
                "无法自动识别上游工具格式；请用 --from <tool> 显式指定"
                f"（可用：{registry.REGISTRY.names()}）。"
            )
        return resolved, registry.convert(resolved, tree, cons_list, **conv_kwargs)
    return tool, registry.convert(tool, tree, cons_list, **conv_kwargs)


def _as_dict(value) -> Dict:
    """把可能是 ``diagnostics`` 字典 / 异常携带值的东西规整为 dict。"""
    return dict(value) if isinstance(value, dict) else {}


def _adapter_diagnostics_of(cset) -> Dict:
    """取 ``ConstraintSet.diagnostics`` 的**浅拷贝**（不得只留在内存里）。"""
    return _as_dict(getattr(cset, "diagnostics", None))


def _multi_file_prefix_note(cons_list: List[str]) -> str:
    """多约束文件时的"产物前缀只体现第一个文件"声明。"""
    return (
        f"NOTE: 输入了 {len(cons_list)} 个约束文件，"
        f"三输出文件与 HTML 报告均以第一个文件为前缀：{cons_list[0]}"
        f"（其余 {len(cons_list) - 1} 个文件已参与计算：{cons_list[1:]}）。"
        "如需为不同来源分别留存产物，请对每次运行显式指定 -p/--output-prefix。"
    )
