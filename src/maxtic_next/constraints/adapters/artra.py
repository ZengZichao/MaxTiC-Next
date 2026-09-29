"""ARTra 调和结果适配器。

ARTra（Additive and Replacing Transfer classifier）在 RANGER-DTL 风格的 DTL 调和
基础上，用规则启发式 + 机器学习把每个转移事件分类为**加性转移（Additive
Transfer）**或**替换转移（Replacing Transfer）**。其最终输出（``ARTra-程序/output.txt``
实证）与 RANGER-DTLx 调和报告**同形**，仅在事件名上多出 ``Replacing`` / ``Additive``
前缀::

    ------------ Reconciliation for Gene Tree 1 (rooted) -------------
    Species Tree:
    (...);
    Gene Tree:
    (...);
    Reconciliation:
    H117_0: Leaf Node
     = LCA[H119, H121]: Speciation, Mapping --> H122
     = LCA[H117, H14]: Replacing Transfer, Mapping --> H15, Recipient --> H125
     = LCA[H159, H17]: Additive Transfer, Mapping --> H17, Recipient --> H159
    ...
    The minimum reconciliation cost is: 65 (Duplications: 4, Transfers: 13, Losses: 18)
    Total number of optimal solutions: 24

因此本适配器与 ``RangerDTLAdapter`` 共用同一"调和报告"解析核心
（``_recon_report``）：转移事件的 ``Mapping --> <donor>`` 为**供体物种**、
``Recipient --> <receptor>`` 为**受体物种**，产出统一 ``(donor, receptor)`` 有向约束。

语义说明：MaxTiC 只关心"donor 应早于/不晚于 receptor"的时间约束，
与转移是加性还是替换无关，故本适配器**默认统计所有转移事件**；如需仅取某一类别，
可用 ``transfer_kind`` 参数（``"all"`` / ``"replacing"`` / ``"additive"``）。

口径与可观测性
--------------------------------------------
* **供体端点约定** ``donor_endpoint_convention = "donor_itself"``（供体本身；ALE 取
  "供体的父" ``parent_of_donor``，两者相差一个物种层级，混用会被告警）；
* **支持度分母**＝报告**自己声明**的 ``Total number of optimal solutions: N``
  （真实 ``output.txt:206`` 声明 24 个最优解却只打印 1 个 → 分母是 24，不是 1）。
  旧实现把"块数"当分母，于是 13 条约束权重全为 1.0、``--min-support 0.95`` 成
  为空操作。声明值缺失时按诊断口径退化为**整数计数**并明说；
* 每一次静默丢弃（标签不匹配 / 支持度 / 家族规模 / 自环）都被计数，随
  ``ConstraintSet.diagnostics`` 返回并打印一行 stderr 摘要。

接口与其它适配器一致（``species_tree`` 在 ``convert`` 时传入），支持
``ProcessPoolExecutor`` 并行解析和文件级 pickle 缓存（断点续传）。
"""

from typing import Dict, List, Optional, Set, Tuple

from maxtic_next.constraints.constraint import Constraint, ConstraintSet
from maxtic_next.constraints.adapters._diagnostics import (
    CONVENTION_DONOR_ITSELF,
    DEFAULT_MIN_ENDPOINT_HIT_RATE,
    finalize_adapter_diagnostics,
    map_payloads,
    new_diagnostics,
    passes_support,
    reciprocal_subtract,
    weight_and_support,
)
from maxtic_next.constraints.adapters._recon_report import (
    parse_reconciliation_report_detailed,
    process_report_file,
    resolve_family_size_filter,
)
from maxtic_next.constraints.adapters._species import (
    build_valid_labels as _build_valid_labels,
    distance_from as _distance_from,
    species_cache_identity as _species_cache_identity,
)
from maxtic_next.io.compression import logical_stem, read_text_lines
from maxtic_next.tree.tree import Tree


# 与 ALE 一致的默认阈值
MINIMUM_SUPPORT_WITHIN_A_FAMILY = 0.05
MINIMUM_FAMILY_SIZE = 5


# ----------------------------------------------------------------------
# 模块级 worker（可 pickle，供 ProcessPoolExecutor 子进程调用）
# 物种树辅助（build_valid_labels / distance_from）见 ``_species`` 共享模块。
# 调和报告解析核心见 ``_recon_report``（与 RANGER-DTLx 共用）。
# ----------------------------------------------------------------------


def _process_artra_worker(args: Tuple) -> Tuple[Dict[str, int], int, int, str, Dict]:
    """模块级 worker（可 pickle）：解析单个 ARTra ``output.txt`` 文件。

    委托共享 ``_recon_report.process_report_file``（ARTra 与 RANGER-DTLx 同形，
    共用同一解析核心，保证行为一致）。

    Args:
        args: ``(file_path, min_family_size, cache_dir, min_support,
               valid_labels[, species_identity])``。

    Returns:
        ``(transfers, blocks, family_size, family, stats)`` 五元组。
    """
    return process_report_file(args)


class ARTraAdapter:
    """ARTra 调和结果 -> MaxTiC-Next 统一 ``ConstraintSet``。

    与 ``RangerDTLAdapter`` 接口一致（``species_tree`` 在 ``convert`` 时传入），
    支持 ``ProcessPoolExecutor`` 并行解析和文件级 pickle 缓存。
    """

    def __init__(
        self,
        min_support: float = MINIMUM_SUPPORT_WITHIN_A_FAMILY,
        min_family_size: int = MINIMUM_FAMILY_SIZE,
        cache_dir: Optional[str] = None,
        transfer_kind: str = "all",
        min_endpoint_hit_rate: float = DEFAULT_MIN_ENDPOINT_HIT_RATE,
        quiet: bool = False,
    ) -> None:
        """初始化适配器。

        Args:
            min_support: 单基因家族内最小支持度阈值（默认 0.05，仅在权重口径为
                [0,1] 支持度时有意义；口径由报告声明的样本数决定）。
            min_family_size: 参与计算的最小基因家族规模（叶子数，默认 5；
                探测失败时跳过过滤并告警）。
            cache_dir: 基因家族文件级缓存目录；``None`` 表示不缓存。缓存键含
                物种树拓扑指纹。
            transfer_kind: 统计的转移类别，``"all"``（默认，加性+替换）/
                ``"replacing"``（仅替换转移）/ ``"additive"``（仅加性转移）。
            min_endpoint_hit_rate: 端点命中率告警阈值（默认 0.5；0 命中报错）。
            quiet: 不打印 stderr 摘要。
        """
        if transfer_kind not in ("all", "replacing", "additive"):
            raise ValueError(
                f"transfer_kind 必须为 'all' / 'replacing' / 'additive'，收到 {transfer_kind!r}"
            )
        self.min_support = min_support
        self.min_family_size = min_family_size
        self.cache_dir = cache_dir
        self.transfer_kind = transfer_kind
        self.min_endpoint_hit_rate = min_endpoint_hit_rate
        self.quiet = quiet

    @staticmethod
    def _reciprocal_subtract(counts: Dict[str, int]) -> None:
        """互反减法（委托共用实现，含自环不抵消的修正，见 ``_diagnostics``）。"""
        reciprocal_subtract(counts)

    def _process_family_kind(
        self, file_path: str, valid_labels: Set[str]
    ) -> Tuple[Dict[str, int], int, int, str, Dict]:
        """按 ``transfer_kind`` 过滤类别的解析（仅 ``transfer_kind != "all"`` 时使用）。

        ``"all"`` 时直接走共享 ``process_report_file``（含缓存）；否则本方法在主进程
        内按类别（Replacing / Additive）逐行过滤后再统计，不使用文件级缓存
        （类别过滤属少数场景，避免污染共享缓存键）。诊断与计数路径与 ``"all"``
        完全一致。
        """
        # ：共享读取器，.gz / 单成员 .tar.gz 上游输出可直接解析
        raw_lines = read_text_lines(file_path, encoding="utf-8")
        # 按类别保留转移行：非目标类别的转移行剔除（其计数进 stats.malformed_lines）
        want = "Replacing" if self.transfer_kind == "replacing" else "Additive"
        kept: List[str] = []
        dropped_by_kind = 0
        for line in raw_lines:
            if "Recipient" in line and "Transfer" in line and want not in line:
                dropped_by_kind += 1
                continue
            kept.append(line)
        transfers, blocks, family_size, stats = parse_reconciliation_report_detailed(
            kept, valid_labels
        )
        stats["transfer_kind_filtered"] = dropped_by_kind
        stats["files_seen"] = 1
        family = logical_stem(file_path)
        if resolve_family_size_filter(stats, family_size, self.min_family_size):
            stats["families_skipped"] = 1
            transfers = {}
        return transfers, blocks, family_size, family, stats

    def _merge_payload(
        self,
        cset: ConstraintSet,
        species_tree: Tree,
        payload: Tuple,
    ) -> None:
        """把一个家族的解析结果合并进 ``cset``（互反减法 + 支持度过滤 + 距离标注）。"""
        if len(payload) >= 5 and isinstance(payload[4], dict):
            transfers, blocks, family_size, family, stats = payload[:5]
        else:
            transfers, blocks, family_size, family = payload[:4]
            stats = new_diagnostics("artra", CONVENTION_DONOR_ITSELF)
            stats["blocks_seen"] = blocks
            stats["sample_denominator"] = blocks
        stats["files_seen"] = 1
        cset.diagnostics.setdefault("files", []).append(stats)
        self._reciprocal_subtract(transfers)
        for key, cnt in transfers.items():
            if cnt <= 0:
                stats["dropped_by_support"] = int(stats.get("dropped_by_support", 0)) + 1
                continue
            weight, support = weight_and_support(stats, cnt)
            if not passes_support(stats, weight, cnt, self.min_support):
                stats["dropped_by_support"] = int(stats.get("dropped_by_support", 0)) + 1
                continue
            donor, receptor = key.split(",")
            # ：自环事件（donor == receptor）的距离置 None 而非 0.0，
            # 否则会被 --min-transfer-distance 静默滤掉，无法进入 "to itself" 统计。
            distance = None if donor == receptor else _distance_from(species_tree, receptor, donor)
            cset.add(
                Constraint(
                    donor=donor,
                    receptor=receptor,
                    weight=float(weight),
                    metadata={
                        "family": family,
                        "support": support,
                        "distance": distance,  # None = 未知
                        "weight_semantics": stats.get("weight_semantics"),
                        "sample_denominator": stats.get("sample_denominator"),
                    },
                )
            )
            stats["constraints_kept"] = int(stats.get("constraints_kept", 0)) + 1

    def convert(
        self,
        species_tree: Tree,
        artra_files: List[str],
        output_path: Optional[str] = None,
        max_workers: Optional[int] = None,
        parallel: str = "process",
    ) -> ConstraintSet:
        """把 ARTra ``output.txt`` 文件列表转换为 ``ConstraintSet``。

        各文件解析相互独立，默认用进程池并行加速（``ProcessPoolExecutor``）；
        可选 ``parallel="thread"`` 回退到 ``ThreadPoolExecutor``。结果按文件顺序
        合并，由于 ``ConstraintSet.add`` 按边键聚合，合并顺序不影响最终结果。

        Args:
            species_tree: 物种树。
            artra_files: ARTra 输出文件路径列表。
            output_path: 可选，若提供则把约束以 CSV 格式写出。
            max_workers: 并行进程/线程数；``None`` 时取默认。
            parallel: 并行模式，``"process"``（默认）或 ``"thread"``。

        Returns:
            统一 ``ConstraintSet``，``diagnostics`` 记录权重口径、分母与全部
            丢弃计数。
        """
        if not artra_files:
            cset = ConstraintSet(diagnostics=new_diagnostics("artra", CONVENTION_DONOR_ITSELF))
            finalize_adapter_diagnostics(
                cset, self.min_endpoint_hit_rate, set(), "artra", quiet=True, files=0
            )
            if output_path is not None:
                self._write_constraints_file(cset, output_path)
            return cset

        valid_labels = _build_valid_labels(species_tree)
        identity = _species_cache_identity(species_tree)

        cset = ConstraintSet(diagnostics=new_diagnostics("artra", CONVENTION_DONOR_ITSELF))

        # transfer_kind != "all" 时走主进程类别过滤（不并行、不缓存）
        if self.transfer_kind != "all":
            for rf in artra_files:
                payload = self._process_family_kind(rf, valid_labels)
                self._merge_payload(cset, species_tree, payload)
        else:
            worker_args = [
                (rf, self.min_family_size, self.cache_dir, self.min_support, valid_labels, identity)
                for rf in artra_files
            ]
            payloads = map_payloads(
                _process_artra_worker, worker_args, max_workers=max_workers, parallel=parallel
            )
            for payload in payloads:
                if payload is not None:
                    self._merge_payload(cset, species_tree, payload)

        # 先落盘再跑守卫（守卫在端点 0 命中时抛错；-o 产物仍须可核对）
        if output_path is not None:
            self._write_constraints_file(cset, output_path)
        finalize_adapter_diagnostics(
            cset,
            self.min_endpoint_hit_rate,
            valid_labels,
            "artra",
            quiet=self.quiet,
            files=len(artra_files),
        )
        return cset

    @staticmethod
    def _write_constraints_file(cset: ConstraintSet, path: str) -> None:
        """把 ``ConstraintSet`` 以逗号格式（family,donor,receptor,weight,distance）写出。"""
        with open(path, "w", encoding="utf-8") as fh:
            fh.write("#family,donor,receptor,weight,distance\n")
            for c in cset.constraints:
                dist = c.metadata.get("distance")
                dist_str = "" if dist is None else str(dist)
                fam = c.metadata.get("family") or ""
                fh.write(f"{fam},{c.donor},{c.receptor},{c.weight},{dist_str}\n")


def convert_from_artra(
    species_tree: Tree,
    artra_files: List[str],
    min_support: float = MINIMUM_SUPPORT_WITHIN_A_FAMILY,
    min_family_size: int = MINIMUM_FAMILY_SIZE,
    cache_dir: Optional[str] = None,
    transfer_kind: str = "all",
    output_path: Optional[str] = None,
    max_workers: Optional[int] = None,
    parallel: str = "process",
    min_endpoint_hit_rate: float = DEFAULT_MIN_ENDPOINT_HIT_RATE,
    quiet: bool = False,
) -> ConstraintSet:
    """便捷函数：把 ARTra 输出文件列表转换为 ``ConstraintSet``。

    Args:
        species_tree: 物种树。
        artra_files: ARTra ``output.txt`` 文件路径列表。
        min_support: 单家族内最小支持度阈值。
        min_family_size: 最小基因家族规模。
        cache_dir: 文件级缓存目录（断点续传）。
        transfer_kind: 转移类别 ``"all"`` / ``"replacing"`` / ``"additive"``。
        output_path: 可选约束写出路径。
        max_workers: 并行解析进程/线程数。
        parallel: 并行模式，``"process"``（默认）或 ``"thread"``。
        min_endpoint_hit_rate: 端点命中率告警阈值（默认 0.5；0 命中必报错）。
        quiet: 不打印 stderr 摘要。

    Returns:
        统一 ``ConstraintSet``。
    """
    adapter = ARTraAdapter(
        min_support=min_support,
        min_family_size=min_family_size,
        cache_dir=cache_dir,
        transfer_kind=transfer_kind,
        min_endpoint_hit_rate=min_endpoint_hit_rate,
        quiet=quiet,
    )
    return adapter.convert(
        species_tree,
        artra_files,
        output_path=output_path,
        max_workers=max_workers,
        parallel=parallel,
    )
