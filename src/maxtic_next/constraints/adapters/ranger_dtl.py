"""Ranger-DTL 调和结果适配器。

Ranger-DTL（来自 Ranger 软件包）输出 DTL 调和结果，其转移事件直接给出
``(donor species, recipient species)`` 对，与 ALE 的 ``T@donor->receptor`` 语义
相似但格式不同。

本适配器把 Ranger-DTL 输出转换为统一 ``ConstraintSet``，核心逻辑：

* 解析 RANGER-DTLx 调和报告（``m<idx> = LCA[..]: Transfer, Mapping --> D, ...,
  Recipient --> R``），抽取每个转移事件的 ``(donor, receptor)`` 物种对（donor 取
  ``Mapping -->``，receptor 取 ``Recipient -->``）；委托 ``_recon_report`` 共享解析核心；
* 按基因家族（一个文件 = 一个家族）聚合转移事件计数，做互反减法
  （reciprocal subtraction），再按 ``min_support`` 过滤；
* 按 ``min_family_size``（最小基因家族规模，叶子数）跳过过小的家族；
* 产出统一 ``ConstraintSet``，每条约束 ``metadata`` 含 ``family`` / ``support`` /
  ``distance``（按物种树拓扑距离，等价原版把物种树所有枝长设为 1 后的距离）。

接口与 ``ALEAdapter`` 一致（但 ``species_tree`` 在 ``convert`` 时传入而非构造器），
支持 ``ProcessPoolExecutor`` 并行解析和文件级 pickle 缓存（断点续传）。

口径与可观测性
---------------------------------------------------------
* **供体端点约定** ``donor_endpoint_convention = "donor_itself"``：本适配器取
  ``Mapping -->`` 给出的**供体本身**为 donor；ALE 路径取"供体的父"
  （``parent_of_donor``）。同一转移事件在两条路径下相差**一个物种层级**，
  端点不可直接互换；混用时 ``registry.merge_constraint_sets`` / ``--from-auto``
  会给出警告。原生约定保持不变，本适配器不做层级换算。
* **权重口径**：分母优先取报告**自己声明**的样本数
  （``Total number of optimal solutions: N``，求证自 ``DTL-algorithm.h:1785/:2224``
  与 ARTra ``output.txt:206``）；块数仅在各块叶子集合一致时作近似分母；
  否则权重退化为**整数计数**并在诊断/stderr 中声明（见 ``_recon_report``）。
* **每次静默丢弃都被计数**：标签不匹配、支持度不足、家族规模不足/探测失败、
  自环（donor == receptor）；汇总在 ``ConstraintSet.diagnostics``，并在
  ``convert`` 结束时打印一行 stderr 摘要。端点命中率低于
  ``min_endpoint_hit_rate``（默认 50%）告警、为 0 时报错。

RANGER-DTLx 真实输出格式示例（求证自 ``Ranger-DTLx/DTL-algorithm.h:1840``）::

    Reconciliation for Gene Tree 1:
    H117: Leaf Node
     = LCA[H117, H121]: Speciation, Mapping --> H123
     = LCA[H117, H14]: Transfer, Mapping --> H15, Edge, Parent = P, Recipient --> H125
     = LCA[H18, H21]: Duplication, Mapping --> H22

其中转移事件的 ``Mapping -->`` 为供体物种、``Recipient -->`` 为受体物种。
上述格式已按 RANGER-DTLx 源码与真实输出逐项求证。
"""

from typing import Dict, List, Optional, Tuple

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
    process_report_file,
)
from maxtic_next.constraints.adapters._species import (
    build_valid_labels as _build_valid_labels,
    distance_from as _distance_from,
    species_cache_identity as _species_cache_identity,
)
from maxtic_next.tree.tree import Tree


# 与 ALE 一致的默认阈值
MINIMUM_SUPPORT_WITHIN_A_FAMILY = 0.05
MINIMUM_FAMILY_SIZE = 5


# ----------------------------------------------------------------------
# 模块级 worker（可 pickle，供 ProcessPoolExecutor 子进程调用）
# 物种树辅助（build_valid_labels / distance_from）见 ``_species`` 共享模块。
# ----------------------------------------------------------------------


def _process_ranger_dtl_worker(args: Tuple) -> Tuple[Dict[str, int], int, int, str, Dict]:
    """模块级 worker（可 pickle）：解析单个 RANGER-DTLx 调和报告文件。

    委托共享 ``_recon_report.process_report_file``（RANGER-DTLx 与 ARTra 同形，
    共用同一解析核心，保证行为一致）。

    Args:
        args: ``(file_path, min_family_size, cache_dir, min_support,
               valid_labels[, species_identity])``。

    Returns:
        ``(transfers, blocks, family_size, family, stats)`` 五元组。
    """
    return process_report_file(args)


# ----------------------------------------------------------------------
# 适配器类
# ----------------------------------------------------------------------


class RangerDTLAdapter:
    """Ranger-DTL 调和结果 -> MaxTiC-Next 统一 ``ConstraintSet``。

    与 ``ALEAdapter`` 接口一致（但 ``species_tree`` 在 ``convert`` 时传入），
    支持 ``ProcessPoolExecutor`` 并行解析和文件级 pickle 缓存。
    """

    def __init__(
        self,
        min_support: float = MINIMUM_SUPPORT_WITHIN_A_FAMILY,
        min_family_size: int = MINIMUM_FAMILY_SIZE,
        cache_dir: Optional[str] = None,
        source: str = "transfer",
        min_endpoint_hit_rate: float = DEFAULT_MIN_ENDPOINT_HIT_RATE,
        quiet: bool = False,
    ) -> None:
        """初始化适配器。

        Args:
            min_support: 单基因家族内最小支持度阈值（默认 0.05）。仅当权重口径为
                [0,1] 支持度时才有该语义；整数计数口径下退化为"计数 > 0"，
                该事实会写进诊断与 stderr。
            min_family_size: 参与计算的最小基因家族规模（叶子数，默认 5）。
                家族规模**探测失败**时跳过过滤并告警，而不是清空全部转移
                。
            cache_dir: 基因家族文件级缓存目录；``None`` 表示不缓存。缓存键含
                物种树**拓扑**指纹与版本号。
            source: 约束来源（接口一致性，始终使用转移事件）。
            min_endpoint_hit_rate: 端点命中率低于该值时告警（默认 0.5）；
                命中率为 0 直接报错。``<= 0`` 关闭守卫。
            quiet: 不打印 stderr 摘要（诊断仍随结果返回）。
        """
        self.min_support = min_support
        self.min_family_size = min_family_size
        self.cache_dir = cache_dir
        self.source = source
        self.min_endpoint_hit_rate = min_endpoint_hit_rate
        self.quiet = quiet

    @staticmethod
    def _reciprocal_subtract(counts: Dict[str, int]) -> None:
        """互反减法（委托共用实现，含自环不抵消的修正，见 ``_diagnostics``）。"""
        reciprocal_subtract(counts)

    def _merge_payload(
        self,
        cset: ConstraintSet,
        species_tree: Tree,
        payload: Tuple,
    ) -> None:
        """把一个家族的解析结果合并进 ``cset``。

        权重口径取自该文件的 ``stats``（声明样本数 / 一致块数 / 整数计数），
        互反抵消与支持度不足都在诊断中计数。

        Args:
            cset: 待累加的约束集。
            species_tree: 物种树（用于距离计算）。
            payload: ``(transfers, blocks, family_size, family, stats)``。
        """
        if len(payload) >= 5 and isinstance(payload[4], dict):
            transfers, blocks, family_size, family, stats = payload[:5]
        else:
            transfers, blocks, family_size, family = payload[:4]
            stats = new_diagnostics("ranger", CONVENTION_DONOR_ITSELF)
            stats["blocks_seen"] = blocks
            stats["sample_denominator"] = blocks
            stats["family_size_probe_ok"] = family_size > 0
        stats["files_seen"] = 1
        cset.diagnostics.setdefault("files", []).append(stats)
        self._reciprocal_subtract(transfers)
        for key, cnt in transfers.items():
            if cnt <= 0:
                # 互反对完全抵消（原版口径），单独计数不静默
                stats["dropped_by_support"] = int(stats.get("dropped_by_support", 0)) + 1
                continue
            weight, support = weight_and_support(stats, cnt)
            if not passes_support(stats, weight, cnt, self.min_support):
                stats["dropped_by_support"] = int(stats.get("dropped_by_support", 0)) + 1
                continue
            donor, receptor = key.split(",")
            # ：donor == receptor 的自环事件拓扑距离恒为 0；写成 0.0 会被
            # --min-transfer-distance 当成"真距离 0"滤掉，使该事件从原版口径的
            # "to itself" 统计中消失（文本输入路径则会保留）。置 None = 无距离信息，
            # 与文本路径同权处理，过滤与统计都如实进行。
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
        dtl_files: List[str],
        output_path: Optional[str] = None,
        max_workers: Optional[int] = None,
        parallel: str = "process",
    ) -> ConstraintSet:
        """把 Ranger-DTL 文件列表转换为 ``ConstraintSet``。

        各文件解析相互独立，默认用进程池并行加速
        （``ProcessPoolExecutor``）；可选 ``parallel="thread"`` 回退到
        ``ThreadPoolExecutor``。结果按文件顺序合并，由于
        ``ConstraintSet.add`` 按边键聚合，合并顺序不影响最终结果。

        Args:
            species_tree: 物种树。
            dtl_files: Ranger-DTL 输出文件路径列表。
            output_path: 可选，若提供则把约束以 CSV 格式写出。
            max_workers: 并行进程/线程数；``None`` 时取默认。
            parallel: 并行模式，``"process"``（默认）或 ``"thread"``。

        Returns:
            统一 ``ConstraintSet``：每条约束 ``metadata`` 含 ``family`` /
            ``support`` / ``distance`` / ``weight_semantics``；集合的
            ``diagnostics`` 携带口径与全部丢弃计数。
        """
        if not dtl_files:
            cset = ConstraintSet(diagnostics=new_diagnostics("ranger", CONVENTION_DONOR_ITSELF))
            finalize_adapter_diagnostics(
                cset, self.min_endpoint_hit_rate, set(), "ranger", quiet=True, files=0
            )
            if output_path is not None:
                self._write_constraints_file(cset, output_path)
            return cset

        valid_labels = _build_valid_labels(species_tree)
        identity = _species_cache_identity(species_tree)
        worker_args = [
            (rf, self.min_family_size, self.cache_dir, self.min_support, valid_labels, identity)
            for rf in dtl_files
        ]

        payloads = map_payloads(
            _process_ranger_dtl_worker, worker_args, max_workers=max_workers, parallel=parallel
        )

        cset = ConstraintSet(diagnostics=new_diagnostics("ranger", CONVENTION_DONOR_ITSELF))
        for payload in payloads:
            if payload is not None:
                self._merge_payload(cset, species_tree, payload)
        # 先写出已产出的约束，再跑守卫：守卫可能抛错（0 命中），但两阶段流程
        # 的产物（-o）必须仍然落盘，用户才能拿它去核对命名。
        if output_path is not None:
            self._write_constraints_file(cset, output_path)
        finalize_adapter_diagnostics(
            cset,
            self.min_endpoint_hit_rate,
            valid_labels,
            "ranger",
            quiet=self.quiet,
            files=len(dtl_files),
        )
        return cset

    @staticmethod
    def _write_constraints_file(cset: ConstraintSet, path: str) -> None:
        """把 ``ConstraintSet`` 以 CSV 格式写出。"""
        with open(path, "w", encoding="utf-8") as fh:
            fh.write("#family,donor,receptor,weight,distance\n")
            for c in cset.constraints:
                dist = c.metadata.get("distance")
                dist_str = "" if dist is None else str(dist)
                fam = c.metadata.get("family") or ""
                fh.write(f"{fam},{c.donor},{c.receptor},{c.weight},{dist_str}\n")


def convert_from_ranger_dtl(
    species_tree: Tree,
    dtl_files: List[str],
    min_support: float = MINIMUM_SUPPORT_WITHIN_A_FAMILY,
    min_family_size: int = MINIMUM_FAMILY_SIZE,
    cache_dir: Optional[str] = None,
    source: str = "transfer",
    output_path: Optional[str] = None,
    max_workers: Optional[int] = None,
    parallel: str = "process",
    min_endpoint_hit_rate: float = DEFAULT_MIN_ENDPOINT_HIT_RATE,
    quiet: bool = False,
) -> ConstraintSet:
    """便捷函数：把 Ranger-DTL 文件列表转换为 ``ConstraintSet``。

    Args:
        species_tree: 物种树。
        dtl_files: Ranger-DTL 输出文件路径列表。
        min_support: 单家族内最小支持度阈值。
        min_family_size: 最小基因家族规模。
        cache_dir: 文件级缓存目录（断点续传）。
        source: 约束来源（接口一致性）。
        output_path: 可选约束写出路径。
        max_workers: 并行解析进程/线程数。
        parallel: 并行模式，``"process"``（默认）或 ``"thread"``。
        min_endpoint_hit_rate: 端点命中率告警阈值（默认 0.5；0 命中必报错）。
        quiet: 不打印 stderr 摘要。

    Returns:
        统一 ``ConstraintSet``。
    """
    adapter = RangerDTLAdapter(
        min_support=min_support,
        min_family_size=min_family_size,
        cache_dir=cache_dir,
        source=source,
        min_endpoint_hit_rate=min_endpoint_hit_rate,
        quiet=quiet,
    )
    return adapter.convert(
        species_tree,
        dtl_files,
        output_path=output_path,
        max_workers=max_workers,
        parallel=parallel,
    )
