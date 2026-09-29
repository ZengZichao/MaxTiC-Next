"""AleRax 调和结果适配器。

AleRax 在贝叶斯框架下对基因家族做 DTL 调和与物种树打分，其输出目录内每个基因家族
在 ``reconciliations/summaries/`` 下的转移频率文件（AleRax 早期版本为
``<family>_transfers.txt``，现行版本为 ``<family>_meanTransfers.txt``，求证自 ``AleOptimizer.cpp`` 写 summaries 处）
汇总了**物种间转移频率**。经 ``AleRax-代码/scripts/extract_families_transfer.py``
求证，该文件每行形如::

    <donor_species> <receptor_species> <frequency>

其中 ``frequency`` 为该转移在后验样本中的期望频率（浮点，通常 0–1 或期望计数）。
与 RANGER/ARTra/ALE 的"计数 / 采样数"不同，AleRax 已在家族内**预聚合为频率**，
因此本适配器**直接以频率为权重**（不再除以块数）。

输入既可是若干 ``*_transfers.txt`` 文件，也可是一个 AleRax 输出目录或其
``reconciliations/summaries`` 目录——后者会被自动展开为其中的全部
``*_transfers.txt``（ 自动检测友好）。

约束语义：``donor -> receptor`` 表示 donor 应早于/不晚于 receptor。
端点合法性以物种树标签集合校验。每条约束 ``metadata`` 含
``family`` / ``support``（= 频率）/ ``distance``（物种树拓扑距离）。

说明：``_transfers.txt`` 无基因树叶子信息，故本适配器**不做 min_family_size 过滤**
（AleRax 已在家族内聚合）；``min_support`` 仍按频率阈值生效。

口径与可观测性：本适配器权重口径为
``weight_semantics = "posterior_frequency"``（上游已预聚合，**不再除以任何分母**，
故它与 ALE 的"count/采样数"支持度量纲不同，不可与后者直接相加比较）；供体端点约定为
``donor_itself``（表中第一列即供体本身；ALE 路径取"供体的父"，相差一个物种层级）。
每个被跳过的行（字段不足 / 频率非数字 / 端点不在物种树）都计入
``ConstraintSet.diagnostics`` 并在 ``convert`` 末尾打印一行 stderr 摘要；端点命中率为
0 时抛 ``LabelMismatchError``，不再静默产出空约束集。
"""

import os
from typing import Dict, List, Optional, Tuple

from maxtic_next.constraints.constraint import Constraint, ConstraintSet
from maxtic_next.constraints.adapters._diagnostics import (
    CONVENTION_DONOR_ITSELF,
    DEFAULT_MIN_ENDPOINT_HIT_RATE,
    WEIGHT_POSTERIOR_FREQUENCY,
    finalize_adapter_diagnostics,
    map_payloads,
    new_diagnostics,
    passes_support,
    reciprocal_subtract,
    weight_and_support,
)
from maxtic_next.constraints.adapters._species import (
    build_valid_labels as _build_valid_labels,
    distance_from as _distance_from,
    resolve_species_label as _resolve_species_label,
)
from maxtic_next.io.compression import logical_path, open_text
from maxtic_next.tree.tree import Tree


# 与其它适配器一致的默认阈值（min_family_size 对 AleRax 不生效，仅保接口一致）
MINIMUM_SUPPORT_WITHIN_A_FAMILY = 0.05
MINIMUM_FAMILY_SIZE = 5

# 兼容两代 AleRax 输出命名：早期 summaries 为 ``<family>_transfers.txt``
# （scripts/extract_families_transfer.py 约定）；现行 AleOptimizer.cpp 写
# ``<family>_meanTransfers.txt``（reconciliations/summaries/ 下，另有逐样本
# ``reconciliations/all/<family>_transfers_<sample>.txt``，此处不消费逐样本文件）。
_TRANSFER_SUFFIXES = ("_transfers.txt", "_meanTransfers.txt")
# 小写化后缀（用于大小写不敏感的匹配；如 "FAM1_Transfers.txt"）
_TRANSFER_SUFFIXES_LOWER = tuple(s.lower() for s in _TRANSFER_SUFFIXES)


def _has_transfer_suffix(name: str) -> bool:
    """大小写不敏感地判断文件名是否带 AleRax 转移频率文件后缀。

    ：先剥掉压缩后缀，故 ``FAM_transfers.txt.gz`` 同样被识别为转移文件。
    """
    lowered = logical_path(name).lower()
    return any(lowered.endswith(s) for s in _TRANSFER_SUFFIXES_LOWER)


def _family_name(basename: str) -> str:
    """从文件名提取家族名（剥离任一已知的转移文件后缀，大小写不敏感）。

    ：先剥压缩后缀，故 ``FAM_transfers.txt.gz`` 与 ``FAM_transfers.txt``
    得到**同一个**家族名。
    """
    basename = os.path.basename(logical_path(basename))
    for suffix in _TRANSFER_SUFFIXES_LOWER:
        if basename.lower().endswith(suffix):
            return basename[: -len(suffix)]
    return os.path.splitext(basename)[0]


def expand_alerax_inputs(paths: List[str]) -> List[str]:
    """把输入路径展开为 AleRax 转移频率文件（``*_transfers.txt`` /
    ``*_meanTransfers.txt``）列表。

    支持三种输入元素：

    * 直接的转移频率文件路径（原样保留）；
    * ``reconciliations/summaries`` 目录（glob 其中的转移频率文件）；
    * AleRax 输出根目录（自动拼接 ``reconciliations/summaries`` 后再 glob）。

    Args:
        paths: 输入路径列表（文件或目录混合）。

    Returns:
        排序去重后的转移频率文件路径列表。
    """
    result: List[str] = []
    seen = set()

    def _add(fp: str) -> None:
        ap = os.path.abspath(fp)
        if ap not in seen:
            seen.add(ap)
            result.append(fp)

    def _scan_dir(d: str) -> None:
        for name in sorted(os.listdir(d)):
            if _has_transfer_suffix(name):
                _add(os.path.join(d, name))

    for p in paths:
        if os.path.isdir(p):
            summaries = os.path.join(p, "reconciliations", "summaries")
            if os.path.isdir(summaries):
                _scan_dir(summaries)
            else:
                # 传入的可能已是 summaries 目录本身
                _scan_dir(p)
        elif os.path.isfile(p) and _has_transfer_suffix(os.path.basename(p)):
            _add(p)
        else:
            # 非约定后缀的文件：仍尝试作为转移文件解析（宽松）
            _add(p)
    return result


def _process_alerax_worker(args: Tuple) -> Tuple[Dict[str, float], str, Dict]:
    """模块级 worker（可 pickle）：解析单个 AleRax ``*_transfers.txt`` 文件。

    每一次 ``continue``（字段不足 / 频率非数字 / 端点不在物种树）都被计数，
    不再静默。

    Args:
        args: ``(file_path, valid_labels)``，可选第 3 元素 ``numeric_map``。

    Returns:
        ``(transfers, family, stats)``：``transfers`` 为
        ``{"donor,receptor": frequency}``；``family`` 为家族名（文件名去掉
        ``_transfers.txt`` 后缀）；``stats`` 为诊断字典。
    """
    file_path, valid_labels = tuple(args[:2])
    numeric_map = args[2] if len(args) > 2 and args[2] else None
    transfers: Dict[str, float] = {}
    stats = new_diagnostics("alerax", CONVENTION_DONOR_ITSELF, WEIGHT_POSTERIOR_FREQUENCY)
    stats["files_seen"] = 1
    stats["blocks_seen"] = 1  # AleRax 已在家族内预聚合：一文件=一份汇总
    family = _family_name(os.path.basename(file_path))

    # ：经共享读取器逐行解析，``*_transfers.txt.gz`` 可直接作为输入
    with open_text(file_path, encoding="utf-8") as fh:
        for line in fh:
            sp = line.split()
            if not sp:
                continue
            if len(sp) < 3:
                stats["malformed_lines"] += 1
                continue
            stats["transfers_seen"] += 1
            donor, receptor = sp[0], sp[1]
            try:
                freq = float(sp[2])
            except ValueError:
                stats["malformed_lines"] += 1
                continue
            d_lab, d_num = _resolve_species_label(donor, valid_labels, numeric_map)
            r_lab, r_num = _resolve_species_label(receptor, valid_labels, numeric_map)
            if d_lab is None or r_lab is None:
                stats["dropped_by_label_miss"] += 1
                samples = stats["unresolved_endpoint_samples"]
                for miss in (donor if d_lab is None else None, receptor if r_lab is None else None):
                    if miss and miss not in samples and len(samples) < 8:
                        samples.append(miss)
                continue
            stats["transfers_resolved"] += 1
            stats["numeric_id_resolutions"] += int(d_num) + int(r_num)
            if d_lab == r_lab:
                stats["self_loop_transfers"] += 1
            key = f"{d_lab},{r_lab}"
            transfers[key] = transfers.get(key, 0.0) + freq
    return transfers, family, stats


class AleRaxAdapter:
    """AleRax 输出 -> MaxTiC-Next 统一 ``ConstraintSet``。

    与其它适配器接口一致（``species_tree`` 在 ``convert`` 时传入）。因 AleRax
    每文件解析极轻量，默认顺序解析（可选进程池）。
    """

    def __init__(
        self,
        min_support: float = MINIMUM_SUPPORT_WITHIN_A_FAMILY,
        min_family_size: int = MINIMUM_FAMILY_SIZE,
        cache_dir: Optional[str] = None,
        min_endpoint_hit_rate: float = DEFAULT_MIN_ENDPOINT_HIT_RATE,
        quiet: bool = False,
    ) -> None:
        """初始化适配器。

        Args:
            min_support: 转移频率阈值，保留 ``frequency > min_support`` 的约束。
            min_family_size: 接口一致性占位（AleRax 已家族内聚合，本参数不生效）。
            cache_dir: 接口一致性占位（AleRax 单文件解析极快，暂不缓存）。
            min_endpoint_hit_rate: 端点命中率低于该值告警（默认 0.5）；命中率为 0
                抛 ``LabelMismatchError`` 而不是静默产出空约束集。
            quiet: 不打印 stderr 摘要。
        """
        self.min_support = min_support
        self.min_family_size = min_family_size
        self.cache_dir = cache_dir
        self.min_endpoint_hit_rate = min_endpoint_hit_rate
        self.quiet = quiet

    @staticmethod
    def _reciprocal_subtract(counts: Dict[str, float]) -> None:
        """互反减法（浮点版，委托共用实现；自环 ``"X,X"`` 不自我抵消，见 _diagnostics）。"""
        reciprocal_subtract(counts)

    def _merge_payload(
        self,
        cset: ConstraintSet,
        species_tree: Tree,
        payload: Tuple,
    ) -> None:
        """把一个家族的解析结果合并进 ``cset``（互反减法 + 频率阈值 + 距离标注）。"""
        transfers, family = payload[0], payload[1]
        if len(payload) >= 3 and isinstance(payload[2], dict):
            stats = payload[2]
        else:
            stats = new_diagnostics("alerax", CONVENTION_DONOR_ITSELF, WEIGHT_POSTERIOR_FREQUENCY)
            stats["transfers_resolved"] = len(transfers)
            stats["transfers_seen"] = len(transfers)
        cset.diagnostics.setdefault("files", []).append(stats)
        self._reciprocal_subtract(transfers)
        for key, freq in transfers.items():
            if freq <= 0.0:
                stats["dropped_by_support"] += 1
                continue
            weight, support = weight_and_support(stats, freq)
            if not passes_support(stats, weight, freq, self.min_support):
                stats["dropped_by_support"] += 1
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
                    },
                )
            )
            stats["constraints_kept"] += 1

    def convert(
        self,
        species_tree: Tree,
        inputs: List[str],
        output_path: Optional[str] = None,
        max_workers: Optional[int] = None,
        parallel: str = "process",
    ) -> ConstraintSet:
        """把 AleRax 输出（文件或目录）转换为 ``ConstraintSet``。

        Args:
            species_tree: 物种树。
            inputs: ``*_transfers.txt`` 文件或 AleRax 目录（自动展开）。
            output_path: 可选，若提供则把约束以 CSV 格式写出。
            max_workers: 并行进程数；``None`` 时取默认。
            parallel: 并行模式，``"process"``（默认）或 ``"thread"``；文件数少时自动顺序。

        Returns:
            统一 ``ConstraintSet``，每条约束 ``metadata`` 含
            ``family`` / ``support`` / ``distance``。
        """
        transfer_files = expand_alerax_inputs(inputs)
        if not transfer_files:
            if inputs:
                # 输入被识别为 AleRax 输出却展开不出任何转移频率文件（如现行版本的
                # ``_meanTransfers.txt`` 未被旧逻辑覆盖、或目录为空）：显式报错，
                # 绝不静默返回空约束集（否则下游会误以为"无转移"）。
                raise ValueError(
                    "未在输入中找到任何 AleRax 转移频率文件（*_transfers.txt / "
                    "*_meanTransfers.txt）：{!r}。请确认输入为 AleRax 输出目录"
                    "（含 reconciliations/summaries/）或直接指向转移频率文件。".format(inputs)
                )
            cset = ConstraintSet(
                diagnostics=new_diagnostics(
                    "alerax", CONVENTION_DONOR_ITSELF, WEIGHT_POSTERIOR_FREQUENCY
                )
            )
            finalize_adapter_diagnostics(
                cset, self.min_endpoint_hit_rate, set(), "alerax", quiet=True, files=0
            )
            if output_path is not None:
                self._write_constraints_file(cset, output_path)
            return cset

        valid_labels = _build_valid_labels(species_tree)
        worker_args = [(tf, valid_labels) for tf in transfer_files]

        payloads = map_payloads(
            _process_alerax_worker,
            worker_args,
            max_workers=max_workers,
            parallel=parallel if len(transfer_files) >= 4 else "sequential",
        )

        cset = ConstraintSet(
            diagnostics=new_diagnostics(
                "alerax", CONVENTION_DONOR_ITSELF, WEIGHT_POSTERIOR_FREQUENCY
            )
        )
        for payload in payloads:
            if payload is not None:
                self._merge_payload(cset, species_tree, payload)

        # 先落盘再跑守卫：0 命中会抛错，但 -o 产物仍须可核对
        if output_path is not None:
            self._write_constraints_file(cset, output_path)
        finalize_adapter_diagnostics(
            cset,
            self.min_endpoint_hit_rate,
            valid_labels,
            "alerax",
            quiet=self.quiet,
            files=len(transfer_files),
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


def convert_from_alerax(
    species_tree: Tree,
    inputs: List[str],
    min_support: float = MINIMUM_SUPPORT_WITHIN_A_FAMILY,
    min_family_size: int = MINIMUM_FAMILY_SIZE,
    cache_dir: Optional[str] = None,
    output_path: Optional[str] = None,
    max_workers: Optional[int] = None,
    parallel: str = "process",
    min_endpoint_hit_rate: float = DEFAULT_MIN_ENDPOINT_HIT_RATE,
    quiet: bool = False,
) -> ConstraintSet:
    """便捷函数：把 AleRax 输出（文件或目录）转换为 ``ConstraintSet``。

    Args:
        species_tree: 物种树。
        inputs: ``*_transfers.txt`` 文件或 AleRax 目录（自动展开）。
        min_support: 转移频率阈值。
        min_family_size: 接口一致性占位（不生效）。
        cache_dir: 接口一致性占位（暂不缓存）。
        output_path: 可选约束写出路径。
        max_workers: 并行解析进程/线程数。
        parallel: 并行模式，``"process"``（默认）或 ``"thread"``。
        min_endpoint_hit_rate: 端点命中率告警阈值（默认 0.5；0 命中必报错）。
        quiet: 不打印 stderr 摘要。

    Returns:
        统一 ``ConstraintSet``。
    """
    adapter = AleRaxAdapter(
        min_support=min_support,
        min_family_size=min_family_size,
        cache_dir=cache_dir,
        min_endpoint_hit_rate=min_endpoint_hit_rate,
        quiet=quiet,
    )
    return adapter.convert(
        species_tree,
        inputs,
        output_path=output_path,
        max_workers=max_workers,
        parallel=parallel,
    )
