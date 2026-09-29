"""适配器注册表与输出格式自动检测（多软件兼容框架的插件式接口）。

本模块把 5 个上游调和工具的适配器统一登记为**插件式条目**，并提供：

* ``AdapterRegistry``：工具名 -> 适配器元信息（``convert`` 便捷函数、输入描述、
  检测优先级）的注册表，支持 ``register`` 动态扩展；
* ``detect_format(path)`` / ``detect_formats(paths)``：按**文件扩展名 + 内容嗅探**
  自动识别上游工具（对应 ``--from-auto`` 与 ``convert_auto``）；
* ``convert(tool, species_tree, inputs, **kwargs)``：按工具名分发到对应适配器；
* ``convert_auto(species_tree, inputs, **kwargs)``：先检测再分发（``--from-auto``）。

统一约定：所有适配器最终产出统一 ``ConstraintSet``（``Constraint(donor, receptor,
weight, metadata)``），下游排序器与文本约束完全同构，实现"一次检测、统一管道"。

两个不可通约的口径（，机器可读）：每个 ``ConstraintSet.diagnostics``
都带 ``donor_endpoint_convention``（ALE = ``parent_of_donor``，其余四工具 =
``donor_itself``，同一事件相差一个物种层级）与 ``weight_semantics``
（``support_over_declared_samples`` / ``support_over_replicate_blocks`` /
``integer_count`` / ``posterior_frequency``）。``merge_constraint_sets`` 用于合并
多来源约束集，并在**混合口径**时打印警告 + 记入诊断；``convert_auto(...,
allow_mixed=True)`` 走的就是这条路径。语义详解见
``constraints/adapters/README_semantics.md``。

检测优先级（由具体到宽松）：
    ale (.uml_rec / T@,D@) > eccetera (recPhyloXML) > alerax (*_transfers.txt / 目录)
    > artra (Replacing/Additive Transfer) > ranger (Transfer + Recipient-->)。
"""

import inspect
import os
from typing import Callable, Dict, List, Optional

from maxtic_next.constraints.constraint import ConstraintSet
from maxtic_next.io.compression import iter_text_chunks, logical_path
from maxtic_next.tree.tree import Tree

from maxtic_next.constraints.adapters._diagnostics import (
    CONVENTION_DONOR_ITSELF,
    CONVENTION_PARENT_OF_DONOR,
    warn,
)
from maxtic_next.constraints.adapters.ale import convert_from_ale
from maxtic_next.constraints.adapters.ranger_dtl import convert_from_ranger_dtl
from maxtic_next.constraints.adapters.eccetera import convert_from_eccetera
from maxtic_next.constraints.adapters.artra import convert_from_artra
from maxtic_next.constraints.adapters.alerax import (
    convert_from_alerax,
    expand_alerax_inputs,  # noqa: F401  有意再导出：外部按 registry 统一入口取适配器符号
    _TRANSFER_SUFFIXES,
)


#: 各工具的供体端点层级约定：混用时必须警告，不可静默相加
ENDPOINT_CONVENTIONS: Dict[str, str] = {
    "ale": CONVENTION_PARENT_OF_DONOR,  # parent(donnor) / donnor_search
    "ranger": CONVENTION_DONOR_ITSELF,  # Mapping -->
    "eccetera": CONVENTION_DONOR_ITSELF,  # branchingOut speciesLocation
    "artra": CONVENTION_DONOR_ITSELF,  # Mapping -->（与 RANGER 同形）
    "alerax": CONVENTION_DONOR_ITSELF,  # 转移表第一列
}


class MixedFormatError(ValueError):
    """输入被嗅探出**多种**上游格式。

    除消息外还携带 ``groups``（工具名 -> 该工具的路径列表）：调用方若打算"分组各自
    解析后合并并警告"，可直接用它，**不必**把同一批文件再嗅探一遍。
    """

    def __init__(self, message: str, groups: Dict[str, List[str]]) -> None:
        super().__init__(message)
        self.groups: Dict[str, List[str]] = {k: list(v) for k, v in (groups or {}).items()}


def endpoint_convention(tool: str) -> Optional[str]:
    """返回工具的原生供体端点约定（未知工具返回 ``None``）。"""
    return ENDPOINT_CONVENTIONS.get(tool)


# 内容嗅探读取的最大字节数（避免整读超大文件）
_SNIFF_BYTES = 65536
# 内容嗅探的最大扫描字节数（分块扫描大文件；超过后按扩展名/目录结构判定）
_SNIFF_SCAN_LIMIT = 4 * 1024 * 1024


def _endswith_transfer_suffix(name: str) -> bool:
    """大小写不敏感地判断文件名是否为 AleRax 转移频率文件。"""
    lowered = name.lower()
    return any(lowered.endswith(s.lower()) for s in _TRANSFER_SUFFIXES)


class AdapterEntry:
    """单个适配器注册条目。"""

    def __init__(self, name: str, convert: Callable, description: str, input_desc: str) -> None:
        self.name = name
        self.convert = convert  # convert(species_tree, inputs, **kwargs)
        self.description = description  # 工具简介
        self.input_desc = input_desc  # 期望输入描述（CLI help / 文档）


class AdapterRegistry:
    """适配器注册表：工具名 -> ``AdapterEntry``（支持动态注册）。"""

    def __init__(self) -> None:
        self._entries: Dict[str, AdapterEntry] = {}

    def register(self, entry: AdapterEntry) -> None:
        """登记一个适配器条目（同名覆盖，支持插件式扩展）。"""
        self._entries[entry.name] = entry

    def get(self, name: str) -> AdapterEntry:
        """按工具名取条目；未知工具抛 ``KeyError``（附可用列表）。"""
        if name not in self._entries:
            raise KeyError(f"未知上游工具 {name!r}；可用：{sorted(self._entries)}")
        return self._entries[name]

    def names(self) -> List[str]:
        """返回已注册工具名（排序）。"""
        return sorted(self._entries)

    def describe(self) -> str:
        """返回可读的工具清单（供 CLI help / 文档），含各工具的供体端点约定。"""
        lines = []
        for name in self.names():
            e = self._entries[name]
            conv = endpoint_convention(name) or "unknown"
            lines.append(
                f"  {name:<9} {e.description}（输入：{e.input_desc}；donor 端点约定：{conv}）"
            )
        return "\n".join(lines)


# ----------------------------------------------------------------------
# 便捷分发包装：把各 convert_from_* 统一为 convert(species_tree, inputs, **kwargs)
# 各适配器签名略有差异（ale 有 source；artra 有 transfer_kind；alerax 无 family_size
# 语义），此处用公共参数子集 + 工具特定可选参数分发。
# ----------------------------------------------------------------------

#: 各适配器**共同**的可选参数名（并非每个适配器都接受全部：例如 ALE 忠实沿用原版
#: 的端点校验，没有 ``min_endpoint_hit_rate`` 守卫，见 ``_public_kwargs``）。
_COMMON_PARAM_NAMES = (
    "min_support",
    "min_family_size",
    "cache_dir",
    "output_path",
    "max_workers",
    "parallel",
    "min_endpoint_hit_rate",
    "quiet",
)


def _public_kwargs(
    target: Callable, kwargs: Dict, dropped_sink: Optional[List[str]] = None
) -> Dict:
    """挑出 ``target`` **确实接受**的公共参数。

    ``min_endpoint_hit_rate`` / ``quiet`` 是为  /  新增的开关；但"是否适用"
    取决于各适配器的端点校验口径：ALE 只要求 donor 能解析出父节点、且 receptor 不是
    现存物种（与原版 ``constraints_from_reconciliations.py`` 一致），因此它没有
    ``min_endpoint_hit_rate`` 参数。按**签名**过滤而不是盲转，避免 ``TypeError``，
    同时新增参数无需改动分发代码。

    签名过滤不得静默：给 ``api.rank`` 或某个
    适配器新增形参会悄悄改变其它适配器的入参集合。因此"被丢弃的公共参数"经
    ``dropped_sink`` 回报给调用方，由 ``_convert_*`` 写进 ``ConstraintSet.diagnostics``
    （与"一切丢弃都要可见"的原则一致），再由 CLI / 报告 / run_metadata 呈现。
    """
    params: Dict[str, inspect.Parameter]
    try:
        params = dict(inspect.signature(target).parameters)
    except (TypeError, ValueError):  # 内建 / 无签名对象：原样放行
        params = {}
    has_var_kw = any(p.kind is inspect.Parameter.VAR_KEYWORD for p in params.values())
    accepted = set(kwargs) if (not params or has_var_kw) else set(params)
    out = {}
    for k, v in kwargs.items():
        if k not in _COMMON_PARAM_NAMES or v is None:
            continue
        if k in accepted:
            out[k] = v
        elif dropped_sink is not None:
            dropped_sink.append(k)
    return out


def _note_dropped(cset, adapter: str, dropped: List[str]):
    """把"该适配器不适用的公共参数"记入诊断。"""
    if dropped:
        diag = getattr(cset, "diagnostics", None)
        if isinstance(diag, dict):
            prev = diag.get("adapter_params_not_applicable")
            note = f"{adapter} 不接受并已忽略公共参数：{', '.join(sorted(set(dropped)))}"
            diag["adapter_params_not_applicable"] = f"{prev}; {note}" if prev else note
    return cset


def _convert_ale(species_tree: Tree, inputs: List[str], **kwargs) -> ConstraintSet:
    dropped: List[str] = []
    ck = _public_kwargs(convert_from_ale, kwargs, dropped)
    if "ale_source" in kwargs and kwargs["ale_source"] is not None:
        ck["source"] = kwargs["ale_source"]
    return _note_dropped(convert_from_ale(species_tree, inputs, **ck), "ale", dropped)


def _convert_ranger(species_tree: Tree, inputs: List[str], **kwargs) -> ConstraintSet:
    dropped: List[str] = []
    ck = _public_kwargs(convert_from_ranger_dtl, kwargs, dropped)
    return _note_dropped(convert_from_ranger_dtl(species_tree, inputs, **ck), "ranger", dropped)


def _convert_eccetera(species_tree: Tree, inputs: List[str], **kwargs) -> ConstraintSet:
    dropped: List[str] = []
    ck = _public_kwargs(convert_from_eccetera, kwargs, dropped)
    return _note_dropped(convert_from_eccetera(species_tree, inputs, **ck), "eccetera", dropped)


def _convert_artra(species_tree: Tree, inputs: List[str], **kwargs) -> ConstraintSet:
    dropped: List[str] = []
    ck = _public_kwargs(convert_from_artra, kwargs, dropped)
    if "artra_transfer_kind" in kwargs and kwargs["artra_transfer_kind"] is not None:
        ck["transfer_kind"] = kwargs["artra_transfer_kind"]
    return _note_dropped(convert_from_artra(species_tree, inputs, **ck), "artra", dropped)


def _convert_alerax(species_tree: Tree, inputs: List[str], **kwargs) -> ConstraintSet:
    dropped: List[str] = []
    ck = _public_kwargs(convert_from_alerax, kwargs, dropped)
    return _note_dropped(convert_from_alerax(species_tree, inputs, **ck), "alerax", dropped)


# 默认全局注册表（登记 5 个内置适配器）
REGISTRY = AdapterRegistry()
REGISTRY.register(
    AdapterEntry(
        "ale",
        _convert_ale,
        "ALE_undated 分摊调和",
        ".uml_rec 采样调和文件（T@/D@ 事件；默认 --ale-source trf＝转移事件口径，"
        "即 ALE 官方 constraints_from_transfers；rec＝调和事件，可选）",
    )
)
REGISTRY.register(
    AdapterEntry(
        "ranger", _convert_ranger, "RANGER-DTLx DTL 调和", "调和报告文本（Mapping-->/Recipient-->）"
    )
)
REGISTRY.register(
    AdapterEntry("eccetera", _convert_eccetera, "ecceTERA DTL 调和", "recPhyloXML 文件")
)
REGISTRY.register(
    AdapterEntry(
        "artra",
        _convert_artra,
        "ARTra 加性/替换转移分类",
        "output.txt 调和报告（Replacing/Additive Transfer）",
    )
)
REGISTRY.register(
    AdapterEntry(
        "alerax",
        _convert_alerax,
        "AleRax 贝叶斯调和",
        "reconciliations/summaries/*_transfers.txt（或 *_meanTransfers.txt）或 AleRax 目录",
    )
)


# ----------------------------------------------------------------------
# 自动检测
# ----------------------------------------------------------------------


def detect_format(path: str) -> Optional[str]:
    """检测单个路径（文件或目录）对应的上游工具名；无法判定返回 ``None``。

    检测顺序：目录（AleRax 结构）> 扩展名 > 内容嗅探。检测**不消费**内容语义，
    仅做工具识别；真正解析由对应适配器完成。

    Args:
        path: 文件或目录路径。

    Returns:
        工具名（``"ale"`` / ``"ranger"`` / ``"eccetera"`` / ``"artra"`` /
        ``"alerax"``）或 ``None``（未识别，通常为文本约束文件）。
    """
    # 1) 目录：AleRax 输出根目录或 summaries 目录
    # （注：Ranger / ecceTERA 的结果目录无统一结构约定，暂不识别为目录输入。）
    if os.path.isdir(path):
        summaries = os.path.join(path, "reconciliations", "summaries")
        if os.path.isdir(summaries):
            return "alerax"
        try:
            if any(_endswith_transfer_suffix(n) for n in os.listdir(path)):
                return "alerax"
        except OSError:
            return None
        return None

    # 2) 扩展名（先剥压缩后缀，故 ``x.uml_rec.gz`` 仍判为 ale）
    lower = logical_path(path).lower()
    if lower.endswith(".uml_rec"):
        return "ale"
    if _endswith_transfer_suffix(os.path.basename(logical_path(path))):
        return "alerax"
    if lower.endswith(".recphyloxml"):
        return "eccetera"

    # 3) 内容嗅探（分块扫描，避免 64KB 截断漏检文件中后部的特征标记）
    for text in _sniff_chunks(path):
        if "<recPhylo" in text or "<recGeneTree" in text or "recPhyloXML" in text:
            return "eccetera"
        if "T@" in text and ("reconciled" in text or "D@" in text):
            return "ale"
        if "Recipient" in text and "-->" in text:
            if "Replacing Transfer" in text or "Additive Transfer" in text:
                return "artra"
            return "ranger"
    return None


def _sniff_chunks(path: str):
    """按块产出**解压后**的文件内容用于嗅探；总扫描量不超过 ``_SNIFF_SCAN_LIMIT``。

    块间有少量重叠（保留上一块尾部 64 字节），保证跨块边界的标记不被截断。
    读失败时不产出任何内容。

    ：走 :func:`maxtic_next.io.compression.iter_text_chunks`，故
    ``x.uml_rec.gz`` / gzip 流嗅探到的是**解压后**的文本（旧实现会把压缩字节读成
    乱码而判为"未识别"）；多成员 ``.tgz`` 归档则抛出带 ``tar -xzOf`` 指引的
    ``CompressedArchiveError``（``ValueError`` 子类），而不是静默判为未识别。
    """
    tail = ""
    scanned = 0
    for chunk in iter_text_chunks(path, _SNIFF_BYTES, encoding="utf-8", errors="replace"):
        scanned += len(chunk)
        yield tail + chunk
        tail = chunk[-64:]
        if scanned >= _SNIFF_SCAN_LIMIT:
            break


def detect_formats(paths: List[str]) -> Optional[str]:
    """对多个路径检测统一工具名；若各路径检测结果冲突则抛 ``ValueError``。

    Args:
        paths: 输入路径列表。

    Returns:
        统一工具名，或 ``None``（全部未识别）。

    Raises:
        MixedFormatError: 检测到多种不同工具（``--from-auto`` 会据此分组解析后
            **警告并合并**；显式 ``--from <tool>`` 的单工具路径不会走到这里）。
            它是 ``ValueError`` 的子类，故既有 ``except ValueError`` 语义不变。
    """
    detected: Dict[str, List[str]] = {}
    for p in paths:
        tool = detect_format(p)
        if tool is not None:
            detected.setdefault(tool, []).append(p)
    if not detected:
        return None
    if len(detected) > 1:
        raise MixedFormatError(
            f"输入混合了多种上游格式：{ {k: len(v) for k, v in detected.items()} }；"
            f"请用 --from <tool> 显式指定单一工具。",
            detected,
        )
    return next(iter(detected))


def convert(tool: str, species_tree: Tree, inputs: List[str], **kwargs) -> ConstraintSet:
    """按工具名分发到对应适配器，返回统一 ``ConstraintSet``。

    返回集合的 ``diagnostics`` 会补上 ``tool`` 与 ``donor_endpoint_convention``
    ，供上层（``api`` / ``--from-auto``）判断口径是否一致。

    Args:
        tool: 工具名（见 ``REGISTRY.names()``）。
        species_tree: 物种树。
        inputs: 输入文件/目录路径列表。
        **kwargs: 通用参数（``min_support`` / ``min_family_size`` / ``cache_dir`` /
            ``output_path`` / ``max_workers`` / ``parallel`` /
            ``min_endpoint_hit_rate`` / ``quiet``）与工具特定参数
            （``ale_source`` / ``artra_transfer_kind``）。

    Returns:
        统一 ``ConstraintSet``。
    """
    entry = REGISTRY.get(tool)
    cset = entry.convert(species_tree, inputs, **kwargs)
    diag = getattr(cset, "diagnostics", None)
    if diag is None:
        diag = {}
        cset.diagnostics = diag
    diag.setdefault("tool", tool)
    diag["donor_endpoint_convention"] = (
        diag.get("donor_endpoint_convention")
        or endpoint_convention(tool)
        or CONVENTION_DONOR_ITSELF
    )
    diag.setdefault("source", kwargs.get("ale_source") if tool == "ale" else None)
    return cset


def merge_constraint_sets(
    constraint_sets: List[ConstraintSet], quiet: bool = False
) -> ConstraintSet:
    """把多个来源的 ``ConstraintSet`` 合并为一个（供 ``--from-auto`` 混合输入使用）。

    合并前做两项**口径一致性**检查，任一不一致都会在 stderr 警告并写入返回集合的
    ``diagnostics["convention_conflict"]`` / ``["semantics_conflict"]``：

    * ``donor_endpoint_convention`` 不同：ALE 的 donor 是"供体的父"，
      其余工具是"供体本身"——同一边键在两种约定下指向**不同**的物种节点，
      直接相加会让同一 edge 混入两代端点；
    * ``weight_semantics`` 不同：[0,1] 支持度、整数计数、后验频率
      量纲不可通约，相加后的权重不再可比。

    Args:
        constraint_sets: 待合并的约束集列表。
        quiet: 只记诊断、不打印 stderr。

    Returns:
        合并后的 ``ConstraintSet``（``diagnostics["merged_from"]`` 记录各来源）。
    """
    merged = ConstraintSet()
    conventions: Dict[str, int] = {}
    semantics: Dict[str, int] = {}
    for cset in constraint_sets:
        for c in cset.constraints:
            merged.add(c)
        diag = getattr(cset, "diagnostics", {}) or {}
        tool = str(diag.get("tool", "unknown"))
        conv = diag.get("donor_endpoint_convention") or endpoint_convention(tool)
        sem = diag.get("weight_semantics")
        conventions[str(conv)] = conventions.get(str(conv), 0) + 1
        semantics[str(sem)] = semantics.get(str(sem), 0) + 1
        merged.diagnostics.setdefault("merged_from", []).append(
            {
                "tool": tool,
                "donor_endpoint_convention": conv,
                "weight_semantics": sem,
                "constraints": len(cset.constraints),
            }
        )
        merge_into = merged.diagnostics.setdefault("aggregated", {})
        _accumulate(merge_into, diag)

    if len(conventions) > 1:
        msg = (
            "待合并的约束集来自**不同的供体端点层级约定**："
            f"{sorted(conventions)}（ALE 取 parent(donor)，"
            "RANGER-DTLx / ARTra / ecceTERA / AleRax 取 donor 本身）。"
            "同一边键在两代端点上含义不同，相加后的约束集不再自洽："
            "请分别运行后按工具比较，或只喂同一种约定的上游输出。"
        )
        merged.diagnostics["convention_conflict"] = sorted(conventions)
        if not quiet:
            warn(msg, "auto")
    if len(semantics) > 1:
        msg = (
            "待合并的约束集权重口径不一致："
            f"{sorted(semantics)}——[0,1] 支持度、整数计数与后验频率不可通约，"
            "相加后的权重没有概率含义。"
        )
        merged.diagnostics["semantics_conflict"] = sorted(semantics)
        if not quiet:
            warn(msg, "auto")
    merged.diagnostics.setdefault(
        "donor_endpoint_convention", next(iter(conventions)) if len(conventions) == 1 else "mixed"
    )
    merged.diagnostics.setdefault(
        "weight_semantics", next(iter(semantics)) if len(semantics) == 1 else "mixed"
    )
    return merged


def _accumulate(dst: Dict, src: Dict) -> Dict:
    """把 ``src`` 的计数器累加进 ``dst``（合并多来源时只保留可加项）。"""
    for k, v in src.items():
        if isinstance(v, bool):
            continue
        if isinstance(v, int):
            dst[k] = int(dst.get(k, 0)) + v
        elif isinstance(v, float):
            dst[k] = max(float(dst.get(k, 0.0)), v)
        elif isinstance(v, list):
            cur = dst.setdefault(k, [])
            for item in v:
                if item not in cur:
                    cur.append(item)
    return dst


def convert_auto(species_tree: Tree, inputs: List[str], **kwargs) -> ConstraintSet:
    """自动检测上游格式后分发（``--from-auto`` / ``--from auto``）。

    默认要求全部输入同格式（与历史行为一致：混合格式抛 ``ValueError``）。
    传 ``allow_mixed=True`` 时改为**按检测到的工具分组**、逐组用各自适配器解析、
    再经 ``merge_constraint_sets`` 合并——此时若两组的供体端点层级或权重口径不同
    （例如 ALE + RANGER），会打印警告并记入诊断。

    Args:
        species_tree: 物种树。
        inputs: 输入文件/目录路径列表。
        **kwargs: 同 ``convert``；额外支持 ``allow_mixed``。

    Returns:
        统一 ``ConstraintSet``。

    Raises:
        ValueError: 无法识别任一输入；或（``allow_mixed`` 未开时）输入混合多种格式。
    """
    allow_mixed = bool(kwargs.pop("allow_mixed", False))
    if allow_mixed:
        groups = _group_by_format(inputs)
        if not groups:
            raise ValueError(
                "无法自动识别上游工具格式；请用 --from <tool> 显式指定"
                f"（可用：{REGISTRY.names()}）。"
            )
        sets = [
            convert(tool, species_tree, paths, **kwargs) for tool, paths in sorted(groups.items())
        ]
        # 冲突（口径混用）不受 quiet 影响：quiet 只抑制逐适配器的例行摘要。
        return merge_constraint_sets(sets)
    tool = detect_formats(inputs)
    if tool is None:
        raise ValueError(
            f"无法自动识别上游工具格式；请用 --from <tool> 显式指定（可用：{REGISTRY.names()}）。"
        )
    return convert(tool, species_tree, inputs, **kwargs)


def _group_by_format(paths: List[str]) -> Dict[str, List[str]]:
    """按 ``detect_format`` 把输入分组（未识别者忽略）；不重复嗅探同一目录。"""
    groups: Dict[str, List[str]] = {}
    for p in paths:
        tool = detect_format(p)
        if tool:
            groups.setdefault(tool, []).append(p)
    return groups
