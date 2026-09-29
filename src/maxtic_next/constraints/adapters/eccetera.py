"""ecceTERA 调和结果适配器。

ecceTERA 的结构化调和输出为 **recPhyloXML**（Duchemin et al. 2018 标准；
``--recPhyloXML.reconciliation=true``）。经 ecceTERA 源码
``src/DTLGraph.cpp:getRecPhyloXMLReconciliation`` 求证，其转移事件编码为：

* **供体**节点写 ``<branchingOut speciesLocation="D">``（D = 供体物种，转移离开处）；
* 被转移到的**受体**子节点，其 ``<eventsRec>`` 以
  ``<transferBack destinationSpecies="R">`` 起始（R = 受体物种，转移到达处），
  随后是该受体节点的真实事件（speciation / leaf / ...）。

因此抽取规则为：对每个含 ``transferBack destinationSpecies="R"`` 的 clade，其**父
clade** 的 ``branchingOut speciesLocation="D"`` 即供体；产出统一 ``(donor, receptor)
= (D, R)`` 有向约束（donor 应早于/不晚于 receptor）。

recPhyloXML 片段示例::

    <recPhylo>
     <recGeneTree><phylogeny rooted="true">
      <clade><name>5</name>
       <eventsRec><branchingOut speciesLocation="H12"></branchingOut></eventsRec>
       <clade><name>6</name>
        <eventsRec><speciation speciesLocation="H12"></speciation></eventsRec>
       </clade>
       <clade><name>7</name>
        <eventsRec><transferBack destinationSpecies="H45"></transferBack>
                   <leaf speciesLocation="H45" geneName="g7"></leaf></eventsRec>
       </clade>
      </clade>
     </phylogeny></recGeneTree>
    </recPhylo>

上例产出约束 ``H12 -> H45``（供体 H12、受体 H45）。

内部物种节点的数字 ID
----------------------------------------------------------------
ecceTERA 对**内部**物种节点写的是它自己的 ``node->getId()``：
``DTLGraph.cpp:2818-2821`` 对叶子用 ``->getName()``、对内部节点用
``->getId()``；``DTLGraph.cpp:1820`` 的 ``x->getFather()->getId()`` 同理。该编号
规则求证自 ``MySpeciesTree.cpp:112-136``（``assignPostOrderIds``，注释
"Assign ids using breadth-first postorder"）：**排序后的叶子先取 0..L-1，
其后各代节点按"自底向上、逐层广度优先"依次编号**。

因此本适配器：

1. 先按标签直接匹配用户物种树（叶子名 + 内部 bootstrap 标签）；
2. 匹配不上且端点是纯数字时，按上述规则把数字 ID **反解**回物种树标签
   （``_species.numeric_species_id_map``），使真实 ecceTERA 输出可用；
3. 反解仍失败 ⇒ 计数并计算"端点命中率"：**命中率为 0 直接抛**
   ``LabelMismatchError``（错误信息同时给出两侧的标签样本），低于可配置阈值
   （默认 50%）则在诊断 + stderr 中告警——绝不再静默返回 ``{}``；
4. ``min_family_size`` 依赖 ``<leaf>`` 元素计数：探测失败（0 个叶子）时**跳过**
   过滤并告警，而不是把所有转移清空。

XML 解析失败
------------------------
旧实现 ``except ET.ParseError: return {}, 0, 0``，把"文件损坏"与"该工具没检出
转移"混为一谈。现由 ``_process_eccetera_worker`` 抛 ``UpstreamParseError``，
异常信息含**文件名、行列号与原始异常**；纯函数 ``parse_recphyloxml`` 保持旧签名
（默认仍返回空三元组，但会向 stderr 打出**响亮的警告**，``strict=True`` 则抛错）。

口径：``donor_endpoint_convention = "donor_itself"``（供体本身，
ALE 路径则取"供体的父"）；权重口径为"该文件的 ``<recGeneTree>`` 块数 / 文件声明的
场景数"，无法确立时退化为**整数计数**，全部记录在 ``ConstraintSet.diagnostics``。

上述格式已按 ecceTERA 源码与真实输出（recPhyloXML）逐项求证。

接口与其它适配器一致（``species_tree`` 在 ``convert`` 时传入），支持
``ProcessPoolExecutor`` 并行解析和文件级 pickle 缓存（断点续传）。
"""

import hashlib
import os
import pickle
import xml.etree.ElementTree as ET
from typing import Dict, List, Optional, Set, Tuple

from maxtic_next.constraints.constraint import Constraint, ConstraintSet
from maxtic_next.constraints.adapters._diagnostics import (
    CACHE_KEY_VERSION,
    CONVENTION_DONOR_ITSELF,
    DEFAULT_MIN_ENDPOINT_HIT_RATE,
    WEIGHT_BLOCK_SUPPORT,
    WEIGHT_DECLARED_SUPPORT,
    WEIGHT_INTEGER_COUNT,
    UpstreamParseError,
    finalize_adapter_diagnostics,
    map_payloads,
    new_diagnostics,
    passes_support,
    reciprocal_subtract,
    warn,
    weight_and_support,
)
from maxtic_next.constraints.adapters._recon_report import (
    declared_sample_counts,
    resolve_family_size_filter,
)
from maxtic_next.constraints.adapters._species import (
    build_valid_labels as _build_valid_labels,
    distance_from as _distance_from,
    numeric_species_id_map as _numeric_species_id_map,
    resolve_species_label as _resolve_species_label,
    species_cache_identity as _species_cache_identity,
    species_labels_digest,
)
from maxtic_next.io.compression import logical_stem, read_text
from maxtic_next.tree.tree import Tree


# 与 ALE 一致的默认阈值
MINIMUM_SUPPORT_WITHIN_A_FAMILY = 0.05
MINIMUM_FAMILY_SIZE = 5


# ----------------------------------------------------------------------
# recPhyloXML 解析（模块级纯函数，可 pickle）
# ----------------------------------------------------------------------


def _local(tag: str) -> str:
    """去除 XML 命名空间前缀，返回本地标签名。"""
    return tag.rsplit("}", 1)[-1] if "}" in tag else tag


def _events_of(clade: ET.Element) -> List[Tuple[str, Dict[str, str]]]:
    """返回 clade 直属 ``<eventsRec>`` 内的事件列表 ``[(localtag, attrib), ...]``。"""
    events: List[Tuple[str, Dict[str, str]]] = []
    for child in clade:
        if _local(child.tag) == "eventsRec":
            for ev in child:
                events.append((_local(ev.tag), dict(ev.attrib)))
    return events


def _child_clades(clade: ET.Element) -> List[ET.Element]:
    """返回 clade 的直属子 ``<clade>`` 元素（不含 ``eventsRec`` 内的元素）。"""
    return [c for c in clade if _local(c.tag) == "clade"]


def _walk_clade(
    clade: ET.Element,
    parent_donor: Optional[str],
    valid_labels: Set[str],
    transfers: Dict[str, int],
    leaf_ids: Set[str],
    stats: Optional[Dict] = None,
    numeric_map: Optional[Dict[str, str]] = None,
) -> None:
    """递归遍历 clade 树，抽取转移约束、统计叶子数并**逐项计数**。

    Args:
        clade: 当前 ``<clade>`` 元素。
        parent_donor: 父 clade 的 ``branchingOut`` 供体物种（无则 ``None``）。
        valid_labels: 物种树合法标签集合（端点校验）。
        transfers: 累加的 ``{"donor,receptor": count}``（原地修改）。
        leaf_ids: 累加的基因叶子标识集合（家族规模，原地修改）。
        stats: 诊断字典（原地计数；``None`` 时不计数，便于兼容旧调用）。
        numeric_map: ``数字ID -> 物种树标签`` 反查表。
    """
    donor_for_children: Optional[str] = None
    transferback_dest: Optional[str] = None
    for tag, attrib in _events_of(clade):
        if tag == "leaf":
            gid = attrib.get("geneName") or attrib.get("speciesLocation")
            if gid:
                leaf_ids.add(gid)
        elif tag == "branchingOut":
            donor_for_children = attrib.get("speciesLocation")
        elif tag == "transferBack":
            transferback_dest = attrib.get("destinationSpecies")

    # 受体 clade（含 transferBack）与父 clade 的供体（branchingOut）配对
    if transferback_dest is not None:
        if parent_donor is None:
            # 根 clade 即为受体：上游未给出供体谱系，无法产出约束（计数而非静默）
            if stats is not None:
                stats["transfers_unpaired"] = int(stats.get("transfers_unpaired", 0)) + 1
        else:
            if stats is not None:
                stats["transfers_seen"] = int(stats.get("transfers_seen", 0)) + 1
            d_lab, d_num = _resolve_species_label(parent_donor, valid_labels, numeric_map)
            r_lab, r_num = _resolve_species_label(transferback_dest, valid_labels, numeric_map)
            if d_lab is None or r_lab is None:
                if stats is not None:
                    stats["dropped_by_label_miss"] = int(stats.get("dropped_by_label_miss", 0)) + 1
                    samples = stats.setdefault("unresolved_endpoint_samples", [])
                    for miss in (
                        parent_donor if d_lab is None else None,
                        transferback_dest if r_lab is None else None,
                    ):
                        if miss and miss not in samples and len(samples) < 8:
                            samples.append(miss)
            else:
                if stats is not None:
                    stats["transfers_resolved"] = int(stats.get("transfers_resolved", 0)) + 1
                    stats["numeric_id_resolutions"] = (
                        int(stats.get("numeric_id_resolutions", 0)) + int(d_num) + int(r_num)
                    )
                    if d_lab == r_lab:
                        stats["self_loop_transfers"] = int(stats.get("self_loop_transfers", 0)) + 1
                key = f"{d_lab},{r_lab}"
                transfers[key] = transfers.get(key, 0) + 1

    for child in _child_clades(clade):
        _walk_clade(
            child, donor_for_children, valid_labels, transfers, leaf_ids, stats, numeric_map
        )


def _parse_core(
    text: str, valid_labels: Set[str], numeric_map: Optional[Dict[str, str]] = None
) -> Tuple[Dict[str, int], int, int, Dict, Optional[Dict]]:
    """recPhyloXML 解析核心；返回 ``(transfers, blocks, family_size, stats, error)``。

    ``error`` 非 ``None`` 时是一个 ``{"message", "line", "column"}`` 字典：调用方
    必须据此抛错（适配器路径）或至少告警（旧的兼容纯函数路径）——。
    """
    transfers: Dict[str, int] = {}
    blocks = 0
    family_size = 0
    stats = new_diagnostics("eccetera", CONVENTION_DONOR_ITSELF, WEIGHT_BLOCK_SUPPORT)

    # 安全加固：recPhyloXML 无 DTD/实体声明；若出现 DOCTYPE/ENTITY 则拒绝解析，
    # 规避内部实体展开（billion laughs）与外部实体（XXE）风险（stdlib 无需 defusedxml）。
    lowered = text.lstrip()[:512].lower()
    if "<!doctype" in lowered or "<!entity" in lowered:
        stats["rejected_unsafe_markup"] = True
        return (
            transfers,
            blocks,
            family_size,
            stats,
            {
                "message": "文档含 DOCTYPE/ENTITY 声明，已被安全策略拒绝解析"
                "（recPhyloXML 不需要 DTD；防实体展开/XXE）",
                "line": 0,
                "column": 0,
            },
        )

    try:
        root = ET.fromstring(text)
    except ET.ParseError as exc:
        line, column = _extract_position(exc)
        stats["parse_failed"] = True
        return (
            transfers,
            blocks,
            family_size,
            stats,
            {
                "message": "recPhyloXML 解析失败（文件损坏 / 被截断 / 非 XML）",
                "line": line,
                "column": column,
                "cause": exc,
            },
        )

    # 收集所有 recGeneTree（兼容根即 recPhylo 或直接为 recGeneTree）
    rec_gene_trees: List[ET.Element] = []
    if _local(root.tag) == "recGeneTree":
        rec_gene_trees = [root]
    else:
        for el in root.iter():
            if _local(el.tag) == "recGeneTree":
                rec_gene_trees.append(el)

    for rgt in rec_gene_trees:
        # 定位该 recGeneTree 下的根 clade（phylogeny 下的首个 clade）
        root_clades: List[ET.Element] = []
        for el in rgt:
            if _local(el.tag) == "phylogeny":
                root_clades = _child_clades(el)
                break
        if not root_clades:
            root_clades = _child_clades(rgt)
        leaf_ids: Set[str] = set()
        for rc in root_clades:
            _walk_clade(rc, None, valid_labels, transfers, leaf_ids, stats, numeric_map)
        blocks += 1
        if len(leaf_ids) > family_size:
            family_size = len(leaf_ids)

    stats["blocks_seen"] = blocks
    stats["family_size_probe_ok"] = family_size > 0
    stats["files_seen"] = 1

    # ---- 权重口径：优先工具自己声明的场景数，否则用 <recGeneTree> 块数 ----
    declared = declared_sample_counts(text.splitlines())
    stats["declared_sample_counts"] = list(declared)
    if declared:
        stats["sample_denominator"] = max(declared)
        stats["weight_semantics"] = WEIGHT_DECLARED_SUPPORT
    elif blocks > 0:
        stats["sample_denominator"] = blocks
        stats["weight_semantics"] = WEIGHT_BLOCK_SUPPORT
        stats["warnings"].append(
            f"recPhyloXML 未声明场景总数，按 {blocks} 个 <recGeneTree> 块作分母；"
            "若 ecceTERA 实际枚举了更多最优调和，本次支持度会被高估。"
        )
    else:
        stats["sample_denominator"] = None
        stats["weight_semantics"] = WEIGHT_INTEGER_COUNT
    stats["min_support_is_fractional"] = stats["weight_semantics"] in (
        WEIGHT_DECLARED_SUPPORT,
        WEIGHT_BLOCK_SUPPORT,
    ) and bool(stats["sample_denominator"])
    if stats.get("numeric_id_resolutions"):
        stats["warnings"].append(
            f"{stats['numeric_id_resolutions']} 个端点由数字物种节点 ID 反解为"
            "物种树标签（ecceTERA 对内部物种节点写 getId()）。反解依赖"
            "自底向上逐层编号规则；若两棵树的拓扑不同，请核对 ecceTERA 所用"
            "物种树与传入 MaxTiC 的物种树是否同一棵。"
        )
    return transfers, blocks, family_size, stats, None


def _extract_position(exc: BaseException) -> Tuple[int, int]:
    """从 ``ET.ParseError`` 中取 ``(line, column)``（取不到则 0,0）。"""
    pos = getattr(exc, "position", None)
    if isinstance(pos, tuple) and len(pos) == 2:
        try:
            return int(pos[0]), int(pos[1])
        except (TypeError, ValueError):
            return 0, 0
    return 0, 0


def parse_recphyloxml(
    text: str,
    valid_labels: Set[str],
    strict: bool = False,
    numeric_map: Optional[Dict[str, str]] = None,
    source: str = "<字符串>",
) -> Tuple[Dict[str, int], int, int]:
    """解析 recPhyloXML 文本，抽取转移约束（兼容旧的三值返回）。

    每个 ``<recGeneTree>`` 视为一次调和块（支持度分母的候选）；块内按边键计数、
    块间求和。仅保留可解析到 ``valid_labels`` 的转移（数字 ID 端点会先经
    ``numeric_map`` 反解）。

    Args:
        text: recPhyloXML 文件内容。
        valid_labels: 物种树合法标签集合。
        strict: ``True`` 时 XML 解析失败/被安全拒绝会抛 ``UpstreamParseError``；
            ``False``（默认，兼容旧调用）仍返回空三元组，但会向 stderr 打出
            含文件、行列与原始异常的**响亮警告**。
        numeric_map: ``数字ID -> 标签`` 反查表。
        source: 文件名（仅用于告警/异常文案）。

    Returns:
        ``(transfers, blocks, family_size)``。
    """
    transfers, blocks, family_size, stats, error = _parse_core(text, valid_labels, numeric_map)
    if error is not None:
        msg = (
            f"{source}：{error['message']}"
            + (f"（第 {error['line']} 行，第 {error['column']} 列）" if error.get("line") else "")
            + (
                f"；原始异常：{type(error['cause']).__name__}: {error['cause']}"
                if error.get("cause") is not None
                else ""
            )
            + "。此情形**不等于**该工具没检出转移：损坏/被截断的 XML 旧实现会"
            "静默当作空结果。"
        )
        if strict:
            raise UpstreamParseError(
                f"{source}：{error['message']}",
                source=source,
                line=int(error.get("line") or 0),
                column=int(error.get("column") or 0),
                cause=error.get("cause"),
                diagnostics=stats,
            )
        warn(msg, "eccetera")
    return transfers, blocks, family_size


def _process_eccetera_worker(args: Tuple) -> Tuple[Dict[str, int], int, int, str, Dict]:
    """模块级 worker（可 pickle）：解析单个 ecceTERA recPhyloXML 文件（含缓存）。

    Args:
        args: ``(file_path, min_family_size, cache_dir, min_support,
               valid_labels[, species_identity[, numeric_map]])``。

    Returns:
        ``(transfers, blocks, family_size, family, stats)`` 五元组。

    Raises:
        UpstreamParseError: XML 解析失败或被安全策略拒绝。
    """
    (file_path, min_family_size, cache_dir, min_support, valid_labels) = tuple(args[:5])
    species_identity = args[5] if len(args) > 5 else None
    numeric_map = args[6] if len(args) > 6 and args[6] else None

    cache_file: Optional[str] = None
    if cache_dir:
        st = os.stat(file_path)
        # 缓存载荷依赖物种树（端点过滤 + 数字 ID 反查表 = 拓扑），故键须同时含
        # 标签摘要与**拓扑指纹**；版本号使旧条目永不被读到。
        identity = species_identity or f"labels:{species_labels_digest(valid_labels)}"
        raw = (
            f"v{CACHE_KEY_VERSION}|{os.path.abspath(file_path)}|{st.st_mtime}|"
            f"{st.st_size}|{min_support}|{min_family_size}|{identity}"
        )
        key = hashlib.sha256(raw.encode("utf-8")).hexdigest()
        cache_file = os.path.join(cache_dir, key + ".pkl")
        if os.path.exists(cache_file):
            try:
                with open(cache_file, "rb") as fh:
                    payload = pickle.load(fh)
                if len(payload) >= 5:
                    return payload
                stats = new_diagnostics("eccetera", CONVENTION_DONOR_ITSELF)
                stats["blocks_seen"] = payload[1]
                stats["cache_legacy_payload"] = True
                return tuple(payload) + (stats,)
            except Exception:
                pass

    # ：recPhyloXML 常以 .xml.gz / .recphyloxml.gz 分发（如 ecceTERA 示例包），
    # 经共享读取器直接解压；家族名取"逻辑词干"，与未压缩输入一致。
    text = read_text(file_path, encoding="utf-8")
    family = logical_stem(file_path)

    transfers, blocks, family_size, stats, error = _parse_core(text, valid_labels, numeric_map)
    if error is not None:
        # ：解析失败绝不表现为"没有约束"
        raise UpstreamParseError(
            f"{file_path}：{error['message']}",
            source=file_path,
            line=int(error.get("line") or 0),
            column=int(error.get("column") or 0),
            cause=error.get("cause"),
            diagnostics=stats,
        )
    stats["files_seen"] = 1

    if resolve_family_size_filter(
        stats, family_size, min_family_size, probe_desc="``<leaf …>`` 元素"
    ):
        stats["families_skipped"] = 1
        transfers = {}

    payload = (transfers, blocks, family_size, family, stats)

    if cache_dir and cache_file is not None:
        os.makedirs(cache_dir, exist_ok=True)
        with open(cache_file, "wb") as fh:
            pickle.dump(payload, fh)

    return payload


class EcceTERAAdapter:
    """ecceTERA recPhyloXML 调和结果 -> MaxTiC-Next 统一 ``ConstraintSet``。

    与 ``RangerDTLAdapter`` 接口一致（``species_tree`` 在 ``convert`` 时传入），
    支持 ``ProcessPoolExecutor`` 并行解析和文件级 pickle 缓存。
    """

    def __init__(
        self,
        min_support: float = MINIMUM_SUPPORT_WITHIN_A_FAMILY,
        min_family_size: int = MINIMUM_FAMILY_SIZE,
        cache_dir: Optional[str] = None,
        min_endpoint_hit_rate: float = DEFAULT_MIN_ENDPOINT_HIT_RATE,
        resolve_numeric_ids: bool = True,
        quiet: bool = False,
    ) -> None:
        """初始化适配器。

        Args:
            min_support: 单基因家族内最小支持度阈值（默认 0.05）。
            min_family_size: 参与计算的最小基因家族规模（叶子数，默认 5；
                ``<leaf>`` 探测失败时跳过过滤并告警）。
            cache_dir: 基因家族文件级缓存目录；``None`` 表示不缓存。缓存键含
                物种树拓扑指纹与版本号。
            min_endpoint_hit_rate: 端点命中率低于该值告警（默认 0.5）；命中率为
                0 抛 ``LabelMismatchError``。``<= 0`` 只关闭告警，
                0 命中仍报错。
            resolve_numeric_ids: 是否把数字物种节点 ID 反解为物种树标签
                （默认 ``True``；关掉即旧行为：数字端点一律丢弃）。
            quiet: 不打印 stderr 摘要。
        """
        self.min_support = min_support
        self.min_family_size = min_family_size
        self.cache_dir = cache_dir
        self.min_endpoint_hit_rate = min_endpoint_hit_rate
        self.resolve_numeric_ids = resolve_numeric_ids
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
        """把一个家族的解析结果合并进 ``cset``（互反减法 + 支持度过滤 + 距离标注）。"""
        if len(payload) >= 5 and isinstance(payload[4], dict):
            transfers, blocks, family_size, family, stats = payload[:5]
        else:
            transfers, blocks, family_size, family = payload[:4]
            stats = new_diagnostics("eccetera", CONVENTION_DONOR_ITSELF)
            stats["blocks_seen"] = blocks
            stats["sample_denominator"] = blocks
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
        eccetera_files: List[str],
        output_path: Optional[str] = None,
        max_workers: Optional[int] = None,
        parallel: str = "process",
    ) -> ConstraintSet:
        """把 ecceTERA recPhyloXML 文件列表转换为 ``ConstraintSet``。

        Args:
            species_tree: 物种树。
            eccetera_files: ecceTERA recPhyloXML 输出文件路径列表。
            output_path: 可选，若提供则把约束以 CSV 格式写出。
            max_workers: 并行进程/线程数；``None`` 时取默认。
            parallel: 并行模式，``"process"``（默认）或 ``"thread"``。

        Returns:
            统一 ``ConstraintSet``，每条约束 ``metadata`` 含
            ``family`` / ``support`` / ``distance``；``diagnostics`` 含口径、
            数字 ID 反解数、命中率与各类丢弃计数。

        Raises:
            UpstreamParseError: 任一输入 XML 解析失败。
            LabelMismatchError: 端点命中率为 0。
        """
        if not eccetera_files:
            cset = ConstraintSet(diagnostics=new_diagnostics("eccetera", CONVENTION_DONOR_ITSELF))
            finalize_adapter_diagnostics(
                cset, self.min_endpoint_hit_rate, set(), "eccetera", quiet=True, files=0
            )
            if output_path is not None:
                self._write_constraints_file(cset, output_path)
            return cset

        valid_labels = _build_valid_labels(species_tree)
        numeric_map = _numeric_species_id_map(species_tree) if self.resolve_numeric_ids else None
        identity = _species_cache_identity(species_tree)
        worker_args = [
            (
                rf,
                self.min_family_size,
                self.cache_dir,
                self.min_support,
                valid_labels,
                identity,
                numeric_map,
            )
            for rf in eccetera_files
        ]

        payloads = map_payloads(
            _process_eccetera_worker, worker_args, max_workers=max_workers, parallel=parallel
        )

        cset = ConstraintSet(diagnostics=new_diagnostics("eccetera", CONVENTION_DONOR_ITSELF))
        for payload in payloads:
            if payload is not None:
                self._merge_payload(cset, species_tree, payload)

        # 先落盘再跑守卫：0 命中会抛错，但 -o 产物仍须写出以便用户核对命名
        if output_path is not None:
            self._write_constraints_file(cset, output_path)
        finalize_adapter_diagnostics(
            cset,
            self.min_endpoint_hit_rate,
            valid_labels,
            "eccetera",
            quiet=self.quiet,
            files=len(eccetera_files),
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


def convert_from_eccetera(
    species_tree: Tree,
    eccetera_files: List[str],
    min_support: float = MINIMUM_SUPPORT_WITHIN_A_FAMILY,
    min_family_size: int = MINIMUM_FAMILY_SIZE,
    cache_dir: Optional[str] = None,
    output_path: Optional[str] = None,
    max_workers: Optional[int] = None,
    parallel: str = "process",
    min_endpoint_hit_rate: float = DEFAULT_MIN_ENDPOINT_HIT_RATE,
    resolve_numeric_ids: bool = True,
    quiet: bool = False,
) -> ConstraintSet:
    """便捷函数：把 ecceTERA recPhyloXML 文件列表转换为 ``ConstraintSet``。

    Args:
        species_tree: 物种树。
        eccetera_files: ecceTERA recPhyloXML 输出文件路径列表。
        min_support: 单家族内最小支持度阈值。
        min_family_size: 最小基因家族规模。
        cache_dir: 文件级缓存目录（断点续传）。
        output_path: 可选约束写出路径。
        max_workers: 并行解析进程/线程数。
        parallel: 并行模式，``"process"``（默认）或 ``"thread"``。
        min_endpoint_hit_rate: 端点命中率告警阈值（默认 0.5；0 命中必报错）。
        resolve_numeric_ids: 是否把数字物种节点 ID 反解回物种树标签（默认开）。
        quiet: 不打印 stderr 摘要。

    Returns:
        统一 ``ConstraintSet``。
    """
    adapter = EcceTERAAdapter(
        min_support=min_support,
        min_family_size=min_family_size,
        cache_dir=cache_dir,
        min_endpoint_hit_rate=min_endpoint_hit_rate,
        resolve_numeric_ids=resolve_numeric_ids,
        quiet=quiet,
    )
    return adapter.convert(
        species_tree,
        eccetera_files,
        output_path=output_path,
        max_workers=max_workers,
        parallel=parallel,
    )
