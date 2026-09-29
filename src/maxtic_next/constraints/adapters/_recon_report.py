"""共享的"调和报告"文本解析核心（RANGER-DTLx / ARTra 共用）。

经对参考软件源码求证，RANGER-DTLx 与 ARTra 输出的是**同一形状**的人类可读调和
报告，每个内部基因节点一行，转移事件形如::

     = LCA[H117, H14]: Transfer, Mapping --> H15, Edge, Parent = P, Recipient --> H125
     = LCA[H117, H14]: Replacing Transfer, Mapping --> H15, Recipient --> H125
     = LCA[H117, H14]: Additive Transfer, Mapping --> H15, Recipient --> H125

关键求证来源：

* RANGER-DTLx ``Ranger-DTLx/DTL-algorithm.h:1840``（``os << ... "Transfer, Mapping --> "
  << donor << ... ", Recipient --> " << recipient``）；
* ARTra ``ARTra-程序/output.txt``（``Replacing Transfer`` / ``Additive Transfer`` 前缀，
  其余结构一致）。

因此本模块把转移事件的语义统一抽象为：``Mapping --> <donor>`` 给出**供体物种**
（转移离开的谱系，映射所在），``Recipient --> <receptor>`` 给出**受体物种**（转移
到达的谱系），产出统一的 ``(donor, receptor)`` 有向约束（donor 应早于/不晚于
receptor）。

**供体端点层级约定（，机器可读于诊断的**
``donor_endpoint_convention = "donor_itself"``）**：本模块取"供体本身"
（``Mapping -->`` 给出的那个谱系）为约束 donor；而 ALE 路径取"供体的父"
（``parent(donnor)``，见 ``adapters/ale.py`` 与原版
``constraints_from_reconciliations.py:160-161``）。同一生物学转移事件在两条路径下
的端点相差**一个物种层级**，端点不可直接互换；``--from-auto`` 混用两种约定时
``registry`` 会给出警告。本模块保持各工具的原生约定，不做层级换算。

权重口径
-------------------------------------
**"输出块数"不是采样分母**。RANGER-DTLx 的块头（``DTL-algorithm.h:2170``）按**输入
基因树**逐棵打印，一块 = 另一棵树，不是另一次采样；真实 ARTra 输出
ARTra 源码的 ``output.txt:206`` 打印 ``Total number of optimal solutions: 24``
却只输出 1 个最优解。故本模块按下列优先级确立分母：

1. 工具**自己声明**的样本数：``Total number of optimal solutions: N``
   （``DTL-algorithm.h:1785`` 与 ``:2224``、ARTra 同名行）等；
   权重 = ``count / N``，口径 ``support_over_declared_samples``；
2. 无声明值，但各块叶子集合**完全一致**（确为同一家族的重复调和）：
   权重 = ``count / 块数``，口径 ``support_over_replicate_blocks``，并记警告；
3. 其余情形（无声明值且块数不可信 / 只有 0 块）：**权重退化为整数计数**
   （口径 ``integer_count``），并在返回的诊断 + stderr 警告里明说——
   此时 ``--min-support`` 之类别再假装自己是 [0,1] 支持度阈值。

所有静默 ``continue`` 均被计数：标签不匹配、自环
（donor == receptor）、块/转移总数，经 ``process_report_file`` 返回的 ``stats``
交给适配器并入 ``ConstraintSet.diagnostics``，并由适配器打印一行 stderr 摘要。

设计约束：本模块仅含**纯函数**（无类、无实例状态、无第三方依赖），可被
``ProcessPoolExecutor`` 子进程安全 pickle 调用；打印一律放在主进程的适配器里。
"""

import hashlib
import os
import pickle
import re
from typing import Dict, List, Optional, Set, Tuple

from maxtic_next.constraints.adapters._diagnostics import (
    CACHE_KEY_VERSION,
    CONVENTION_DONOR_ITSELF,
    WEIGHT_BLOCK_SUPPORT,
    WEIGHT_DECLARED_SUPPORT,
    WEIGHT_INTEGER_COUNT,
    new_diagnostics,
)
from maxtic_next.constraints.adapters._species import species_labels_digest
from maxtic_next.io.compression import logical_stem, read_text_lines


# 转移事件：只有转移才带 "Recipient -->"，据此触发（兼容 Transfer / Replacing
# Transfer / Additive Transfer）。donor 取 "Mapping -->" 后到逗号/行尾之间的内容；
# receptor 取 "Recipient -->" 后的第一个 token。
_MAPPING_RE = re.compile(r"Mapping\s*-->\s*([^,\n]+?)\s*(?:,|$)")
_RECIPIENT_RE = re.compile(r"Recipient\s*-->\s*(\S+)")

# 调和块分隔（RANGER-DTLx: "Reconciliation for Gene Tree 1:"；
# ARTra: "------ Reconciliation for Gene Tree 1 (rooted) ------"）。
_BLOCK_RE = re.compile(r"Reconciliation\s+for\s+Gene\s+Tree", re.IGNORECASE)

# 叶子行（家族规模估算）：形如 "H117_0: Leaf Node" 或 "H117: Leaf Node"。
_LEAF_RE = re.compile(r"^\s*(\S+)\s*:\s*Leaf\s+Node", re.IGNORECASE)

# 工具**自己声明**的样本/解数量行：分母的第一优先来源。
# 逐条求证：
#   RANGER-DTLx DTL-algorithm.h:1785 与 :2224  "Total number of optimal solutions: N"
#   ARTra       output.txt:206                 "Total number of optimal solutions: 24"
#   RANGER-DTL  "Total number of reconciliations: N" 等变体一并兼容。
_DECLARED_SAMPLE_RES = (
    re.compile(r"Total number of optimal solutions\s*:\s*(\d+)", re.IGNORECASE),
    re.compile(
        r"Number of (?:equally )?optimal (?:solutions|reconciliations)"
        r"\s*:\s*(\d+)",
        re.IGNORECASE,
    ),
    re.compile(r"Total number of (?:reconciliations|scenarios)\s*:\s*(\d+)", re.IGNORECASE),
)


def _clean_label(label: str) -> str:
    """清理物种标签：去首尾空白与可能的尾随标点。"""
    return label.strip().strip(".;,")


def parse_transfer_line(line: str) -> Tuple[str, str]:
    """从单行调和报告中抽取 ``(donor, receptor)``；非转移行返回 ``("", "")``。

    仅当行同时含 ``Mapping -->`` 与 ``Recipient -->`` 时才判定为转移事件
    （损失/复制/物种形成事件不含 ``Recipient -->``）。

    Args:
        line: 调和报告中的一行。

    Returns:
        ``(donor, receptor)`` 物种标签对；若非转移行则为 ``("", "")``。
    """
    if "Recipient" not in line or "Mapping" not in line:
        return ("", "")
    m_don = _MAPPING_RE.search(line)
    m_rec = _RECIPIENT_RE.search(line)
    if not m_don or not m_rec:
        return ("", "")
    donor = _clean_label(m_don.group(1))
    receptor = _clean_label(m_rec.group(1))
    return (donor, receptor)


def declared_sample_counts(lines: List[str]) -> List[int]:
    """收集报告里工具**自己声明**的样本/最优解数量（出现次序，去重保序）。"""
    out: List[int] = []
    for raw in lines:
        for rx in _DECLARED_SAMPLE_RES:
            m = rx.search(raw)
            if m:
                try:
                    v = int(m.group(1))
                except ValueError:  # pragma: no cover - 正则已限定数字
                    continue
                if v not in out:
                    out.append(v)
    return out


def parse_reconciliation_report_detailed(
    lines: List[str],
    valid_labels: Set[str],
    numeric_map: Optional[Dict[str, str]] = None,
) -> Tuple[Dict[str, int], int, int, Dict]:
    """解析 RANGER-DTLx / ARTra 调和报告文本，并产出完整计数诊断。

    把文件按 ``Reconciliation for Gene Tree`` 分块，块内收集转移事件计数，
    块间按边键求和。端点须能解析进 ``valid_labels``（物种树标签集合）才保留；
    每一次丢弃都被计数，不再静默 ``continue``。

    Args:
        lines: 文件的所有行。
        valid_labels: 物种树中所有合法标签（叶子名 + 内部节点 bootstrap）。
        numeric_map: 可选的 ``数字ID -> 标签`` 反查表（；
            RANGER/ARTra 原生写真实标签，故默认 ``None``）。

    Returns:
        ``(transfers, blocks, family_size, stats)``：

        * ``transfers``：``{"donor,receptor": count}``（块间求和后的计数）；
        * ``blocks``：``Reconciliation for Gene Tree`` 块数（无显式块头但有
          转移/叶子时记为 1）；
        * ``family_size``：家族规模估算（单块内叶子行数的最大值；探测失败为 0）；
        * ``stats``：诊断字典（``transfers_seen`` / ``dropped_by_label_miss`` /
          ``self_loop_transfers`` / ``declared_sample_counts`` /
          ``sample_denominator`` / ``weight_semantics`` /
          ``family_size_probe_ok`` 等）。
    """
    transfers: Dict[str, int] = {}
    blocks = 0
    started = False
    current: Dict[str, int] = {}
    leaves_in_block: Set[str] = set()
    block_leaf_sets: List[Set[str]] = []
    max_family_size = 0
    stats = new_diagnostics("recon_report", CONVENTION_DONOR_ITSELF)
    stats["weight_semantics"] = WEIGHT_INTEGER_COUNT

    def _flush() -> None:
        """把当前块的转移计数并入总计，并结算块数与家族规模。"""
        nonlocal blocks, max_family_size, current, leaves_in_block
        for k, v in current.items():
            transfers[k] = transfers.get(k, 0) + v
        blocks += 1
        block_leaf_sets.append(set(leaves_in_block))
        if len(leaves_in_block) > max_family_size:
            max_family_size = len(leaves_in_block)
        current = {}
        leaves_in_block = set()

    for raw in lines:
        line = raw.rstrip("\n")
        if _BLOCK_RE.search(line):
            # 新块开始：先结算上一块（若已有内容）
            if started:
                _flush()
            started = True
            continue

        leaf_m = _LEAF_RE.match(line)
        if leaf_m:
            leaves_in_block.add(leaf_m.group(1))
            continue

        donor, receptor = parse_transfer_line(line)
        if not (donor and receptor):
            continue
        stats["transfers_seen"] += 1
        # 端点解析（上游命名可能与物种树不一致）
        d_label, d_numeric = _resolve(donor, valid_labels, numeric_map)
        r_label, r_numeric = _resolve(receptor, valid_labels, numeric_map)
        if d_label is None or r_label is None:
            stats["dropped_by_label_miss"] += 1
            for miss in (donor if d_label is None else None, receptor if r_label is None else None):
                if miss and miss not in stats["unresolved_endpoint_samples"]:
                    if len(stats["unresolved_endpoint_samples"]) < 8:
                        stats["unresolved_endpoint_samples"].append(miss)
            continue
        stats["numeric_id_resolutions"] += int(d_numeric) + int(r_numeric)
        if d_label == r_label:
            # 自环转移（供体=受体谱系）：单独计数，仍保留给下游
            stats["self_loop_transfers"] += 1
        key = f"{d_label},{r_label}"
        current[key] = current.get(key, 0) + 1
        stats["transfers_resolved"] += 1

    # 结算最后一个块（或无显式块头的整文件）
    if started or current or leaves_in_block:
        _flush()

    stats["blocks_seen"] = blocks
    stats["family_size_probe_ok"] = max_family_size > 0

    # ---- 支持度分母----
    declared = declared_sample_counts(lines)
    stats["declared_sample_counts"] = list(declared)
    if declared:
        denominator = max(declared)
        stats["weight_semantics"] = WEIGHT_DECLARED_SUPPORT
        if len(declared) > 1:
            stats["warnings"].append(
                "调和报告声明了多个样本数 "
                f"{declared}（RANGER-DTLx/ARTra 会为每棵基因树各打印一行），"
                f"取最大值 {denominator} 作分母；若确为多家族混合报告，"
                "请按文件拆分后分别传入。"
            )
    elif blocks > 1 and _leaf_sets_consistent(block_leaf_sets):
        denominator = blocks
        stats["weight_semantics"] = WEIGHT_BLOCK_SUPPORT
        stats["warnings"].append(
            f"报告未声明样本数，按 {blocks} 个**叶子集合一致**的调和块作分母；"
            "块数只是重复调和的粗略代理，不是上游采样数。"
        )
    elif blocks >= 1:
        # 单块、或各块叶子集合互不相同（= 多棵不同基因树，绝非同一家族采样）
        denominator = None
        stats["weight_semantics"] = WEIGHT_INTEGER_COUNT
        if blocks > 1:
            stats["warnings"].append(
                f"报告含 {blocks} 个叶子集合**不同**的调和块（RANGER-DTLx/ARTra "
                "的块头是按输入基因树逐棵打印的，一块=另一棵树），故块数不可作"
                "采样分母；权重已退化为**整数计数**而非 [0,1] 支持度。请按文件"
                "拆分家族后分别传入。"
            )
        else:
            stats["warnings"].append(
                "报告未声明样本数（无 ``Total number of optimal solutions: N`` "
                "行）且只有 1 个调和块：权重按**整数计数**输出，"
                "min_support 阈值不具 [0,1] 支持度语义。"
            )
    else:
        denominator = None
        stats["weight_semantics"] = WEIGHT_INTEGER_COUNT
    stats["sample_denominator"] = denominator
    stats["min_support_is_fractional"] = stats["weight_semantics"] in (
        WEIGHT_DECLARED_SUPPORT,
        WEIGHT_BLOCK_SUPPORT,
    ) and bool(denominator)

    # 多块但叶子集合不一致的额外提示（与旧行为兼容：块间被按同一家族求和）
    if blocks > 1 and transfers and not _leaf_sets_consistent(block_leaf_sets):
        stats.setdefault("multi_tree_suspect", True)

    return transfers, blocks, max_family_size, stats


def _resolve(
    label: str, valid_labels: Set[str], numeric_map: Optional[Dict[str, str]]
) -> Tuple[Optional[str], bool]:
    """端点标签解析：直接命中 > 数字 ID 反查 > 失败。"""
    if label in valid_labels:
        return label, False
    if numeric_map and label.isdigit():
        mapped = numeric_map.get(label)
        if mapped and mapped in valid_labels:
            return mapped, True
    return None, False


def _leaf_sets_consistent(sets: List[Set[str]]) -> bool:
    """各块叶子集合是否一致（一致才可能是"同一家族的重复调和"）。"""
    if not sets:
        return False
    first = sets[0]
    return all(s == first for s in sets[1:])


def parse_reconciliation_report(
    lines: List[str],
    valid_labels: Set[str],
) -> Tuple[Dict[str, int], int, int]:
    """解析 RANGER-DTLx / ARTra 调和报告文本（兼容旧三值返回）。

    等价 ``parse_reconciliation_report_detailed(...)`` 的前三项；需要权重口径与
    丢弃计数时请直接调用 detailed 版本。

    Returns:
        ``(transfers, blocks, family_size)``。
    """
    transfers, blocks, family_size, _stats = parse_reconciliation_report_detailed(
        lines, valid_labels
    )
    return transfers, blocks, family_size


def resolve_family_size_filter(
    stats: Dict, family_size: int, min_family_size: int, probe_desc: str = "``: Leaf Node`` 行"
) -> bool:
    """``min_family_size`` 过滤判定。

    只有当**家族规模确实探测到**时才过滤。旧实现无条件
    ``if family_size <= min_family_size: transfers = {}``，于是
    ``Leaf Node`` 行探测失败（如报告改版、被截断、或 ARTra 只打印基因树而
    无节点明细）时 ``family_size=0``，导致**全部**转移被静默丢弃；
    ``--ale-min-family-size 0`` 更是必然清空。

    Args:
        stats: 该文件的诊断字典（就地更新）。
        family_size: 探测到的家族规模（0 表示未探测到）。
        min_family_size: 阈值。
        probe_desc: 探测源的描述（报告用 ``: Leaf Node`` 行，
            recPhyloXML 用 ``<leaf>`` 元素）。

    Returns:
        ``True`` 表示应丢弃该家族的转移（并已在 ``stats`` 计数）。
    """
    if stats.get("family_size_probe_ok"):
        if family_size <= min_family_size:
            stats["dropped_by_family_size"] += int(stats.get("transfers_resolved", 0))
            return True
        return False
    # 探测失败：跳过过滤 + 记警告（绝不静默清空）
    stats["family_size_probe_ok"] = False
    stats["warnings"].append(
        f"无法从报告中探测基因家族规模（未见 {probe_desc}），"
        f"已**跳过** min_family_size={min_family_size} 过滤而不是清空全部转移；"
        "请核对上游报告格式。"
    )
    return False


def process_report_file(args: Tuple) -> Tuple[Dict[str, int], int, int, str, Dict]:
    """文件级 worker（可 pickle，供 ``ProcessPoolExecutor`` 子进程调用）。

    解析单个调和报告文件（RANGER-DTLx / ARTra 共用），含文件级 pickle
    缓存（断点续传）。RANGER-DTLx 与 ARTra 适配器均委托本函数，确保两者
    解析行为逐字节一致（因两工具输出同形）。

    Args:
        args: ``(file_path, min_family_size, cache_dir, min_support,
               valid_labels)``，可选第 6/7 元素
               ``species_identity``（缓存键中的物种树标识，含拓扑指纹）与
               ``numeric_map``。

    Returns:
        ``(transfers, blocks, family_size, family, stats)``：``stats`` 为诊断字典
        （含权重口径与各类丢弃计数）。
    """
    (file_path, min_family_size, cache_dir, min_support, valid_labels) = tuple(args[:5])
    species_identity = args[5] if len(args) > 5 else None
    numeric_map = args[6] if len(args) > 6 and args[6] else None

    # ---- 文件级缓存 ----
    cache_file: Optional[str] = None
    if cache_dir:
        st = os.stat(file_path)
        # 缓存载荷依赖物种树标签集合（端点过滤）**与拓扑**（仅标签摘要
        # 会让拓扑不同的两棵树复用同一缓存）；键另带版本号，旧条目不会被读到。
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
            except Exception:
                payload = None
            if payload is not None:
                return _normalise_payload(payload)

    # ---- 读取 + 解析 ----
    # ：经 io.compression 共享读取器，.gz / 单成员 .tar.gz 上游输出可直接解析；
    # 家族名取"逻辑词干"，故 gene_1.phy.gz 与 gene_1.phy 得到同一家族标识。
    lines = read_text_lines(file_path, encoding="utf-8")
    family = logical_stem(file_path)

    transfers, blocks, family_size, stats = parse_reconciliation_report_detailed(
        lines, valid_labels, numeric_map
    )

    # 家族规模过滤（探测失败时跳过而不是清空，见 resolve_family_size_filter）
    if resolve_family_size_filter(stats, family_size, min_family_size):
        stats["families_skipped"] = 1
        transfers = {}
    else:
        stats["families_skipped"] = 0

    payload = (transfers, blocks, family_size, family, stats)

    # ---- 缓存写入 ----
    if cache_dir and cache_file is not None:
        os.makedirs(cache_dir, exist_ok=True)
        with open(cache_file, "wb") as fh:
            pickle.dump(payload, fh)

    return payload


def _normalise_payload(payload) -> Tuple[Dict[str, int], int, int, str, Dict]:
    """把缺字段的（4 元组）缓存载荷规整成 5 元组。"""
    if len(payload) >= 5:
        transfers, blocks, family_size, family, stats = payload[:5]
        if not isinstance(stats, dict):
            stats = {}
        return transfers, blocks, family_size, family, stats
    transfers, blocks, family_size, family = payload
    stats = new_diagnostics("recon_report", CONVENTION_DONOR_ITSELF)
    stats["blocks_seen"] = blocks
    stats["transfers_resolved"] = sum(transfers.values())
    stats["transfers_seen"] = stats["transfers_resolved"]
    stats["weight_semantics"] = WEIGHT_BLOCK_SUPPORT
    stats["sample_denominator"] = blocks
    stats["cache_legacy_payload"] = True
    stats["warnings"].append(
        "读到的是**缺诊断信息的缓存载荷**：其分母按块数推定，"
        "口径未必与本次一致。建议换一个 cache 目录重跑。"
    )
    return transfers, blocks, family_size, family, stats
