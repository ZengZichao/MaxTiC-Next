"""共享的物种树辅助函数（各调和适配器共用）。

把 ``_build_valid_labels`` / ``_find_label_id`` / ``_distance_from`` 从各适配器
（ale / ranger_dtl / eccetera / artra / alerax）中抽出，避免多处重复实现导致
的语义漂移。标签类助手（``build_valid_labels`` / ``find_label_id`` /
``numeric_species_id_map`` / ``resolve_species_label``）均在**主进程**中调用
（``convert`` 建 ``valid_labels``、建数字 ID 映射），``distance_from`` 亦在主进程的
``_merge_payload`` 中调用，不进入 ``ProcessPoolExecutor`` worker，故无需 pickle。

距离语义：等价于原版把物种树所有枝长设为 1.0 后的拓扑距离
（``distance_from``），供 ``--min-transfer-distance`` 按 phylogenetic distance 过滤。
：**查不到标签时返回 ``None``**（而非 0.0），否则"未知"与"真距离 0"
（供体=受体）不可区分，且该假 0 会被 ``--min-transfer-distance`` 当成合法距离
而把约束静默删除。

缓存键语义：适配器文件级缓存的载荷依赖物种树的**拓扑**（parent_map /
端点解析 / 距离），故键必须含拓扑指纹，而非只含标签集合摘要；见
``species_topology_digest`` 与 ``species_cache_identity``。
"""

import hashlib
from typing import Dict, List, Optional, Set, Tuple

from maxtic_next.constraints.adapters._diagnostics import CACHE_KEY_VERSION
from maxtic_next.tree.tree import Tree


def build_valid_labels(species_tree: Tree) -> Set[str]:
    """构建物种树所有合法标签集合（叶子名 + 内部节点 bootstrap）。

    供适配器判断 donor / receptor 是否为物种树中存在的标签（端点合法性），
    避免臆造或跨命名空间的错误约束。缺失内部标签（空字符串）不进入集合，
    防止空端点约束通过校验。
    """
    labels: Set[str] = set()
    for nid in species_tree.get_nodes():
        if species_tree.is_leaf(nid):
            name = species_tree.get_name(nid)
            if name:
                labels.add(name)
        else:
            bootstrap = species_tree.get_bootstrap(nid)
            if bootstrap:
                labels.add(bootstrap)
    return labels


def find_label_id(species_tree: Tree, label: str) -> Optional[int]:
    """按标签（叶子名或内部 bootstrap）在物种树中查找节点 id。

    优先使用 ``label_to_id`` 字典（O(1)），回退到线性搜索。
    """
    if label in species_tree.label_to_id:
        return species_tree.label_to_id[label]
    for nid in species_tree.get_nodes():
        if species_tree.is_leaf(nid):
            if species_tree.get_name(nid) == label:
                return nid
        else:
            if species_tree.get_bootstrap(nid) == label:
                return nid
    return None


def distance_from(species_tree: Tree, x: str, y: str) -> Optional[float]:
    """物种树上 ``x`` 与 ``y`` 的拓扑距离（unit branch length）。

    等价原版 ``distance_from``：先把物种树所有枝长视为 1.0，按边数计距离，
    不修改传入的物种树（避免影响后续 ``order_from_tree`` 的真实距离语义）。

    Args:
        species_tree: 物种树。
        x: 起点标签。
        y: 终点标签。

    Returns:
        拓扑距离（浮点边数）；**任一标签不存在时返回 ``None``**。
        ``None`` 表示"距离未知"，与"距离 0"（``x == y``，供体=受体自环）严格区分：
        前者在 ``ConstraintSet.filter_by_distance`` 中按"无距离列"保留，
        后者会被 ``d >= 0`` 的阈值真实过滤。
    """
    aid = find_label_id(species_tree, x)
    bid = find_label_id(species_tree, y)
    if aid is None or bid is None:
        return None
    ancestor = species_tree.last_common_ancestor(aid, bid)
    distance = 0.0
    cur = aid
    while cur != ancestor:
        distance += 1.0
        cur = species_tree.get_parent(cur)
    cur = bid
    while cur != ancestor:
        distance += 1.0
        cur = species_tree.get_parent(cur)
    return distance


def species_labels_digest(labels) -> str:
    """物种树**标签集合**的短摘要（sha256 前 16 位十六进制）。

    注意：只含标签集合，**不含父子映射**——因此它不足以做缓存标识，
    请改用 ``species_cache_identity``。
    """
    joined = "|".join(sorted(str(x) for x in labels))
    return hashlib.sha256(joined.encode("utf-8")).hexdigest()[:16]


def parent_child_pairs(species_tree: Tree) -> List[Tuple[str, str]]:
    """列出物种树全部 ``(父标签, 子标签)`` 对（叶子以叶名为标签）。

    父/子标签缺失（未标注的内部节点）时跳过该对，保证指纹只由**有名字**的
    关系构成，从而在标签集合相同、拓扑不同的两棵树之间产生不同摘要。
    """
    pairs: List[Tuple[str, str]] = []
    for nid in species_tree.get_nodes():
        if species_tree.is_leaf(nid):
            child = species_tree.get_name(nid)
        else:
            child = species_tree.get_bootstrap(nid)
        if not child:
            continue
        pid = species_tree.get_parent(nid)
        if pid is None or pid == -1:
            continue
        if species_tree.is_leaf(pid):
            parent = species_tree.get_name(pid)
        else:
            parent = species_tree.get_bootstrap(pid)
        if not parent:
            continue
        pairs.append((str(parent), str(child)))
    return pairs


def species_topology_digest(species_tree: Tree) -> str:
    """物种树**拓扑**指纹：对排序后的 ``parent_label->child_label`` 对取 sha256。

    ：两棵标签集合相同、拓扑不同的树必须有不同指纹，否则适配器文件级
    缓存会跨树复用，产出错误的 donor（ALE 尤甚）。
    """
    joined = "|".join(f"{p}->{c}" for p, c in sorted(parent_child_pairs(species_tree)))
    return hashlib.sha256(joined.encode("utf-8")).hexdigest()[:16]


def species_cache_identity(species_tree: Tree) -> str:
    """适配器文件级缓存键中的**物种树标识**：版本号 + 标签摘要 + 拓扑指纹。

    版本前缀（``v{CACHE_KEY_VERSION}``）保证任何格式不符——例如仅含标签摘要——
    的缓存条目永不被读到。
    """
    return (
        f"v{CACHE_KEY_VERSION}|{species_labels_digest(build_valid_labels(species_tree))}"
        f"|{species_topology_digest(species_tree)}"
    )


def numeric_species_id_map(species_tree: Tree) -> Dict[str, str]:
    """``{数字节点 ID 字符串: 物种树标签}``，复刻 ecceTERA 的 ID 编号规则。

    求证自参考源码 ``MySpeciesTree.cpp:112-136``（``assignPostOrderIds``，
    注释 "Assign ids using breadth-first postorder"）：先取排序后的叶子，逐层
    **自底向上广度优先**编号（0 起），同层按"首个被访问的子节点"次序入队。
    ecceTERA 在 recPhyloXML 里用 ``node->getId()``（``DTLGraph.cpp:2821``）标注
    **内部**物种节点（叶子用名字，``DTLGraph.cpp:2819``），故据此可把数字 ID
    反解回用户物种树的标签。

    未被标注（空 bootstrap）的内部节点仍占一个序号，但不进入反查表。
    """
    label_of: Dict[int, str] = {}
    for nid in species_tree.get_nodes():
        if species_tree.is_leaf(nid):
            label_of[nid] = str(species_tree.get_name(nid) or "")
        else:
            label_of[nid] = str(species_tree.get_bootstrap(nid) or "")

    level: List[int] = sorted(
        [n for n in species_tree.get_nodes() if species_tree.is_leaf(n)],
        key=lambda n: (label_of.get(n) or "", n),
    )
    ids: Dict[int, int] = {}
    counter = 0
    while level:
        nxt: List[int] = []
        for nid in level:
            if nid in ids:
                continue
            ids[nid] = counter
            counter += 1
            pid = species_tree.get_parent(nid)
            if pid is not None and pid != -1:
                nxt.append(pid)
        level = nxt
    return {str(num): label_of[nid] for nid, num in ids.items() if label_of.get(nid)}


def resolve_species_label(
    raw: str, valid_labels: Set[str], numeric_map: Optional[Dict[str, str]] = None
) -> Tuple[Optional[str], bool]:
    """把一个上游端点标签解析到物种树标签。

    解析次序（严格，避免误伤已在标签集内的数字 bootstrap）：

    1. ``raw`` 直接命中 ``valid_labels`` → ``(raw, False)``；
    2. ``raw`` 为纯数字且 ``numeric_map`` 能反解出**在标签集内**的标签 →
       ``(label, True)``（ecceTERA 用数字 ID 标注内部物种节点）；
    3. 否则 ``(None, False)``。

    Returns:
        ``(解析后的标签或 None, 是否经由数字 ID 反解)``。
    """
    if raw in valid_labels:
        return raw, False
    if numeric_map and raw.isdigit():
        label = numeric_map.get(raw)
        if label and label in valid_labels:
            return label, True
    return None, False
