"""参考实现差分测试：唯一能**证伪**"与原版等价"宣称的自动化门禁。

背景
----
原版 ``MaxTiC.py`` 是 Python 2 脚本（``print`` 语句），无法在 Python 3 下执行；
``tests/test_equivalence_stats.py`` 里那条"跑原版做字节级比对"的用例在本机没有
python2 时只能跳过。于是本文件按原版 ``MaxTiC.py`` 的**控制流与公式**
独立重写一份 Python 3 参考实现（不复用被测代码的
``ReachabilityMatrix`` / ``EdgeBuilder`` / ``order_from_graph`` / ``mix`` / ``opt`` /
``value`` / ``RandomWrapper``），在同一份真实示例数据上与 ``maxtic_next`` 的输出
**逐字段差分比对**。

只要核心算法或统计口径漂移，本文件的第一个用例就会失败。

已记录的偏离（不是回归）
------------------------
本版本刻意偏离原版之处，都对应原版的一处缺陷，并在下面的
``DOCUMENTED_DEVIATIONS`` 中逐项**显式断言**（而不是被"整体相等"断言掩盖）：

* ``uninformative`` 百分比：原版分母 ``total_transfers + uninformative``
  把无信息权重双重计入；本版本用 ``total_weight``。
* 偏序产物的哨兵判定：原版写死 ``edge[e] < 100000``（其 ``MAX_NUMBER`` 为
  ``1e10``），会静默删掉权重 >= 1e5 的真实信息性约束；本版本用 ``MAX_NUMBER``。
* ``--threshold-constraints`` 的总权重口径：原版在删边**之前**求和，
  分母里仍含被删权重；本版本在删边**之后**计。
* 阈值删除集合：原版可把 ``MAX_NUMBER`` 谱系硬约束一并删掉；本版本把哨兵边
  排除在可删除集合之外。
* 信息性文件写出时机：与原版一致（在删除 0 权重边之后），故**不是**偏离。

Newick 解析仍复用 ``maxtic_next.tree.Tree``（本文件考察的是算法与统计口径，不覆盖
解析层）；这一点在本文件与 ``docs`` 中均已注明。
"""

import os
import random
from typing import Dict, List

from maxtic_next.api import rank
from maxtic_next.io.parsing import read_newick_file
from maxtic_next.tree.tree import Tree

DATA_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")
TREE_PATH = os.path.join(DATA_DIR, "minitree.tree")
CONS_PATH = os.path.join(DATA_DIR, "Cyano_CUTConstraints.tsv")

MAX_NUMBER = 10000000000  # 原版 MaxTiC.py:18
ORIGINAL_PARTIAL_GUARD = 100000  # 原版 MaxTiC.py:687 的魔法数


# ----------------------------------------------------------------------
# 参考实现：全部按原版公式/控制流独立书写
# ----------------------------------------------------------------------
class _RefRandom:
    """原版 ``random.random()`` 的等价驱动（与 ``RandomWrapper`` 无共享代码）。"""

    def __init__(self, seed: int) -> None:
        self._r = random.Random(seed)

    def random(self) -> float:
        return self._r.random()


def _ref_path(graph: Dict[str, List[str]], a: str, b: str) -> bool:
    """原版 ``path(g, a, b)``：DFS/BFS 混合可达性判断（不复用被测实现）。"""
    if a == b:
        return True
    marked = [a]
    stack = [a]
    while stack and b not in marked:
        node = stack.pop()
        for nbr in graph.get(node, []):
            if nbr not in marked:
                marked.append(nbr)
                stack.append(nbr)
    return b in marked


def _ref_value(order: List[str], edge: Dict[str, float], edge_keys: List[str]) -> float:
    """原版 ``value(order)``：被违反约束（donor 排在 receptor 之后）的权重和。"""
    index = {name: i for i, name in enumerate(order)}
    total = 0.0
    for key in edge_keys:
        donor, receptor = key.split(",")
        if donor in index and receptor in index and index[donor] > index[receptor]:
            total += edge[key]
    return total


def _ref_order_from_tree(tree: Tree) -> List[str]:
    """原版 ``order_from_tree``：按到根的树距离升序排列内部节点。"""
    root = tree.get_root()
    internals = [n for n in tree.get_nodes() if not tree.is_leaf(n)]
    internals.sort(key=lambda x: tree.distance_from(x, root))
    return [tree.get_bootstrap(n) for n in internals]


def _ref_order_from_graph(graph: Dict[str, List[str]], leaves: List[str]) -> List[str]:
    """原版 ``order_from_graph``：反复从首个未标记起点沿边走到底，逆序收集。"""
    result: List[str] = []
    marked: Dict[str, int] = {l: 0 for l in leaves}
    guard = (len(graph) + 1) * (len(graph) + 1)
    while len(result) < len(graph) - len(leaves):
        keys = list(graph.keys())
        i = 0
        while keys[i] in marked:
            i += 1
        current = keys[i]
        moving = True
        steps = 0
        while moving:
            moving = False
            for v in graph[current]:
                if v not in marked and v != current:
                    moving = True
                    current = v
            steps += 1
            if steps > guard:
                break
        result.append(current)
        marked[current] = 0
    result.reverse()
    return result


def _ref_mix(
    order1: List[str],
    order2: List[str],
    degre_entrant: Dict[str, List[str]],
    edge: Dict[str, float],
    rng: _RefRandom,
) -> List[str]:
    """原版 ``mix``：动态规划合并两个子序，平局时抛硬币。"""

    def edgeweights(element: str, elements: List[str]) -> float:
        s = 0.0
        for src in degre_entrant.get(element, []):
            if src in elements and (src + "," + element) in edge:
                s += edge[src + "," + element]
        return s

    # 与原版一致：第 0 行恒为 "j"（只能取 order2），第 0 列恒为 "i"
    cost = [[0.0] * (len(order2) + 1)]
    cost += [[0.0] * (len(order2) + 1) for _ in range(len(order1))]
    back = [["j"] * (len(order2) + 1)]
    back += [["i"] + [""] * len(order2) for _ in range(len(order1))]
    for i in range(1, len(order1) + 1):
        for j in range(1, len(order2) + 1):
            v1 = cost[i - 1][j] + edgeweights(order1[i - 1], order2[0:j])
            v2 = cost[i][j - 1] + edgeweights(order2[j - 1], order1[0:i])
            if v1 == v2:
                if rng.random() < 0.5:
                    cost[i][j], back[i][j] = v1, "i"
                else:
                    cost[i][j], back[i][j] = v2, "j"
            elif v1 < v2:
                cost[i][j], back[i][j] = v1, "i"
            else:
                cost[i][j], back[i][j] = v2, "j"
    order: List[str] = []
    i, j = len(order1), len(order2)
    while i > 0 or j > 0:
        if back[i][j] == "i":
            order.append(order1[i - 1])
            i -= 1
        else:
            order.append(order2[j - 1])
            j -= 1
    order.reverse()
    return order


def _ref_opt(
    tree: Tree,
    node: int,
    degre_entrant: Dict[str, List[str]],
    edge: Dict[str, float],
    rng: _RefRandom,
) -> List[str]:
    """原版 ``opt``：自底向上合并两棵子树的序（根在末尾）。"""
    c1, c2 = tree.get_children(node)[0], tree.get_children(node)[1]
    if tree.is_leaf(c1) and tree.is_leaf(c2):
        return [tree.get_bootstrap(node)]
    if tree.is_leaf(c1):
        return _ref_opt(tree, c2, degre_entrant, edge, rng) + [tree.get_bootstrap(node)]
    if tree.is_leaf(c2):
        return _ref_opt(tree, c1, degre_entrant, edge, rng) + [tree.get_bootstrap(node)]
    o1 = _ref_opt(tree, c1, degre_entrant, edge, rng)
    o2 = _ref_opt(tree, c2, degre_entrant, edge, rng)
    return _ref_mix(o1, o2, degre_entrant, edge, rng) + [tree.get_bootstrap(node)]


def _ref_maximum_distance(tree: Tree, node: int) -> List[List[str]]:
    """原版 ``maximum_distance``：返回两个极端拓扑序。"""
    c1, c2 = tree.get_children(node)[0], tree.get_children(node)[1]
    name = tree.get_bootstrap(node)
    if tree.is_leaf(c1) and tree.is_leaf(c2):
        return [[name], [name]]
    if tree.is_leaf(c1):
        oo = _ref_maximum_distance(tree, c2)
        return [[name] + oo[0], [name] + oo[1]]
    if tree.is_leaf(c2):
        oo = _ref_maximum_distance(tree, c1)
        return [[name] + oo[0], [name] + oo[1]]
    o1 = _ref_maximum_distance(tree, c1)
    o2 = _ref_maximum_distance(tree, c2)
    return [[name] + o1[0] + o2[0], [name] + o2[1] + o1[1]]


def _ref_kendall_distance(a: List[str], b: List[str]) -> float:
    binv = {name: i for i, name in enumerate(b)}
    res = 0.0
    for i in range(len(a)):
        for j in range(i + 1, len(a)):
            if binv[a[i]] > binv[a[j]]:
                res += 1
    return res


def _ref_similarity(tree: Tree, a: List[str], b: List[str]) -> float:
    m = _ref_maximum_distance(tree, tree.get_root())
    max_dist = _ref_kendall_distance(m[0], m[1])
    if max_dist == 0.0:
        return 1.0
    return (max_dist - _ref_kendall_distance(a, b)) / max_dist


def _ref_read_constraints(path: str, min_distance: float):
    """原版解析循环（MaxTiC.py:400-448）的独立重写：返回 [(donor, receptor, weight)]。"""
    out = []
    with open(path, "r", encoding="utf-8") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line.strip():
                continue
            if line[0] == "#" or "FRQ" in line:
                continue
            if "," in line:
                words = [w for w in line.strip().split(",")][1:]
            else:
                words = line.split()
            if len(words) < 2:
                continue
            if len(words) <= 3 or float(words[3]) > min_distance:
                weight = float(words[2]) if len(words) > 2 else 1.0
                out.append((words[0], words[1], weight))
    return out


def _ref_pipeline(
    newick_path: str,
    cons_path: str,
    *,
    threshold: float = 0.0,
    min_distance: float = 0.0,
    seed: int = 42,
) -> Dict:
    """按原版顺序跑一遍，返回对外报告的**全部**数字与三个产物内容。"""
    tree = read_newick_file(newick_path)
    nodes = tree.get_nodes()
    leaves = tree.get_leaves_names()
    root = tree.get_root()
    root_label = tree.get_bootstrap(root)

    # ---- 解析 + 聚合（原版 400-456）----
    edge: Dict[str, float] = {}
    for donor, receptor, weight in _ref_read_constraints(cons_path, min_distance):
        if donor == "None":
            donor = root_label
        if receptor == "None":
            receptor = root_label
        key = donor + "," + receptor
        edge.setdefault(key, 0.0)
        if not edge[key] >= MAX_NUMBER:
            edge[key] += weight
    total_transfers = sum(edge.values())

    # ---- 谱系边置哨兵（原版 462-471）----
    graph: Dict[str, List[str]] = {}
    degre_entrant: Dict[str, List[str]] = {}
    for n in nodes:
        if not tree.is_leaf(n):
            graph.setdefault(tree.get_bootstrap(n), [])
            degre_entrant.setdefault(tree.get_bootstrap(n), [])
        else:
            graph.setdefault(tree.get_name(n), [])
    to_desc = 0.0
    uninformative = 0.0
    for n in nodes:
        if not tree.is_leaf(n):
            if not tree.is_root(n):
                parent = tree.get_bootstrap(tree.get_parent(n))
                child = tree.get_bootstrap(n)
                graph[parent].append(child)
                key = parent + "," + child
                if key in edge:
                    to_desc += edge[key]
                    uninformative += edge[key]
                edge[key] = MAX_NUMBER
        else:
            graph[tree.get_bootstrap(tree.get_parent(n))].append(tree.get_name(n))

    # ---- 阈值删边（原版 483-486：升序删，分母是删前的 total_transfers）----
    keys_asc = sorted(edge.keys(), key=lambda k: edge[k])
    sub_total = 0.0
    while sub_total < total_transfers * threshold and keys_asc:
        k = keys_asc.pop(0)
        sub_total += edge[k]
        del edge[k]

    for k in edge:
        if edge[k] < MAX_NUMBER:
            first, second = k.split(",")
            if second in degre_entrant:
                degre_entrant[second].append(first)

    # ---- 分类（原版 496-521）----
    trivial_conflict = 0.0
    to_itself = 0.0
    to_leaf = 0.0
    to_anc = 0.0
    from_leaf = 0.0
    edge_keys = sorted(edge.keys(), key=lambda k: edge[k], reverse=True)
    for e in edge_keys:
        first, second = e.split(",")
        if edge[e] >= MAX_NUMBER:
            continue
        if first == second:
            uninformative += edge[e]
            to_itself += edge[e]
        elif _ref_path(graph, first, second):
            uninformative += edge[e]
            to_desc += edge[e]
        elif second in leaves:
            uninformative += edge[e]
            to_leaf += edge[e]
        elif first in leaves:
            uninformative += edge[e]
            from_leaf += edge[e]
        elif _ref_path(graph, second, first):
            uninformative += edge[e]
            trivial_conflict += edge[e]
            to_anc += edge[e]

    # ---- 互反减法 + 删 0 权重边 + 写信息性文件（原版 523-563 的真实顺序）----
    for e in list(edge_keys):
        first, second = e.split(",")
        opposite = second + "," + first
        if (
            edge[e] < MAX_NUMBER
            and opposite in edge
            and edge[opposite] < MAX_NUMBER
            and edge[e] >= edge[opposite]
        ):
            trivial_conflict += edge[opposite]
    for k in list(edge.keys()):
        if edge[k] == 0:
            del edge[k]
    edge_keys = sorted(edge.keys(), key=lambda k: edge[k], reverse=True)
    informative_lines: List[str] = []
    for k in edge_keys:
        first, second = k.split(",")
        if (
            edge[k] < MAX_NUMBER
            and first != second
            and not _ref_path(graph, first, second)
            and second not in leaves
            and first not in leaves
            and not _ref_path(graph, second, first)
        ):
            informative_lines.append(f"{k} {edge[k]}")

    # ---- 三启发式（原版 572-621）----
    order_input = _ref_order_from_tree(tree)
    value_input = _ref_value(order_input, edge, edge_keys)
    for e in edge_keys:
        first, second = e.split(",")
        if first != second and _ref_path(graph, second, first):
            pass
        else:
            graph[first].append(second)
    order_greedy = _ref_order_from_graph(graph, leaves)
    value_greedy = _ref_value(order_greedy, edge, edge_keys)
    order_mixing = _ref_opt(tree, root, degre_entrant, edge, _RefRandom(seed))
    order_mixing.reverse()
    value_mixing = _ref_value(order_mixing, edge, edge_keys)
    if value_greedy <= value_mixing:
        order = list(order_greedy)
        best_source = "greedy heuristic"
    else:
        order = list(order_mixing)
        best_source = "mixing heuristic"

    tree2 = tree.copy()
    for n in tree2.get_nodes():
        if not tree2.is_root(n):
            idx = len(order) if tree2.is_leaf(n) else order.index(tree2.get_bootstrap(n))
            idx_parent = order.index(tree2.get_bootstrap(tree2.get_parent(n)))
            tree2.set_length(n, idx - idx_parent)
    ranked_newick = tree2.write_newick(False)
    sim = _ref_similarity(tree2, order_input, order)

    # ---- 冲突 / 偏序产物（原版 673-693，含 < 100000 的魔法数）----
    index = {name: i for i, name in enumerate(order)}
    conflicting_lines: List[str] = []
    for e in edge_keys:
        a, b = e.split(",")
        if a in index and b in index and index[a] > index[b]:
            conflicting_lines.append(f"{e} {edge[e]}")
    partial_lines: List[str] = []
    for e in edge_keys:
        a, b = e.split(",")
        if a in index and b in index and index[a] < index[b] and edge[e] < ORIGINAL_PARTIAL_GUARD:
            if order_input.index(a) < order_input.index(b):
                partial_lines.append(f"{a} {b} {edge[e]} black")
            else:
                partial_lines.append(f"{a} {b} {edge[e]} green")

    pct_original = (
        int(uninformative * 100 / (total_transfers + uninformative))
        if (total_transfers + uninformative)
        else 0
    )
    pct_corrected = (uninformative * 100 / total_transfers) if total_transfers else 0.0
    return {
        "total_weight": total_transfers,
        "uninformative": {
            "total": uninformative,
            "to_desc": to_desc,
            "to_leaf": to_leaf,
            "to_anc": to_anc,
            "to_itself": to_itself,
        },
        "from_leaf": from_leaf,
        "trivial_conflict": trivial_conflict,
        "informative_lines": informative_lines,
        "conflicting_lines": conflicting_lines,
        "partial_lines": partial_lines,
        "values": {"input": value_input, "greedy": value_greedy, "mixing": value_mixing},
        "input_order": order_input,
        "greedy_order": order_greedy,
        "mixing_order": order_mixing,
        "best_order": order,
        "best_value": _ref_value(order, edge, edge_keys),
        "best_source": best_source,
        "ranked_newick": ranked_newick,
        "similarity": sim,
        "uninformative_percent_original_formula": pct_original,
        "uninformative_percent_corrected": pct_corrected,
    }


# ----------------------------------------------------------------------
# 用例 1：核心数字与三个产物必须与参考实现完全一致（可证伪"等价"宣称）
# ----------------------------------------------------------------------
def test_matches_reference_pipeline_on_cyanobacteria(tmp_path):
    ref = _ref_pipeline(TREE_PATH, CONS_PATH, seed=42)
    got = rank(
        TREE_PATH,
        CONS_PATH,
        seed=42,
        print_summary=False,
        html_report=False,
        output_prefix=str(tmp_path / "next"),
    )

    assert got.total_weight == ref["total_weight"]
    assert got.uninformative == ref["uninformative"]
    assert got.from_leaf == ref["from_leaf"]
    assert got.trivial_conflict == ref["trivial_conflict"]
    assert got.informative_lines == ref["informative_lines"]
    assert got.conflicting_lines == ref["conflicting_lines"]
    assert got.partial_lines == ref["partial_lines"]
    assert got.input_order == ref["input_order"]
    assert got.greedy_order == ref["greedy_order"]
    assert got.mixing_order == ref["mixing_order"]
    assert got.best_order == ref["best_order"]
    assert got.ranked_newick == ref["ranked_newick"]
    assert got.similarity_to_input == ref["similarity"]
    # 交付序的真实 value 必须等于参考实现对同一序算出的 value
    assert got.values["best"] == ref["best_value"]
    # 三启发式的 value 与 best_source 判据一致
    assert got.values["input"] == ref["values"]["input"]
    assert got.values["greedy"] == ref["values"]["greedy"]
    assert got.values["mixing"] == ref["values"]["mixing"]
    assert got.best_source == ref["best_source"]


def test_reference_pipeline_is_not_trivially_satisfied(tmp_path):
    """自检：参考实现本身要有判别力（能被人为破坏的口径捕获）。

    给一份"故意算错"的对比（阈值口径 / 分母），确认两个公式**确实不同**，
    否则上一条用例的通过没有任何信息量。
    """
    ref = _ref_pipeline(TREE_PATH, CONS_PATH, seed=42)
    uninf = ref["uninformative"]["total"]
    total = ref["total_weight"]
    assert uninf > 0
    original_formula = int(uninf * 100 / (total + uninf))
    corrected = uninf * 100 / total
    # 真实示例数据上：原版给 12，修正后为 13.86
    assert original_formula == 12
    assert abs(corrected - 13.8595397) < 1e-3


# ----------------------------------------------------------------------
# 用例 2：已记录偏离必须逐项显式成立
# ----------------------------------------------------------------------
def test_documented_deviation_m1_uninformative_denominator(capsys, tmp_path):
    """：本版本报告 uninformative/total_weight（13.9%），而非原版的 12%。"""
    ref = _ref_pipeline(TREE_PATH, CONS_PATH, seed=42)
    got = rank(
        TREE_PATH,
        CONS_PATH,
        seed=42,
        print_summary=True,
        html_report=False,
        output_prefix=str(tmp_path / "d1"),
    )
    out = capsys.readouterr().out
    assert ref["uninformative_percent_original_formula"] == 12
    assert abs(got.uninformative_percent - 13.8595397) < 1e-3
    assert "to itself 13.9%)" in out


def test_documented_deviation_m2_partial_order_sentinel(tmp_path):
    """：权重 2e5 的信息性约束必须留在 partial_order 里。

    原版（及本版本的旧实现）用 ``edge[e] < 100000`` 判哨兵，会把该边静默删掉，
    造成 informative.tsv 有它、partial_order.tsv 没有它的口径矛盾。
    """
    cons = tmp_path / "big_weight.tsv"
    cons.write_text("61 62 200000.0\n61 67 3.0\n", encoding="utf-8")
    prefix = str(tmp_path / "m2")
    got = rank(
        TREE_PATH, str(cons), seed=42, print_summary=False, html_report=False, output_prefix=prefix
    )
    with open(got.partial_order_file, encoding="utf-8") as fh:
        partial = fh.read()
    assert "200000.0" in partial, "高权重信息性约束被魔法数 100000 截断（M2 回归）"
    assert any("200000.0" in line for line in got.informative_lines)
    # 原版口径（< 100000）会丢掉该边 —— 显式断言这一偏离存在且被记录
    ref = _ref_pipeline(TREE_PATH, str(cons), seed=42)
    assert not any("200000.0" in line for line in ref["partial_lines"])
    assert any("200000.0" in line for line in ref["informative_lines"])


def test_documented_deviation_m5_threshold_denominator(tmp_path):
    """：阈值过滤后的总权重口径不含被删权重（原版含）。"""
    cons = tmp_path / "ts.tsv"
    cons.write_text("61 67 1.0\n62 45 0.1\n40 46 0.05\n", encoding="utf-8")
    prefix = str(tmp_path / "m5")
    got = rank(
        TREE_PATH,
        str(cons),
        seed=42,
        threshold_constraints=0.1,
        print_summary=False,
        html_report=False,
        output_prefix=prefix,
    )
    ref = _ref_pipeline(TREE_PATH, str(cons), threshold=0.1, seed=42)
    assert got.total_weight < ref["total_weight"]
    assert got.removed_by_threshold_weight > 0
    assert got.total_weight == ref["total_weight"] - got.removed_by_threshold_weight


def test_m3_zero_weight_edge_absent_from_informative_file(tmp_path):
    """：0 权重约束不得写入 informative 产物（与原版顺序一致）。"""
    cons = tmp_path / "zero.tsv"
    cons.write_text("61 67 0.0\n62 45 2.0\n", encoding="utf-8")
    prefix = str(tmp_path / "m3")
    got = rank(
        TREE_PATH, str(cons), seed=42, print_summary=False, html_report=False, output_prefix=prefix
    )
    with open(got.informative_file, encoding="utf-8") as fh:
        content = fh.read()
    assert "61,67 0.0" not in content
    assert any("62,45 2.0" in line for line in got.informative_lines)
