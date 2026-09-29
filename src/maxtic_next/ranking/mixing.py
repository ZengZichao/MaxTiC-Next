"""混合启发式 ``opt`` / ``mix``。

等价原版 ``MaxTiC.opt(tree, root)`` 与 ``MaxTiC.mix(order1, order2)``：

* ``opt``：自底向上把每个子树的两个子序通过 ``mix`` 合并，并在末尾追加当前根
  （显式栈后序遍历，：不用递归，以免极端不平衡的物种树爆栈）；
* ``mix``：用动态规划（编辑距离式的合并代价）合并两个子序，平局时由 ``rng.random() < 0.5``
  选择 i / j（**由原版全局 ``random`` 改为受 ``RandomWrapper`` 驱动**，保证可复现）。
"""

from typing import Dict, List, Set

from maxtic_next.ranking.edge import edgeweights
from maxtic_next.random_ import RandomWrapper
from maxtic_next.tree.tree import Tree


def mix(
    order1: List[str],
    order2: List[str],
    degre_entrant: Dict[str, List[str]],
    edge: Dict[str, float],
    rng: RandomWrapper,
) -> List[str]:
    """等价原版 ``mix(order1, order2)``，合并两子序为最优顺序。"""
    order: List[str] = []
    # 标注成 List[List[float]]：表里存的是 edgeweights 的浮点累加值，
    # 用 int 表初始化会让 mypy 认定 `cost[i][j] = <float>` 是非法赋值。
    cost: List[List[float]] = [[0.0] * (len(order2) + 1)]
    for i in range(1, len(order1) + 1):
        cost.append([0.0] * (len(order2) + 1))
    back: List[List[str]] = [["j"] * (len(order2) + 1)]
    for i in range(1, len(order1) + 1):
        back.append(["i"] + [""] * (len(order2)))
    # 前缀成员判定用**增量维护的集合**替代 `order2[0:j]` 这类列表切片：
    # 原版每格重新切一片列表，使成员判定退化为 O(len) 线性扫描——这是热路径的主要开销。
    # 语义保持不变：seen2 在第 j 列恰为 order2[0:j]，seen1 在第 i 行恰为 order1[0:i]；
    # edgeweights 的累加顺序只取决于 degre_entrant，与容器类型无关 → 输出逐字节一致。
    seen1: Set[str] = set()
    for i in range(1, len(order1) + 1):
        seen1.add(order1[i - 1])
        seen2: Set[str] = set()
        for j in range(1, len(order2) + 1):
            seen2.add(order2[j - 1])
            value1 = cost[i - 1][j] + edgeweights(order1[i - 1], seen2, degre_entrant, edge)
            value2 = cost[i][j - 1] + edgeweights(order2[j - 1], seen1, degre_entrant, edge)
            if value1 == value2:
                x = rng.random()
                if x < 0.5:
                    cost[i][j] = value1
                    back[i][j] = "i"
                else:
                    cost[i][j] = value2
                    back[i][j] = "j"
            elif value1 < value2:
                cost[i][j] = value1
                back[i][j] = "i"
            else:
                cost[i][j] = value2
                back[i][j] = "j"
    i = len(order1)
    j = len(order2)
    while i > 0 or j > 0:
        if back[i][j] == "i":
            order.append(order1[i - 1])
            i = i - 1
        else:
            order.append(order2[j - 1])
            j = j - 1
    order.reverse()
    return order


def _children_or_raise(tree: Tree, node: int) -> tuple:
    """取内部节点的两个子节点；非二叉时抛出与 ``Ranker`` 同口径的错误。"""
    children = tree.get_children(node)
    if len(children) != 2:
        raise ValueError(
            f"opt/mix 只支持二叉物种树：内部节点 {tree.get_bootstrap(node)!r} 有 "
            f"{len(children)} 个子节点，多余的子树会被静默丢弃。请先解消多歧或重定根。"
        )
    return children[0], children[1]


def opt(
    tree: Tree,
    root: int,
    degre_entrant: Dict[str, List[str]],
    edge: Dict[str, float],
    rng: RandomWrapper,
) -> List[str]:
    """等价原版 ``opt(tree, root)``，返回以 ``root`` 为根的子树排序（根在末尾）。

    强前置条件：物种树必须是**二叉树**。原版与移植版都只取
    ``children[0] / children[1]``，多歧节点的多余子树会被**静默丢弃**，导致残缺排序
    漏算约束、目标值人为变优。``Ranker.run()`` 已在入口无条件校验；此处再设一道
    显式守卫，使直接调用 ``opt`` 的代码也能得到指名道姓的错误。

    实现为**显式栈的后序遍历**：原版与本函数的递归版本深度等于树高，
    在极端不平衡（caterpillar，1200 个内部节点）的物种树上会 ``RecursionError``。
    子节点按"先左后右"的顺序处理，因此 ``mix`` 的调用顺序、从而 ``rng`` 的消耗序列
    与递归实现**逐次一致**（同种子同结果），正常树上的输出与递归版完全相同。
    """
    results: Dict[int, List[str]] = {}
    stack: List[tuple] = [(root, False)]
    while stack:
        node, expanded = stack.pop()
        c1, c2 = _children_or_raise(tree, node)
        if not expanded:
            if tree.is_leaf(c1) and tree.is_leaf(c2):
                results[node] = [tree.get_bootstrap(node)]
            elif tree.is_leaf(c1):
                stack.append((node, True))
                stack.append((c2, False))
            elif tree.is_leaf(c2):
                stack.append((node, True))
                stack.append((c1, False))
            else:
                stack.append((node, True))
                stack.append((c2, False))  # c1 先出栈 => 与递归版同序
                stack.append((c1, False))
            continue
        name = tree.get_bootstrap(node)
        if tree.is_leaf(c1) and tree.is_leaf(c2):
            results[node] = [name]
        elif tree.is_leaf(c1):
            results[node] = results[c2] + [name]
        elif tree.is_leaf(c2):
            results[node] = results[c1] + [name]
        else:
            results[node] = mix(results[c1], results[c2], degre_entrant, edge, rng) + [name]
    return results[root]
