"""可达性判断 ``path``。

等价原版 ``MaxTiC.path(g, a, b)``：在**有向图** ``g`` 中判断是否存在从 ``a`` 到 ``b`` 的路径，
``a == b`` 时返回 ``True``。该图是算法层的增广图（树种系边 + 被选中的约束边），
键与值均为算法层字符串标识（内部节点 bootstrap 标签或叶子名）。
"""

from typing import Dict, List


def path(graph: Dict[str, List[str]], a: str, b: str) -> bool:
    """判断有向图 ``graph`` 中是否存在从 ``a`` 到 ``b`` 的路径（``a==b`` 返回 ``True``）。

    等价原版 ``path()`` 的 DFS/BFS 混合实现，逐行移植。
    """
    marques: List[str] = [a]
    pile: List[str] = [a]
    while len(pile) > 0 and (b not in marques):
        sommet = pile[-1]
        del pile[-1]
        voisins = graph[sommet]
        for v in voisins:
            if v not in marques:
                marques.append(v)
                pile.append(v)
    if b in marques:
        return True
    else:
        return False
