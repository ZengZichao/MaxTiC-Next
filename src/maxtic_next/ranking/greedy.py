"""贪婪启发式 ``order_from_graph``。

等价原版 ``MaxTiC.order_from_graph(graph)``：在给定的有向无环图（树种系边 + 被选中的约束边）
上做拓扑式遍历，依次取"当前未标记起点 → 沿边尽可能向下"的节点，逆序后得到内部节点排序。

注意：结果**依赖 ``graph`` 的插入顺序**（原版 ``dict`` 保序 + Py3.7+ 保序），
因此本实现严格保持原版的建图顺序（内部节点先行、叶子随后；约束边按权重降序追加）。
"""

from typing import Dict, List


def order_from_graph(graph: Dict[str, List[str]], leaves: List[str]) -> List[str]:
    """等价原版 ``order_from_graph``，返回内部节点排序（bootstrap 标签列表）。"""
    result: List[str] = []
    marques: Dict[str, int] = {}
    for l in leaves:
        marques[l] = 0
    # 防御性兑底：迭代上限 = 节点数 × (节点数 + 1)，
    # 防止自环边等意外情况导致死循环。正常情况下不会触及此上限。
    max_iter = (len(graph) + 1) * (len(graph) + 1)
    while len(result) < len(graph.keys()) - len(leaves):
        keys = list(graph.keys())
        i = 0
        while keys[i] in marques:
            i = i + 1
        current = keys[i]
        suivant = True
        iter_count = 0
        while suivant:
            suivant = False
            for v in graph[current]:
                if v not in marques and v != current:  # v != current 防止自环边挂死
                    suivant = True
                    current = v
            iter_count += 1
            if iter_count > max_iter:
                break
        result.append(current)
        marques[current] = 0
    result.reverse()
    return result
