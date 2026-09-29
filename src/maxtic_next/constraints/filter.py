"""约束过滤工具。

包含两处与原版 ``MaxTiC.py`` 一一对应的过滤逻辑：

* ``filter_by_distance``：在 ``ConstraintSet`` 上按 phylogenetic distance 过滤
  （见 ``constraint.ConstraintSet.filter_by_distance``）；
* ``filter_by_threshold``：按权重阈值比例从聚合后的 ``edge`` 字典中移除低权重边，
  等价原版 ``while sub_total < total_transfers * THRESHOLD_CONSTRAINTS`` 循环。
"""

from typing import Dict


def filter_by_threshold(
    edge: Dict[str, float], total_transfers: float, threshold_constraints: float
) -> Dict[str, float]:
    """按权重阈值比例从 ``edge`` 中移除最低权重的边。

    等价原版逻辑：累加最小权重边直到累计值 >= ``total_transfers * threshold_constraints``，
    并将这些边从 ``edge`` 中删除（原地修改并返回）。``threshold_constraints=0`` 时不移除任何边。

    Args:
        edge: 聚合后的边权重字典（键 ``"donor,receptor"`` -> 权重）。
        total_transfers: 约束总权重（过滤前）。
        threshold_constraints: 阈值比例（0.0 ~ 1.0）。

    Returns:
        过滤后的 ``edge`` 字典（同一对象）。
    """
    edge_keys = sorted(edge.keys(), key=lambda k: edge[k])  # 升序
    sub_total = 0.0
    while sub_total < total_transfers * threshold_constraints and edge_keys:
        key = edge_keys[0]
        sub_total += edge[key]
        del edge[key]
        del edge_keys[0]
    return edge
