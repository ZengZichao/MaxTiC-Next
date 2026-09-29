"""``value`` 目标函数测试的单元测试。"""

from maxtic_next.ranking.value import value


def _edge():
    return {"a,b": 2.0, "b,c": 3.0, "a,c": 1.0}


def test_value_topological_order_zero():
    edge = _edge()
    keys = list(edge.keys())
    order = ["a", "b", "c"]
    assert value(order, edge, keys) == 0.0


def test_value_reversed_counts_all():
    edge = _edge()
    keys = list(edge.keys())
    order = ["c", "b", "a"]
    # a=2,b=1,c=0: a,b(2>1)->2; b,c(1>0)->3; a,c(2>0)->1 => 6
    assert value(order, edge, keys) == 6.0


def test_value_partial_violation():
    edge = _edge()
    keys = list(edge.keys())
    order = ["a", "c", "b"]
    # a=0,c=1,b=2: a,b(0>2)F; b,c(2>1)T->3; a,c(0>1)F => 3
    assert value(order, edge, keys) == 3.0


def test_value_ignores_leaf_endpoints():
    edge = {"a,L": 5.0, "a,b": 2.0}
    keys = list(edge.keys())
    order = ["a", "b"]
    # L 不在 order 的 index 中，故 a,L 不计入
    assert value(order, edge, keys) == 0.0
