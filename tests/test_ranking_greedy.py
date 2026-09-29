"""``order_from_graph`` 贪婪启发式测试的单元测试。"""

from maxtic_next.ranking.greedy import order_from_graph


def test_order_from_graph_simple_chain():
    graph = {"a": ["b"], "b": ["c"], "c": [], "L": []}
    leaves = ["L"]
    order = order_from_graph(graph, leaves)
    assert order == ["a", "b", "c"]


def test_order_from_graph_diamond():
    graph = {"a": ["b", "c"], "b": ["d"], "c": ["d"], "d": [], "L": []}
    leaves = ["L"]
    order = order_from_graph(graph, leaves)
    # 合法拓扑序：a 在 b、c 之前；b、c 在 d 之前
    idx = {n: i for i, n in enumerate(order)}
    assert idx["a"] < idx["b"]
    assert idx["a"] < idx["c"]
    assert idx["b"] < idx["d"]
    assert idx["c"] < idx["d"]
    assert len(order) == 4


def test_order_from_graph_only_internal():
    graph = {"a": ["b"], "b": ["L1", "L2"], "L1": [], "L2": []}
    leaves = ["L1", "L2"]
    order = order_from_graph(graph, leaves)
    assert set(order) == {"a", "b"}
    assert "L1" not in order and "L2" not in order
