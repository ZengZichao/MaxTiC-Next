"""``path`` 可达性测试的单元测试。"""

from maxtic_next.ranking.path import path


def test_path_simple_chain():
    graph = {"a": ["b"], "b": ["c"], "c": [], "L": []}
    assert path(graph, "a", "c") is True
    assert path(graph, "a", "b") is True
    assert path(graph, "c", "a") is False
    assert path(graph, "b", "a") is False


def test_path_self():
    graph = {"a": ["b"], "b": []}
    assert path(graph, "a", "a") is True
    assert path(graph, "b", "b") is True


def test_path_diamond():
    graph = {"a": ["b", "c"], "b": ["d"], "c": ["d"], "d": []}
    assert path(graph, "a", "d") is True
    assert path(graph, "b", "c") is False
