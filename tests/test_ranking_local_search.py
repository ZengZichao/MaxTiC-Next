"""``optimisation_locale`` 局部搜索测试的单元测试。"""

from maxtic_next.random_ import RandomWrapper
from maxtic_next.ranking.local_search import optimisation_locale
from maxtic_next.ranking.value import value


def test_local_search_zero_duration_returns_same():
    edge = {"a,b": 2.0, "b,c": 3.0}
    edge_keys = list(edge.keys())
    rng = RandomWrapper(42)
    order = ["a", "b", "c"]
    # duration=0 -> 循环不执行，返回原序副本
    result = optimisation_locale(order, edge, edge_keys, rng, 0.001, 0.0)
    assert result == ["a", "b", "c"]


def test_local_search_keeps_valid_permutation():
    edge = {"a,b": 2.0, "b,c": 3.0, "a,c": 1.0}
    edge_keys = list(edge.keys())
    rng = RandomWrapper(42)
    order = ["a", "b", "c"]
    result = optimisation_locale(order, edge, edge_keys, rng, 0.001, 0.01)
    assert sorted(result) == ["a", "b", "c"]
    assert len(result) == 3


def test_local_search_value_monotonic_non_increasing():
    """局部搜索返回的最佳序目标值不差于初始序（best 单调不增）。"""
    edge = {"a,b": 2.0, "b,c": 3.0, "a,c": 1.0}
    edge_keys = list(edge.keys())
    order = ["a", "b", "c"]
    v0 = value(order, edge, edge_keys)
    rng = RandomWrapper(42)
    result = optimisation_locale(order, edge, edge_keys, rng, 0.001, 0.02)
    assert value(result, edge, edge_keys) <= v0
    assert sorted(result) == ["a", "b", "c"]


def test_local_search_accepts_improving_swap():
    """构造一个任意区间旋转都改进初始序的实例，断言局部搜索接受了改进交换：
    返回序的目标值严格小于初始序（温度极低，改进移动以 ratio=1 被接受）。

    边 (b,a)/(c,a)/(c,b) 在初始序 [a,b,c] 中三者均被违反（value=3），其任意单区间
    旋转都使 value 降到 2 或 1，因此首个被接受的旋转必然改进；以此验证 Metropolis
    接受准则对改进移动的正确处理。
    """
    edge = {"b,a": 1.0, "c,a": 1.0, "c,b": 1.0}
    edge_keys = list(edge.keys())
    order = ["a", "b", "c"]
    v0 = value(order, edge, edge_keys)
    assert v0 == 3.0
    rng = RandomWrapper(42)
    result = optimisation_locale(order, edge, edge_keys, rng, 0.001, 0.02)
    assert value(result, edge, edge_keys) < v0
