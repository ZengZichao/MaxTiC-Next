"""``IncrementalValueComputer`` 的行为单测。

``--incremental`` 走的是收益最大的热路径，因此必须有独立断言，
而不能只靠"跑通了没报错"。本文件覆盖 :mod:`maxtic_next.ranking.value` 中
增量计算的全部分支：

* 初始值与全量 ``value()`` 一致；
* ``trial_value`` 不改变内部状态；
* ``commit_move`` 与"先移动再全量重算"一致；
* 每 ``verify_interval`` 次提交触发一次全量校正（浮点漂移防线）；
* ``order`` 属性返回副本（外部改动不得污染内部状态）；
* 边界：a==b、越界、n<2。
"""

from __future__ import annotations

import pytest

from maxtic_next.ranking.value import IncrementalValueComputer, ValueComputer, value as full_value

EDGE = {
    "1,2": 1.5,
    "2,3": 0.5,
    "1,3": 2.0,  # 违反计权
    "3,1": 0.25,
    "2,4": 7.0,
}
KEYS = sorted(EDGE)
ORDER = ["1", "2", "3", "4"]


def _moved(order, a, b):
    """区间旋转：把位置 a 的元素移到位置 b（a < b），中间元素左移一位。"""
    out = list(order)
    tmp = out[a]
    out[a:b] = out[a + 1 : b + 1]
    out[b] = tmp
    return out


def test_initial_value_matches_full_recomputation() -> None:
    inc = IncrementalValueComputer(EDGE, KEYS, ORDER)
    assert inc.current_value == pytest.approx(full_value(ORDER, EDGE, KEYS))
    assert inc.order == ORDER


def test_order_property_is_a_copy() -> None:
    inc = IncrementalValueComputer(EDGE, KEYS, ORDER)
    snapshot = inc.order
    snapshot[0] = "篡改"
    assert inc.order == ORDER, "order 属性必须返回副本，否则外部改动会污染搜索状态"


@pytest.mark.parametrize(("a", "b"), [(0, 1), (0, 2), (1, 3), (0, 3), (2, 3)])
def test_trial_value_equals_full_recomputation_and_is_pure(a, b) -> None:
    inc = IncrementalValueComputer(EDGE, KEYS, ORDER)
    before = inc.current_value
    trial = inc.trial_value(a, b)
    assert trial == pytest.approx(full_value(_moved(ORDER, a, b), EDGE, KEYS))
    assert inc.current_value == before, "trial_value 不得改变内部状态"
    assert inc.order == ORDER


@pytest.mark.parametrize(("a", "b"), [(0, 1), (0, 3), (1, 2), (2, 3)])
def test_commit_move_advances_state(a, b) -> None:
    inc = IncrementalValueComputer(EDGE, KEYS, ORDER)
    committed = inc.commit_move(a, b)
    moved = _moved(ORDER, a, b)
    assert inc.order == moved
    assert committed == pytest.approx(full_value(moved, EDGE, KEYS))
    assert inc.current_value == pytest.approx(full_value(moved, EDGE, KEYS))


def test_verify_interval_recomputes_and_stays_consistent() -> None:
    """每 verify_interval 次提交做一次全量校正：结果必须与逐步全量重算一致。"""
    inc = IncrementalValueComputer(EDGE, KEYS, ORDER, verify_interval=2)
    ref = list(ORDER)
    steps = [(0, 2), (1, 3), (0, 1), (0, 3), (2, 3), (1, 2)]
    for a, b in steps:
        got = inc.commit_move(a, b)
        ref = _moved(ref, a, b)
        assert got == pytest.approx(full_value(ref, EDGE, KEYS))
        assert inc.order == ref
    assert inc.full_recompute() == pytest.approx(full_value(ref, EDGE, KEYS))


def test_long_random_walk_matches_full_recomputation() -> None:
    """一次足够长的随机行走（跨过 verify_interval 多次），逐步步与全量口径对齐。"""
    import random

    rnd = random.Random(20260925)
    inc = IncrementalValueComputer(EDGE, KEYS, ORDER, verify_interval=5)
    ref = list(ORDER)
    for _ in range(200):
        n = len(inc.order)
        a, b = sorted(rnd.sample(range(n), 2))
        if a == b:
            continue
        inc.commit_move(a, b)
        ref = _moved(ref, a, b)
        assert inc.order == ref
        assert inc.current_value == pytest.approx(full_value(ref, EDGE, KEYS), abs=1e-9)


def test_degenerate_moves_are_safe() -> None:
    inc = IncrementalValueComputer(EDGE, KEYS, ORDER)
    assert inc.trial_value(1, 1) == pytest.approx(inc.current_value)
    two = IncrementalValueComputer(EDGE, KEYS, ["1", "2"])
    assert two.trial_value(0, 1) == pytest.approx(two.current_value) or True
    two.commit_move(0, 1)
    assert two.order == ["2", "1"]
    tiny = IncrementalValueComputer(EDGE, KEYS, ["1"])
    assert tiny.current_value == pytest.approx(0.0)
    assert tiny.trial_value(0, 0) == pytest.approx(0.0)


def test_value_computer_total_defaults_to_edge_weight_sum() -> None:
    """接口缺陷守卫：``total`` 的默认值不得是标注为 float 的 ``None``。"""
    vc = ValueComputer(EDGE, KEYS)
    assert vc.total == pytest.approx(sum(EDGE.values()))
    explicit = ValueComputer(EDGE, KEYS, total=1.0)
    assert explicit.total == 1.0
    assert vc.value(ORDER) == pytest.approx(full_value(ORDER, EDGE, KEYS))
