"""约束过滤（按距离 / 按阈值）的单元测试。"""

from maxtic_next.constraints.constraint import Constraint, ConstraintSet
from maxtic_next.constraints.filter import filter_by_threshold


def _cset():
    return ConstraintSet(
        [
            Constraint("61", "62", 10.0, {"distance": 5.0}),
            Constraint("61", "67", 20.0, {"distance": 1.0}),
            Constraint("62", "45", 30.0, {"distance": None}),  # 无距离列
        ]
    )


def test_filter_by_distance_keeps_none_and_above():
    cset = _cset()
    # d=3：distance>3 保留（5.0 的 61->62 保留；1.0 的 61->67 丢弃；无距离保留）
    cset.filter_by_distance(3)
    keys = {c.to_edge_key() for c in cset.constraints}
    assert "61,62" in keys
    assert "62,45" in keys
    assert "61,67" not in keys


def test_filter_by_distance_zero_keeps_all():
    cset = _cset()
    cset.filter_by_distance(0)
    # 所有 distance>0 或无距离都保留
    assert len(cset) == 3


def test_filter_by_threshold_removes_lowest():
    edge = {"a,b": 1.0, "c,d": 2.0, "e,f": 7.0}
    total = 10.0
    # 移除累计 < 10*0.3 = 3.0 的最低权重边：先移除 1.0，再移除 2.0 -> 累计 3.0，不 < 3.0 停止
    result = filter_by_threshold(dict(edge), total, 0.3)
    assert "a,b" not in result
    assert "c,d" not in result
    assert "e,f" in result


def test_filter_by_threshold_zero_noop():
    edge = {"a,b": 1.0, "c,d": 2.0}
    result = filter_by_threshold(dict(edge), 3.0, 0.0)
    assert len(result) == 2
