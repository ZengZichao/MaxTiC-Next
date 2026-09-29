"""约束双格式解析的单元测试。"""

import pytest

from maxtic_next.constraints.parsers import parse_constraints


def test_space_format_with_weight():
    cset = parse_constraints(["61 62 15.04\n", "61 67 15.48\n"])
    assert len(cset) == 2
    c = cset.constraints[0]
    assert c.donor == "61"
    assert c.receptor == "62"
    assert c.weight == 15.04
    assert c.metadata["distance"] is None


def test_space_format_default_weight():
    cset = parse_constraints(["61 62\n"])
    assert len(cset) == 1
    assert cset.constraints[0].weight == 1.0


def test_space_format_with_weight_and_distance():
    # 空格格式第 4 列为 phylogenetic distance（原版 MaxTiC.py:411）
    cset = parse_constraints(["61 62 15.04 3.5\n", "61 67 20.0 0.5\n"])
    assert len(cset) == 2
    c0 = cset.constraints[0]
    assert c0.weight == 15.04
    assert c0.metadata["distance"] == 3.5
    # distance <= 默认值 0 的约束会在 filter_by_distance 中被丢弃
    assert cset.constraints[1].metadata["distance"] == 0.5


def test_space_format_distance_absent_is_none():
    # 仅 3 列（donor receptor weight）时 distance 为 None（不过滤）
    cset = parse_constraints(["61 62 15.04\n"])
    assert cset.constraints[0].metadata["distance"] is None


def test_comma_format_drops_first_column():
    # family,donor,receptor,weight
    cset = parse_constraints(["FAM1,61,62,15.04\n"])
    assert len(cset) == 1
    c = cset.constraints[0]
    assert c.donor == "61"
    assert c.receptor == "62"
    assert c.weight == 15.04
    assert c.metadata["family"] == "FAM1"


def test_comma_three_col_family_donor_receptor_rejected():
    # family,donor,receptor（缺 weight）仍是三列逗号，必须报错。
    # 逗号格式最小合法行为 4 列：family,donor,receptor,weight。
    with pytest.raises(ValueError):
        parse_constraints(["FAM1,61,62\n"])


def test_comma_format_five_columns_distance():
    # family,donor,receptor,weight,distance（ALE 5 列）
    cset = parse_constraints(["FAM1,61,62,15.04,3.5\n"])
    c = cset.constraints[0]
    assert c.weight == 15.04
    assert c.metadata["distance"] == 3.5


def test_three_column_comma_rejected():
    # 三列逗号 = donor,receptor,weight 误用，必须报错
    with pytest.raises(ValueError):
        parse_constraints(["61,62,15.04\n"])


def test_comments_and_frq_ignored():
    cset = parse_constraints(["# comment\n", "FRQ foo bar\n", "61 62 1.0\n"])
    assert len(cset) == 1


def test_blank_lines_ignored():
    cset = parse_constraints(["\n", "   \n", "61 62 1.0\n"])
    assert len(cset) == 1


def test_two_column_comma_raises_clear_error():
    # 回归：两列逗号行（"61,62" 或 "family,donor"）应报清晰 ValueError，
    # 而非 IndexError（首列丢弃规则下结构不完整）
    with pytest.raises(ValueError, match="4 列"):
        parse_constraints(["61,62\n"])
