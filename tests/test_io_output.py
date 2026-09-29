"""输出层（三文件命名 + 摘要）测试的单元测试。"""

import os

from maxtic_next.config import (
    SUFFIX_CONFLICTING,
    SUFFIX_FILTERED,
    SUFFIX_PARTIAL_ORDER,
)
from maxtic_next.io.output import format_summary, write_three_files
from maxtic_next.ranking.ranker import Result


def test_three_file_suffixes():
    inf = "X" + SUFFIX_FILTERED
    con = "X" + SUFFIX_CONFLICTING
    par = "X" + SUFFIX_PARTIAL_ORDER
    assert inf.endswith("_MT_output_filtered_list_of_weighted_informative_constraints")
    assert con.endswith("_MT_output_list_of_constraints_conflicting_with_best_order")
    assert par.endswith("_MT_output_partial_order")


def test_write_three_files(tmp_path):
    prefix = str(tmp_path / "CONS")
    informative, conflicting, partial = write_three_files(
        prefix, ["61,62 15.04"], ["61,62 15.04"], ["61 62 15.04 black"]
    )
    assert os.path.exists(informative)
    assert os.path.exists(conflicting)
    assert os.path.exists(partial)
    with open(informative) as f:
        assert f.read() == "61,62 15.04\n"


def test_default_style_is_short(tmp_path):
    """默认（不传 output_style）为简洁短名：<prefix>.mt.*.tsv。"""
    prefix = str(tmp_path / "CONS")
    inf, con, par = write_three_files(prefix, ["a"], ["b"], ["c"])
    assert inf.endswith(".mt.informative.tsv")
    assert con.endswith(".mt.conflicts.tsv")
    assert par.endswith(".mt.partial_order.tsv")


def test_legacy_style_matches_original(tmp_path):
    """``output_style="legacy"`` 产出与原版逐字节等价的长名。"""
    prefix = str(tmp_path / "CONS")
    inf, con, par = write_three_files(prefix, ["a"], ["b"], ["c"], output_style="legacy")
    assert inf.endswith("_MT_output_filtered_list_of_weighted_informative_constraints")
    assert con.endswith("_MT_output_list_of_constraints_conflicting_with_best_order")
    assert par.endswith("_MT_output_partial_order")


def test_short_and_legacy_same_content(tmp_path):
    """两种命名风格下三文件内容逐字节一致（仅文件名不同）。"""
    lines = (["61,62 15.04"], ["62,65 3.0"], ["61 62 15.04 black"])
    s = write_three_files(str(tmp_path / "S"), *lines, output_style="short")
    lg = write_three_files(str(tmp_path / "L"), *lines, output_style="legacy")
    for a, b in zip(s, lg):
        assert open(a).read() == open(b).read()


def test_format_summary_contains_key_sections():
    r = Result(
        constraint_file="CONS.tsv",
        internal_node_count=13,
        total_weight=100.0,
        uninformative={
            "total": 10.0,
            "to_desc": 4.0,
            "to_leaf": 2.0,
            "to_anc": 3.0,
            "to_itself": 1.0,
        },
        from_leaf=0.0,
        trivial_conflict=5.0,
        values={"input": 50.0, "greedy": 20.0, "mixing": 25.0},
        best_source="greedy heuristic",
        ranked_newick="((...))69;",
        similarity_to_input=0.5,
        conflict_with_input=3.0,
        partial_total=10.0,
    )
    out = format_summary(r)
    assert "tree with  13 internal nodes" in out  # 原版 print 产生两个空格

    assert "total weight of constraints from transfers" in out
    assert "best order is the greedy heuristic" in out
    assert "Similarity of the best order compared with the input order 0.5" in out
    # ：uninformative 百分比分母为 total_weight，不再双重计权
    # 10.0 / 100.0 -> 10.0%（旧实现给 9%）
    assert "to itself 10.0%)" in out


def test_format_summary_local_search_lines():
    """：局部搜索后打印原版两行，且数字取自回写后的 values。"""
    r = Result(
        constraint_file="CONS.tsv",
        internal_node_count=13,
        total_weight=100.0,
        uninformative={
            "total": 10.0,
            "to_desc": 4.0,
            "to_leaf": 2.0,
            "to_anc": 3.0,
            "to_itself": 1.0,
        },
        trivial_conflict=5.0,
        values={"input": 50.0, "greedy": 20.0, "mixing": 25.0, "local_search": 15.0, "best": 15.0},
        best_source="greedy heuristic + local search",
        ranked_newick="((...))69;",
        similarity_to_input=0.5,
        conflict_with_input=3.0,
        partial_total=10.0,
    )
    out = format_summary(r)
    assert "after local search 15.0 rejected" in out
    assert "best found solution 15.0 (15.0%)" in out
    assert "best order is the greedy heuristic + local search" in out
