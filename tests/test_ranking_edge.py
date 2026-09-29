"""``edgeweights`` 原语与 ``EdgeBuilder``（含 RANDOM_TYPE）的单元测试。"""

import os

from maxtic_next.constraints.parsers import parse_constraints
from maxtic_next.io.parsing import read_constraints_file
from maxtic_next.ranking.edge import EdgeBuilder, edgeweights
from maxtic_next.random_ import RandomWrapper
from maxtic_next.tree.tree import Tree

DATA_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")


def _tree_and_cset():
    tree = Tree()
    with open(os.path.join(DATA_DIR, "minitree.tree")) as fh:
        tree.read_newick(fh.readline())
    cset = read_constraints_file(os.path.join(DATA_DIR, "Cyano_CUTConstraints.tsv"))
    return tree, cset


def test_edgeweights_counts_incoming():
    degre_entrant = {"b": ["a"]}
    edge = {"a,b": 5.0}
    assert edgeweights("b", ["a"], degre_entrant, edge) == 5.0


def test_edgeweights_requires_membership():
    degre_entrant = {"b": ["a"]}
    edge = {"a,b": 5.0}
    # a 不在 elements 中 -> 不计入
    assert edgeweights("b", ["x"], degre_entrant, edge) == 0.0


def test_edgeweights_multiple_incoming():
    degre_entrant = {"c": ["a", "b"]}
    edge = {"a,c": 2.0, "b,c": 3.0}
    assert edgeweights("c", ["a", "b"], degre_entrant, edge) == 5.0


def test_edgeweights_missing_key_safe():
    degre_entrant = {"c": ["a"]}
    edge = {}
    # s="a" 在 degre_entrant["c"] 中，但 edge 无 "a,c" 键 -> 安全跳过
    assert edgeweights("c", ["a"], degre_entrant, edge) == 0.0


def test_random_type_0_keeps_original_endpoints():
    """RANDOM_TYPE=0：不改方向、不随机化，边键集合等于原始 (donor,receptor)。"""
    tree, cset = _tree_and_cset()
    edge = EdgeBuilder(tree, cset, RandomWrapper(42), random_type=0).build()
    expected = {f"{c.donor},{c.receptor}" for c in cset.constraints}
    assert set(edge.keys()) == expected


def test_random_type_1_is_deterministic_and_internal():
    """RANDOM_TYPE=1：以 0.5 概率交换端点，但同 seed 可复现，端点均为内部节点标签。"""
    tree, cset = _tree_and_cset()
    e1 = EdgeBuilder(tree, cset, RandomWrapper(42), random_type=1).build()
    e2 = EdgeBuilder(tree, cset, RandomWrapper(42), random_type=1).build()
    assert e1 == e2  # 同 seed 逐字节可复现
    internal = set(tree.internal_node_labels())
    for key in e1:
        d, r = key.split(",")
        assert d in internal and r in internal


def test_random_type_2_randomizes_within_internal_nodes():
    """RANDOM_TYPE=2：两端点均随机化为内部节点（first != second），同 seed 可复现。"""
    tree, cset = _tree_and_cset()
    e1 = EdgeBuilder(tree, cset, RandomWrapper(42), random_type=2).build()
    e2 = EdgeBuilder(tree, cset, RandomWrapper(42), random_type=2).build()
    assert e1 == e2  # 同 seed 逐字节可复现
    internal = set(tree.internal_node_labels())
    for key in e1:
        d, r = key.split(",")
        assert d in internal and r in internal
        assert d != r  # 原版保证 first 与 second 不同


def test_random_type_invalid_value_is_rejected():
    """random_type 非 0/1/2 必须报错。

    旧实现把 99 / -1 / 3 静默当作 0（即"不做任何随机化"），使打错参数的用户拿到
    一份未随机化的"假对照"，并把真实数据当对照分析。
    """
    import pytest

    from maxtic_next.ranking.edge import EdgeBuilder as _EB

    tree, cset = _tree_and_cset()
    for bad in (99, -1, 3):
        with pytest.raises(ValueError, match="random_type"):
            _EB(tree, cset, RandomWrapper(42), random_type=bad).build()


def test_random_type_2_needs_two_internal_nodes():
    """random_type=2 在只有 1 个内部节点时必须报错，而不是死循环。"""
    import pytest

    from maxtic_next.constraints.constraint import Constraint, ConstraintSet
    from maxtic_next.ranking.edge import EdgeBuilder as _EB

    tree = Tree()
    tree.read_newick("(A:1.0,B:1.0)X:1.0;")
    cset = ConstraintSet([Constraint(donor="X", receptor="A", weight=1.0)])
    with pytest.raises(ValueError, match="内部节点"):
        _EB(tree, cset, RandomWrapper(1), random_type=2).build()


def test_edgebuilder_aggregates_weights():
    """同一边键的约束权重被累加。"""
    tree, _ = _tree_and_cset()
    lines = ["61 62 10.0\n", "61 62 5.0\n"]
    cset = parse_constraints(lines)
    edge = EdgeBuilder(tree, cset, RandomWrapper(0), random_type=0).build()
    assert edge["61,62"] == 15.0
