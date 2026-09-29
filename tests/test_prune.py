"""目标类群剪裁单元测试。"""

import os

from maxtic_next.constraints.constraint import Constraint, ConstraintSet
from maxtic_next.io.parsing import read_newick_file
from maxtic_next.prune import prune_constraints

DATA_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")


def _tree():
    return read_newick_file(os.path.join(DATA_DIR, "minitree.tree"))


def _cset():
    # 61/45 都在 clade 61 内；65/56 在 clade 61 之外（65 是 61 的祖先，56 是旁支）
    return ConstraintSet(
        [
            Constraint("61", "45", 10.0),  # 两端点均在 clade 61 内
            Constraint("61", "65", 5.0),  # 跨类群（65 是 61 的祖先）
            Constraint("65", "61", 8.0),  # donor 是祖先，receptor 在类群内
            Constraint("61", "56", 3.0),  # 跨类群（56 是旁支兄弟）
        ]
    )


def test_prune_default_keeps_within_clade_only():
    """默认策略：仅保留两端点都在目标类群内的约束，丢弃一切跨类群约束。"""
    tree = _tree()
    kept = prune_constraints(tree, _cset(), "61")
    keys = {c.to_edge_key() for c in kept.constraints}
    assert "61,45" in keys
    assert "61,65" not in keys
    assert "65,61" not in keys
    assert "61,56" not in keys
    assert len(kept.constraints) == 1


def test_prune_ancestor_map_maps_external_ancestor_to_root():
    """祖先映射模式：仅当外部端点为类群祖先、且落在 donor 侧时才合法保留（映射到根）。"""
    tree = _tree()
    kept = prune_constraints(tree, _cset(), "61", ancestor_map=True)
    keys = {c.to_edge_key() for c in kept.constraints}
    assert "61,45" in keys  # 默认保留的
    # 65->61：donor(65) 是类群祖先 -> 映射到根(69)，receptor(61) 在类群内 -> 合法保留
    assert "69,61" in keys
    # 61->65：receptor(65) 映射到根 -> 会断言"类群内早于根"，丢弃
    assert "61,65" not in keys
    # 61->56：56 是旁支（非祖先） -> 丢弃
    assert "61,56" not in keys


def test_prune_invalid_clade_raises_keyerror():
    """--target-clade 指定不存在的标签时必须抛 KeyError。"""
    tree = _tree()
    try:
        prune_constraints(tree, _cset(), "999")
        assert False, "应当抛出 KeyError"
    except KeyError:
        pass


def test_prune_does_not_mutate_tree_labels():
    """剪裁只过滤约束，绝不修改树的拓扑或节点标识（标签对齐保证）。"""
    tree = _tree()
    before = sorted(tree.internal_node_labels())
    prune_constraints(tree, _cset(), "61", ancestor_map=True)
    after = sorted(tree.internal_node_labels())
    assert before == after


def test_prune_ancestor_map_drops_sibling_branch():
    """祖先映射下，donor 为类群外祖先 + receptor 为旁支（sister clade）
    必须被丢弃，禁止引入「根 -> 旁支」虚假偏序。

    旧实现缺 ``and r_in`` 守卫，会把 ``65->56``（65 是 61 的祖先、56 是旁支兄弟）
    错误保留为 ``根(69)->56``；该约束必须被丢弃（见模块 docstring 的「旁支一律丢弃」条目）。
    """
    tree = _tree()
    cset = ConstraintSet(
        [
            Constraint("65", "56", 5.0),  # donor=类群外祖先, receptor=旁支(姐妹支)
            Constraint("65", "61", 8.0),  # donor=类群外祖先, receptor=类群内(合法保留)
        ]
    )
    kept = prune_constraints(tree, cset, "61", ancestor_map=True)
    keys = {c.to_edge_key() for c in kept.constraints}
    # 合法保留：65->61 映射到 根->61
    assert "69,61" in keys
    # 误保留守卫：根->旁支 必须不存在
    assert "69,56" not in keys

    # 普通约束：donor 在类群内、receptor 为旁支（56）也应被丢弃
    cset2 = ConstraintSet([Constraint("61", "56", 3.0)])
    kept2 = prune_constraints(tree, cset2, "61", ancestor_map=True)
    assert "61,56" not in {c.to_edge_key() for c in kept2.constraints}
