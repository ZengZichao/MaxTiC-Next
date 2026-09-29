"""边界条件测试：空约束、全冲突约束、单内部节点树、多歧分支（polytomy）。

聚焦于输入在极端情况下不崩溃、产出合法结果，以及 ``dry_run`` 对非法树结构（多歧分支）
的拦截。对应的边界项（ 以测试形式落地）。
"""

import os

from maxtic_next.api import rank
from maxtic_next.dry_run import dry_run_check
from maxtic_next.io.parsing import read_newick_file
from maxtic_next.ranking.ranker import kendall_similarity, order_from_tree
from maxtic_next.tree.tree import Tree

DATA_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")
TREE = os.path.join(DATA_DIR, "minitree.tree")


def _write(path: str, content: str) -> None:
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(content)


def test_empty_constraints_runs(tmp_path):
    """空约束文件：rank 正常完成，best_order 为全部内部节点的排列。"""
    tree = read_newick_file(TREE)
    internal_labels = set(tree.internal_node_labels())
    cons = str(tmp_path / "empty_cons.tsv")
    _write(cons, "")
    prefix = str(tmp_path / "r")
    r = rank(TREE, cons, seed=42, output_prefix=prefix, print_summary=False)
    assert sorted(r.best_order) == sorted(internal_labels)
    assert r.total_weight == 0.0
    assert len(r.informative_lines) == 0


def test_all_conflicting_constraints(tmp_path):
    """全为「后代->祖先」类平凡冲突约束：分类为 trivially conflicting，无信息性边。"""
    cons = str(tmp_path / "all_conf.tsv")
    # 45 是 61 的后代、42 是 59 的后代、43 是 45 的后代 —— 均为 descendant->ancestor
    _write(cons, "45 61 5.0\n42 59 6.0\n43 45 7.0\n")
    prefix = str(tmp_path / "r")
    r = rank(TREE, cons, seed=42, output_prefix=prefix, print_summary=False)
    assert r.trivial_conflict > 0
    # 全部被分类为平凡冲突，不应出现信息性约束行
    assert len(r.informative_lines) == 0
    assert r.best_order  # 仍产出合法排序


def test_single_internal_node_tree(tmp_path):
    """单内部节点树 (A,B)X;：退化情形，rank 不崩溃，best_order == [X]。"""
    tree_str = "(A:1.0,B:1.0)X:1.0;"
    tree_path = str(tmp_path / "single.tree")
    _write(tree_path, tree_str)
    cons = str(tmp_path / "single_cons.tsv")
    _write(cons, "")
    prefix = str(tmp_path / "r")
    tree = read_newick_file(tree_path)
    assert len(tree.internal_node_labels()) == 1
    r = rank(tree_path, cons, seed=42, output_prefix=prefix, print_summary=False)
    assert r.internal_node_count == 1
    assert r.best_order == ["X"]


def test_single_internal_node_order_and_similarity():
    """单内部节点树的 order_from_tree / kendall_similarity 边界不崩溃。"""
    tree = Tree()
    tree.read_newick("(A:1.0,B:1.0)X:1.0;")
    order = order_from_tree(tree)
    assert order == ["X"]
    # 两极端拓扑序相同，相似度退化为 1.0（不除零崩溃）
    assert kendall_similarity(tree, order, order) == 1.0


def test_polytomy_flagged_by_dry_run(tmp_path):
    """多歧分支（内部节点子节点数 != 2）应在 dry_run 预检阶段报错（二叉树强假设）。"""
    tree_str = "((A,B,C)X,Y)Z;"
    tree_path = str(tmp_path / "polytomy.tree")
    _write(tree_path, tree_str)
    tree = read_newick_file(tree_path)
    cons = str(tmp_path / "cons.tsv")
    _write(cons, "A B 1.0\n")
    rep = dry_run_check(tree, [cons])
    assert rep["ok"] is False
    assert any("多歧" in it["message"] or "polytomy" in it["message"] for it in rep["issues"])
