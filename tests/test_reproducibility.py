"""排序可复现性测试：``--seed`` 驱动确定性路径。

断言：
* 同 seed 复现完全相同的 ``best_order``（含局部搜索）；
* 不同 seed 不崩溃，均产出合法排序；
* 局部搜索在固定 seed 下逐字节可复现。

说明：与原版 Python 2.7 黄金值的逐字节比对需在容器化 Py2.7 环境进行（本机无
Py2.7），此处仅验证本移植版本在固定 seed 下的确定性（等价自洽，见
test_equivalence.py 声明）。
"""

import os

from maxtic_next.api import rank
from maxtic_next.tree.tree import Tree

DATA_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")
TREE = os.path.join(DATA_DIR, "minitree.tree")
CONS = os.path.join(DATA_DIR, "Cyano_CUTConstraints.tsv")


def test_same_seed_reproduces_identical_order(tmp_path):
    """seed=42 重复运行 -> best_order 与各项 value 完全一致。"""
    p1 = str(tmp_path / "r1")
    p2 = str(tmp_path / "r2")
    r1 = rank(TREE, CONS, seed=42, output_prefix=p1, print_summary=False)
    r2 = rank(TREE, CONS, seed=42, output_prefix=p2, print_summary=False)
    assert r1.best_order == r2.best_order
    assert r1.values == r2.values
    assert r1.ranked_newick == r2.ranked_newick
    assert r1.best_source == r2.best_source


def test_different_seed_no_crash_valid_permutation(tmp_path):
    """不同 seed 均不崩溃，且产出内部节点全集的合法排列。"""
    tree = Tree()
    with open(TREE) as fh:
        tree.read_newick(fh.readline())
    internal_labels = set(tree.internal_node_labels())
    prefix = str(tmp_path / "rd")
    for i, seed in enumerate((1, 2, 7, 123)):
        r = rank(TREE, CONS, seed=seed, output_prefix=f"{prefix}_{i}", print_summary=False)
        assert sorted(r.best_order) == sorted(internal_labels), (
            f"seed={seed} 未产出内部节点全集排列"
        )
        assert r.best_source in ("greedy heuristic", "mixing heuristic")


def test_local_search_deterministic_under_seed(tmp_path):
    """局部搜索开启时，固定 seed 复现相同 best_order（确定性）。"""
    p1 = str(tmp_path / "r1")
    p2 = str(tmp_path / "r2")
    r1 = rank(TREE, CONS, seed=42, local_search=0.05, output_prefix=p1, print_summary=False)
    r2 = rank(TREE, CONS, seed=42, local_search=0.05, output_prefix=p2, print_summary=False)
    assert r1.best_order == r2.best_order
    assert r1.sensitivity_summary is not None
    assert r2.sensitivity_summary is not None
