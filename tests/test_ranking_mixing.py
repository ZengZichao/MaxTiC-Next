"""``mix`` / ``opt`` 混合启发式测试的单元测试。"""

from maxtic_next.random_ import RandomWrapper
from maxtic_next.ranking.mixing import mix, opt
from maxtic_next.tree.tree import Tree


def test_mix_merges_permutation():
    degre_entrant = {"x": [], "y": []}
    edge = {}
    rng = RandomWrapper(42)
    result = mix(["x"], ["y"], degre_entrant, edge, rng)
    assert sorted(result) == ["x", "y"]
    assert len(result) == 2


def test_mix_larger_is_valid_permutation():
    degre_entrant = {f"n{i}": [] for i in range(6)}
    edge = {}
    rng = RandomWrapper(7)
    o1 = [f"n{i}" for i in range(3)]
    o2 = [f"n{i}" for i in range(3, 6)]
    result = mix(o1, o2, degre_entrant, edge, rng)
    assert sorted(result) == [f"n{i}" for i in range(6)]


def test_opt_small_tree():
    newick = "((A:1,B:1)1:2,C:3)2;"
    t = Tree()
    t.read_newick(newick)
    root = t.get_root()
    degre_entrant = {t.get_bootstrap(n): [] for n in t.get_nodes() if not t.is_leaf(n)}
    edge = {}
    rng = RandomWrapper(1)
    result = opt(t, root, degre_entrant, edge, rng)
    # 内部节点应为 {"1","2"} 的某种顺序（根在末尾，反转后根在前）
    assert set(result) == {"1", "2"}
    assert len(result) == 2
