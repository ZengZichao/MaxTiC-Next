"""树模块的单元测试。"""

from maxtic_next.tree.tree import Tree


def _load(path):
    with open(path, "r", encoding="utf-8") as fh:
        line = fh.readline()
    t = Tree()
    t.read_newick(line)
    return t


def test_read_minictree_root_and_leaves(minitree_path):
    t = _load(minitree_path)
    root = t.get_root()
    # 根应为后缀为 69 的内部节点
    assert t.get_bootstrap(root) == "69"
    assert t.is_root(root)
    leaves = t.get_leaves_names()
    # 14 个叶子
    assert len(leaves) == 14
    for name in [
        "CYAP8",
        "CYAP0",
        "CYAA5",
        "UCYNA",
        "MICAN",
        "CYAP2",
        "CYAP7",
        "SYNY3",
        "SYNP2",
        "TRIEI",
        "NOSA0",
        "NOSP7",
        "NOSS1",
        "ANAVT",
    ]:
        assert name in leaves


def test_internal_node_count(minitree_path):
    t = _load(minitree_path)
    internal = [n for n in t.get_nodes() if not t.is_leaf(n)]
    # 实际内部节点：69 67 65 62 61 59 56 46 45 43 42 41 40 = 13
    assert len(internal) == 13
    labels = set(t.internal_node_labels())
    assert labels == {"69", "67", "65", "62", "61", "59", "56", "46", "45", "43", "42", "41", "40"}
    assert len(internal) == 13


def test_children_and_is_leaf(minitree_path):
    t = _load(minitree_path)
    root = t.get_root()
    children = t.get_children(root)
    assert len(children) == 2
    for c in children:
        assert not t.is_leaf(c) or t.is_leaf(c)


def test_is_ancestor(minitree_path):
    t = _load(minitree_path)
    root = t.get_root()
    # 找一个内部节点与它的某个后代
    t.get_bootstrap(root)
    # 节点 "40" 是 "46" 的后代（(NOSP7,(NOSS1,ANAVT)40)46）
    n40 = t.label_to_node_id("40")
    n46 = t.label_to_node_id("46")
    assert t.is_ancestor(n46, n40)
    assert not t.is_ancestor(n40, n46)


def test_distance_from_root(minitree_path):
    t = _load(minitree_path)
    root = t.get_root()
    # 任意叶子到根的距离应 > 0
    leaves = t.get_leaves(root)
    for leaf in leaves:
        d = t.distance_from(leaf, root)
        assert d > 0


def test_write_newick_roundtrip_basic(minitree_path):
    t = _load(minitree_path)
    s = t.write_newick(False)
    assert s.endswith(";")
    assert "69" in s
