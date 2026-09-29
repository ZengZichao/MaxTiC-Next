"""ALE 适配器回归测试。

验证：
1. ``parent`` / ``donnor_search`` / ``receptor_search`` 对 T@/D@ 语义的抽取正确；
2. ``convert_from_ale`` 在合成 ``.uml_rec`` 上产出的约束端点均存在于物种树；
3. 文件级缓存（断点续传）二次运行结果一致；
4. （结构）：产出 ``ConstraintSet`` 非空、端点合法、每条含 family/support/distance
   元数据且 weight>0；
5. （实例结构）：``ALEAdapter`` 的 source 校验、``extant_species``、``distance_from`` 对称。
"""

import os

import pytest

from maxtic_next.constraints.adapters.ale import ALEAdapter, convert_from_ale
from maxtic_next.constraints.constraint import ConstraintSet
from maxtic_next.tree.tree import Tree

DATA_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data", "ale_example")
SPECIES = os.path.join(DATA_DIR, "species.tree")
REC = os.path.join(DATA_DIR, "rec.uml_rec")
REC_REAL = os.path.join(DATA_DIR, "rec_real.uml_rec")


def _read_tree(path):
    t = Tree()
    with open(path) as fh:
        t.read_newick(fh.readline())
    return t


def _find_id_by_bootstrap(tree, label):
    for nid in tree.get_nodes():
        if not tree.is_leaf(nid) and tree.get_bootstrap(nid) == label:
            return nid
    return None


def test_parent_returns_species_parent_label():
    species = _read_tree(SPECIES)
    adapter = ALEAdapter(species)
    # '61' 的父节点在物种树中为 '65'
    assert adapter.parent("61") == "65"
    # 根标签的父为 -1
    assert adapter.parent("69") == -1
    # 不存在的标签父为 -1
    assert adapter.parent("NOPE") == -1


def test_donnor_search_walks_up_to_clade_ancestor():
    # 调和树：根 X65，子节点 R=59.T@61->62；donnor_search(R) 应沿父链上溯到 '65'
    rec = Tree()
    rec.read_newick("((CYAP8:1.0,CYAP0:1.0)59.T@61->62:1.0,CYAA5:1.0)X65:1.0;")
    adapter = ALEAdapter(_read_tree(SPECIES))
    r_id = _find_id_by_bootstrap(rec, "59.T@61->62")
    assert r_id is not None
    assert adapter.donnor_search(rec, r_id) == "65"


def test_receptor_search_collects_d_at_children():
    # 根 P.D@x，两子节点 Q.R 与 S.T -> receptor_search 返回 ['R','T']
    rec = Tree()
    rec.read_newick("((A:1.0,B:1.0)Q.R:1.0,(C:1.0,D:1.0)S.T:1.0)P.D@x:1.0;")
    adapter = ALEAdapter(_read_tree(SPECIES))
    root = rec.get_root()
    assert adapter.receptor_search(rec, root) == ["R", "T"]


def test_receptor_search_leaf_returns_empty():
    rec = Tree()
    rec.read_newick("(A:1.0,B:1.0)X:1.0;")
    adapter = ALEAdapter(_read_tree(SPECIES))
    leaf = rec.get_children(rec.get_root())[0]
    assert adapter.receptor_search(rec, leaf) == []


def test_convert_endpoints_exist_in_species_tree():
    species = _read_tree(SPECIES)
    cset = convert_from_ale(species, [REC], min_family_size=0, source="rec")
    assert isinstance(cset, ConstraintSet)
    assert len(cset) >= 1
    valid_labels = set(species.internal_node_labels()) | set(species.get_leaves_names())
    for c in cset.constraints:
        assert c.donor in valid_labels, f"donor {c.donor} 不在物种树中"
        assert c.receptor in valid_labels, f"receptor {c.receptor} 不在物种树中"
        # 合成样例的调和事件应产出 (65,59)
    assert any(c.donor == "65" and c.receptor == "59" for c in cset.constraints)


def test_convert_transfer_source_produces_65_62():
    species = _read_tree(SPECIES)
    cset = convert_from_ale(species, [REC], min_family_size=0, source="trf")
    assert any(c.donor == "65" and c.receptor == "62" for c in cset.constraints)


def test_ale_file_cache_is_deterministic(tmp_path):
    species = _read_tree(SPECIES)
    cache_dir = str(tmp_path / "cache")
    c1 = convert_from_ale(species, [REC], min_family_size=0, cache_dir=cache_dir)
    c2 = convert_from_ale(species, [REC], min_family_size=0, cache_dir=cache_dir)
    keys1 = sorted(c.to_edge_key() for c in c1.constraints)
    keys2 = sorted(c.to_edge_key() for c in c2.constraints)
    assert keys1 == keys2
    # 缓存文件应已生成
    assert os.listdir(cache_dir)


def test_convert_structured_constraints():
    """ALE 适配器产出的 ConstraintSet 结构良好：非空、端点合法、每条含
    family/support/distance 元数据且 weight>0。"""
    species = _read_tree(SPECIES)
    cset = convert_from_ale(species, [REC], min_family_size=0, source="rec")
    assert len(cset) >= 1
    valid_labels = set(species.internal_node_labels()) | set(species.get_leaves_names())
    for c in cset.constraints:
        assert c.donor in valid_labels, f"donor {c.donor} 不在物种树中"
        assert c.receptor in valid_labels, f"receptor {c.receptor} 不在物种树中"
        assert c.weight > 0, "约束权重应为正"
        assert c.metadata.get("family"), "缺少 family 元数据"
        assert c.metadata.get("support") is not None, "缺少 support 元数据"
        assert c.metadata.get("distance") is not None, "缺少 distance 元数据"


def test_ale_adapter_instance_structure():
    """ALEAdapter 实例结构：source 校验、extant_species、distance_from 对称性。"""
    species = _read_tree(SPECIES)
    adapter = ALEAdapter(species, source="rec")
    assert adapter.source == "rec"
    assert adapter.extant_species == set(species.get_leaves_names())
    # source 非法应报错
    with pytest.raises(ValueError):
        ALEAdapter(species, source="bad")
    # distance_from 对称：d(x,y) == d(y,x)
    assert adapter.distance_from("61", "45") == adapter.distance_from("45", "61")


def test_rec_constraints_gated_by_transfer_events(tmp_path):
    """回归（忠实移植）：rec 约束必须仅由带 T@ 标注的节点产生。

    原版 constraints_from_reconciliations.py 的 "2/ calculer la contrainte selon
    les reconciliations" 块位于 ``if e[:2]=="T@"`` 门控体内：纯物种形成/重复节点
    （标注中无 T@）不得产生 rec 约束。若门控被移到节点循环层，无 T@ 家族也会
    系统性产出 (祖先, 自身映射物种) 约束，使默认 source="rec" 的约束量在真实
    数据上膨胀一个数量级。
    """
    species = _read_tree(SPECIES)
    # 家族 A：全部为纯物种形成节点（真实 ALE 前导点格式），无任何 T@ 事件
    fam_no_t = tmp_path / "FAM_NO_T.uml_rec"
    fam_no_t.write_text(
        "1 reconciled G-s:\nFAM_NO_T\n((CYAP8:1.0,CYAP0:1.0).59:1.0,CYAA5:1.0).65:1.0;\n",
        encoding="utf-8",
    )
    cset_a = convert_from_ale(species, [str(fam_no_t)], min_family_size=0, source="rec")
    assert len(cset_a) == 0, (
        "无 T@ 事件的家族不应产出任何 rec 约束，实际："
        f"{[(c.donor, c.receptor, c.weight) for c in cset_a]}"
    )

    # 家族 B：含一个 T@ 事件 -> rec 约束 (65,59) 与 trf 约束 (65,62) 各一条
    fam_t = tmp_path / "FAM_T.uml_rec"
    fam_t.write_text(
        "1 reconciled G-s:\nFAM_T\n((CYAP8:1.0,CYAP0:1.0).59.T@61->62:1.0,CYAA5:1.0).65:1.0;\n",
        encoding="utf-8",
    )
    cset_b = convert_from_ale(species, [str(fam_t)], min_family_size=0, source="rec")
    pairs_b = {(c.donor, c.receptor) for c in cset_b}
    assert ("65", "59") in pairs_b
    # trf 来源：转移事件给出 (parent(donor)=65, receptor=62)
    cset_c = convert_from_ale(species, [str(fam_t)], min_family_size=0, source="trf")
    pairs_c = {(c.donor, c.receptor) for c in cset_c}
    assert ("65", "62") in pairs_c


def test_real_ale_annotation_format_leading_dots():
    """donnor_search/receptor_search 对真实 ALE 前导点标注（如 ".65"）的解析。"""
    species = _read_tree(SPECIES)
    rec = Tree()
    # 真实 ALE_undated 的内部标注以 "." 开头（undated.cpp 恒输出 "."+estr）
    rec.read_newick("((CYAP8:1.0,CYAP0:1.0).59.T@61->62:1.0,CYAA5:1.0).65:1.0;")
    adapter = ALEAdapter(species)
    r_id = _find_id_by_bootstrap(rec, ".59.T@61->62")
    assert r_id is not None
    # T@ 节点的 donor_search 应上溯到其父节点映射的物种 "65"；且对多数字标签
    # "59"（annot[1:]="59..."）解析完整，不得截断为 "9"
    assert adapter.donnor_search(rec, r_id) == "65"


def test_real_format_sample_end_to_end():
    """回归：真实 ALE 输出格式（前导点标注、``<N> reconciled G-s:`` 头行）
    经完整适配管道应与合成样本产出相同的约束（rec 与 trf 两种来源）。"""
    species = _read_tree(SPECIES)
    for source, expected in (("rec", ("65", "59")), ("trf", ("65", "62"))):
        cset = convert_from_ale(species, [REC_REAL], min_family_size=0, source=source)
        pairs = {(c.donor, c.receptor) for c in cset}
        assert expected in pairs, f"source={source} 缺少 {expected}，实际 {pairs}"
        valid = set(_read_tree(SPECIES).internal_node_labels()) | set(
            _read_tree(SPECIES).get_leaves_names()
        )
        for c in cset:
            assert c.donor in valid and c.receptor in valid
