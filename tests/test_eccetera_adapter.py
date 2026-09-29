# Copyright (C) 2026 MaxTiC-Next rewrite team.
# MaxTiC-Next is free software under the CeCILL 2.1 license (see LICENSE).

"""ecceTERA 适配器回归测试。

针对 ecceTERA **recPhyloXML** 输出（``branchingOut speciesLocation`` = 供体、子
``transferBack destinationSpecies`` = 受体，求证自
``DTLGraph.cpp:getRecPhyloXMLReconciliation``）验证：

1. 真实样本 ``FAM1.recphyloxml`` 的转移被抽取为 (61,62)；
2. 端点合法、元数据完整、距离正确；
3. ``parse_recphyloxml`` 对供体/受体配对与非法端点过滤；
4. min_family_size / min_support 过滤；缓存；DOCTYPE 安全加固；边界。
"""

import os

import pytest

from maxtic_next.constraints.adapters.eccetera import (
    EcceTERAAdapter,
    convert_from_eccetera,
    parse_recphyloxml,
)
from maxtic_next.constraints.constraint import ConstraintSet
from maxtic_next.tree.tree import Tree

DATA_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data", "eccetera_example")
SPECIES = os.path.join(DATA_DIR, "species.tree")
XML = os.path.join(DATA_DIR, "FAM1.recphyloxml")
VALID = {"59", "61", "62", "65", "69", "CYAP8", "CYAP0", "CYAA5", "CYAP2", "NOSP7", "ANAVT"}


def _read_tree(path):
    t = Tree()
    with open(path) as fh:
        t.read_newick(fh.readline())
    return t


class TestParseRecPhyloXML:
    """``parse_recphyloxml`` 单元测试。"""

    def test_transfer_pairing(self):
        text = open(XML).read()
        transfers, number, fam = parse_recphyloxml(text, VALID)
        assert transfers == {"61,62": 1}
        assert number == 1
        assert fam == 6  # 6 个 leaf geneName

    def test_invalid_endpoints_filtered(self):
        text = (
            '<recPhylo><recGeneTree><phylogeny rooted="true">'
            "<clade><name>0</name>"
            '<eventsRec><branchingOut speciesLocation="ZZZ"></branchingOut></eventsRec>'
            "<clade><name>1</name>"
            '<eventsRec><transferBack destinationSpecies="62"></transferBack>'
            '<leaf speciesLocation="62" geneName="g1"></leaf></eventsRec></clade>'
            "</clade></phylogeny></recGeneTree></recPhylo>"
        )
        transfers, _, _ = parse_recphyloxml(text, VALID)
        assert transfers == {}  # donor ZZZ 非法，丢弃

    def test_doctype_rejected(self):
        """含 DOCTYPE 的文档被安全拒绝（防实体展开）。"""
        text = (
            '<!DOCTYPE recPhylo [<!ENTITY x "y">]><recPhylo><recGeneTree></recGeneTree></recPhylo>'
        )
        transfers, number, fam = parse_recphyloxml(text, VALID)
        assert transfers == {} and number == 0

    def test_malformed_xml_returns_empty(self):
        transfers, number, fam = parse_recphyloxml("<recPhylo><clade>", VALID)
        assert transfers == {} and number == 0 and fam == 0


class TestConvertFromEcceTERA:
    """``convert_from_eccetera`` 端到端测试。"""

    def test_expected_transfer(self):
        species = _read_tree(SPECIES)
        cset = convert_from_eccetera(species, [XML], min_family_size=0)
        edges = cset.informative_edges()
        assert edges == {"61,62": pytest.approx(1.0)}

    def test_metadata_complete(self):
        species = _read_tree(SPECIES)
        cset = convert_from_eccetera(species, [XML], min_family_size=0)
        for c in cset.constraints:
            assert c.metadata["family"] == "FAM1"
            assert c.metadata["support"] > 0
            assert c.metadata["distance"] >= 0

    def test_family_size_filter(self):
        species = _read_tree(SPECIES)
        # 6 leaves：min_family_size=6 过滤，=5 保留
        assert len(convert_from_eccetera(species, [XML], min_family_size=6)) == 0
        assert len(convert_from_eccetera(species, [XML], min_family_size=5)) == 1

    def test_high_min_support_filters(self):
        species = _read_tree(SPECIES)
        assert len(convert_from_eccetera(species, [XML], min_family_size=0, min_support=2.0)) == 0

    def test_empty_file_list(self):
        species = _read_tree(SPECIES)
        cset = convert_from_eccetera(species, [], min_family_size=0)
        assert isinstance(cset, ConstraintSet) and len(cset) == 0

    def test_cache_deterministic(self, tmp_path):
        species = _read_tree(SPECIES)
        cache = str(tmp_path / "cache")
        c1 = convert_from_eccetera(species, [XML], min_family_size=0, cache_dir=cache)
        c2 = convert_from_eccetera(species, [XML], min_family_size=0, cache_dir=cache)
        assert sorted(c.to_edge_key() for c in c1) == sorted(c.to_edge_key() for c in c2)
        assert os.listdir(cache)

    def test_default_attributes(self):
        a = EcceTERAAdapter()
        assert a.min_support == 0.05 and a.min_family_size == 5
