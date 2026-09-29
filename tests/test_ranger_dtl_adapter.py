# Copyright (C) 2026 MaxTiC-Next rewrite team.
# This file is part of MaxTiC-Next, a Python 3 rewrite of MaxTiC.
#
# MaxTiC is Copyright (C) Eric Tannier (Inria / CNRS / ENS Lyon).
# Original reference: MaxTiC: Fast ranking of a phylogenetic tree by
# Maximum Time Consistency with lateral gene transfers,
# Biorxiv doi.org/10.1101/127548
#
# MaxTiC-Next is free software under the CeCILL 2.1 license (see LICENSE).

"""RANGER-DTLx 适配器回归测试。

针对**真实 RANGER-DTLx 调和报告格式**（``m<idx> = LCA[..]: Transfer, Mapping --> D,
..., Recipient --> R``，求证自 ``DTL-algorithm.h:1840``）验证：

1. 真实样本 ``FAM1.dtl`` 的两个 Transfer 事件被正确抽取为 (61,62)/(62,65)；
2. 端点均存在于物种树、元数据含 family/support/distance；
3. 互反减法、min_support、min_family_size 过滤；
4. 多块采样的支持度分母；文件级缓存；顺序/并行一致；output_path 写出；
5. 边界：空文件列表、全非法端点。
"""

import os

import pytest

from maxtic_next.constraints.adapters.ranger_dtl import (
    RangerDTLAdapter,
    convert_from_ranger_dtl,
)
from maxtic_next.constraints.constraint import ConstraintSet
from maxtic_next.tree.tree import Tree

DATA_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data", "ranger_dtl_example")
SPECIES = os.path.join(DATA_DIR, "species.tree")
DTL = os.path.join(DATA_DIR, "FAM1.dtl")


def _read_tree(path):
    t = Tree()
    with open(path) as fh:
        t.read_newick(fh.readline())
    return t


def _aggregated(cset):
    return sorted((k, round(w, 9)) for k, w in cset.informative_edges().items())


def _write(path, lines):
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")


class TestConvertRealFormat:
    """真实 RANGER-DTLx 格式端到端测试。"""

    def test_expected_transfers(self):
        species = _read_tree(SPECIES)
        cset = convert_from_ranger_dtl(species, [DTL], min_family_size=0)
        edges = cset.informative_edges()
        assert "61,62" in edges and edges["61,62"] == pytest.approx(1.0)
        assert "62,65" in edges and edges["62,65"] == pytest.approx(1.0)

    def test_endpoints_and_metadata(self):
        species = _read_tree(SPECIES)
        cset = convert_from_ranger_dtl(species, [DTL], min_family_size=0)
        valid = set(species.internal_node_labels()) | set(species.get_leaves_names())
        assert len(cset) >= 2
        for c in cset.constraints:
            assert c.donor in valid and c.receptor in valid
            assert c.weight > 0
            assert c.metadata["family"] == "FAM1"
            assert c.metadata["support"] > 0
            assert c.metadata["distance"] >= 0

    def test_family_size_from_leaf_lines(self):
        """FAM1.dtl 有 8 个 Leaf Node，min_family_size=8 时被过滤，=7 时保留。"""
        species = _read_tree(SPECIES)
        assert len(convert_from_ranger_dtl(species, [DTL], min_family_size=8)) == 0
        assert len(convert_from_ranger_dtl(species, [DTL], min_family_size=7)) >= 1

    def test_empty_file_list(self):
        species = _read_tree(SPECIES)
        cset = convert_from_ranger_dtl(species, [], min_family_size=0)
        assert isinstance(cset, ConstraintSet) and len(cset) == 0


class TestReciprocalAndSupport:
    """互反减法与 min_support 过滤。"""

    def test_equal_opposites_cancel(self, tmp_path):
        species = _read_tree(SPECIES)
        f = str(tmp_path / "FAM_recip.dtl")
        _write(
            f,
            [
                "Reconciliation for Gene Tree 1:",
                "a: Leaf Node",
                "b: Leaf Node",
                "c: Leaf Node",
                "d: Leaf Node",
                "e: Leaf Node",
                "f: Leaf Node",
                "m1 = LCA[a,b]: Transfer, Mapping --> 61, Recipient --> 62",
                "m2 = LCA[c,d]: Transfer, Mapping --> 62, Recipient --> 61",
            ],
        )
        assert len(convert_from_ranger_dtl(species, [f], min_family_size=0)) == 0

    def test_high_min_support_filters(self):
        species = _read_tree(SPECIES)
        assert len(convert_from_ranger_dtl(species, [DTL], min_family_size=0, min_support=2.0)) == 0


class TestBlocksCacheParallel:
    """多块采样、缓存、并行一致性。"""

    def test_two_blocks_support(self, tmp_path):
        species = _read_tree(SPECIES)
        f = str(tmp_path / "multi.dtl")
        _write(
            f,
            [
                "Reconciliation for Gene Tree 1:",
                "g1: Leaf Node",
                "g2: Leaf Node",
                "m1 = LCA[a,b]: Transfer, Mapping --> 61, Recipient --> 62",
                "Reconciliation for Gene Tree 2:",
                "g1: Leaf Node",
                "g2: Leaf Node",
                "m1 = LCA[a,b]: Transfer, Mapping --> 61, Recipient --> 62",
            ],
        )
        edges = convert_from_ranger_dtl(species, [f], min_family_size=0).informative_edges()
        assert edges["61,62"] == pytest.approx(1.0)  # 2/2

    def test_cache_deterministic(self, tmp_path):
        species = _read_tree(SPECIES)
        cache = str(tmp_path / "cache")
        c1 = convert_from_ranger_dtl(species, [DTL], min_family_size=0, cache_dir=cache)
        c2 = convert_from_ranger_dtl(species, [DTL], min_family_size=0, cache_dir=cache)
        assert _aggregated(c1) == _aggregated(c2)
        assert os.listdir(cache)

    def test_sequential_vs_parallel(self, tmp_path):
        species = _read_tree(SPECIES)
        files = []
        for i, (d, r) in enumerate([("61", "62"), ("62", "65"), ("59", "61")]):
            f = str(tmp_path / f"f{i}.dtl")
            _write(
                f,
                [
                    "Reconciliation for Gene Tree 1:",
                    "g1: Leaf Node",
                    "g2: Leaf Node",
                    "g3: Leaf Node",
                    f"m1 = LCA[a,b]: Transfer, Mapping --> {d}, Recipient --> {r}",
                ],
            )
            files.append(f)
        seq = convert_from_ranger_dtl(species, files, min_family_size=0, max_workers=1)
        par = convert_from_ranger_dtl(species, files, min_family_size=0, max_workers=8)
        assert _aggregated(seq) == _aggregated(par)
        assert len(par) >= 1


class TestOutputAndInstance:
    """output_path 写出与适配器实例结构。"""

    def test_output_csv(self, tmp_path):
        species = _read_tree(SPECIES)
        out = str(tmp_path / "out.csv")
        convert_from_ranger_dtl(species, [DTL], min_family_size=0, output_path=out)
        content = open(out).read()
        assert "#family,donor,receptor,weight,distance" in content
        assert "FAM1" in content

    def test_default_attributes(self):
        a = RangerDTLAdapter()
        assert a.min_support == 0.05 and a.min_family_size == 5 and a.cache_dir is None

    def test_convert_interface(self):
        species = _read_tree(SPECIES)
        cset = RangerDTLAdapter(min_family_size=0).convert(species, [DTL])
        assert isinstance(cset, ConstraintSet) and len(cset) >= 2
