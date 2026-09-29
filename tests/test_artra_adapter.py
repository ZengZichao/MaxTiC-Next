# Copyright (C) 2026 MaxTiC-Next rewrite team.
# MaxTiC-Next is free software under the CeCILL 2.1 license (see LICENSE).

"""ARTra 适配器回归测试。

针对 ARTra ``output.txt``（Replacing / Additive Transfer，与 RANGER-DTLx 同形）验证：

1. 真实样本 ``FAM1.txt`` 的两个转移被抽取为 (61,62 Replacing)/(62,65 Additive)；
2. ``transfer_kind`` 过滤：all / replacing / additive；
3. 端点合法、元数据完整；min_support / min_family_size；空输入。
"""

import os

import pytest

from maxtic_next.constraints.adapters.artra import (
    ARTraAdapter,
    convert_from_artra,
)
from maxtic_next.constraints.constraint import ConstraintSet
from maxtic_next.tree.tree import Tree

DATA_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data", "artra_example")
SPECIES = os.path.join(DATA_DIR, "species.tree")
ARTRA = os.path.join(DATA_DIR, "FAM1.txt")


def _read_tree(path):
    t = Tree()
    with open(path) as fh:
        t.read_newick(fh.readline())
    return t


class TestConvertFromARTra:
    """``convert_from_artra`` 端到端测试。"""

    def test_all_transfers(self):
        species = _read_tree(SPECIES)
        edges = convert_from_artra(species, [ARTRA], min_family_size=0).informative_edges()
        assert edges["61,62"] == pytest.approx(1.0)  # Replacing
        assert edges["62,65"] == pytest.approx(1.0)  # Additive

    def test_replacing_only(self):
        species = _read_tree(SPECIES)
        edges = convert_from_artra(
            species, [ARTRA], min_family_size=0, transfer_kind="replacing"
        ).informative_edges()
        assert "61,62" in edges
        assert "62,65" not in edges

    def test_additive_only(self):
        species = _read_tree(SPECIES)
        edges = convert_from_artra(
            species, [ARTRA], min_family_size=0, transfer_kind="additive"
        ).informative_edges()
        assert "62,65" in edges
        assert "61,62" not in edges

    def test_metadata_and_endpoints(self):
        species = _read_tree(SPECIES)
        cset = convert_from_artra(species, [ARTRA], min_family_size=0)
        valid = set(species.internal_node_labels()) | set(species.get_leaves_names())
        for c in cset.constraints:
            assert c.donor in valid and c.receptor in valid
            assert c.metadata["family"] == "FAM1"
            assert c.metadata["distance"] >= 0

    def test_family_size_filter(self):
        species = _read_tree(SPECIES)
        # 8 leaves：min_family_size=8 过滤，=7 保留
        assert len(convert_from_artra(species, [ARTRA], min_family_size=8)) == 0
        assert len(convert_from_artra(species, [ARTRA], min_family_size=7)) >= 1

    def test_empty_file_list(self):
        species = _read_tree(SPECIES)
        cset = convert_from_artra(species, [], min_family_size=0)
        assert isinstance(cset, ConstraintSet) and len(cset) == 0

    def test_invalid_transfer_kind_raises(self):
        with pytest.raises(ValueError):
            ARTraAdapter(transfer_kind="bogus")

    def test_default_attributes(self):
        a = ARTraAdapter()
        assert a.min_support == 0.05 and a.transfer_kind == "all"
