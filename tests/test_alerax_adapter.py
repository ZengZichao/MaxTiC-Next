# Copyright (C) 2026 MaxTiC-Next rewrite team.
# MaxTiC-Next is free software under the CeCILL 2.1 license (see LICENSE).

"""AleRax 适配器回归测试。

针对 AleRax ``reconciliations/summaries/<fam>_transfers.txt``（``donor recipient freq``，
求证自 ``extract_families_transfer.py``）验证：

1. 目录自动展开为 ``*_transfers.txt`` 并跨家族按边键聚合频率；
2. 频率即权重（不再除以块数）；min_support 频率阈值；
3. 直接传文件列表、传 summaries 目录、传根目录三种输入等价；
4. 端点合法、元数据完整；空输入。
"""

import os

import pytest

from maxtic_next.constraints.adapters.alerax import (
    AleRaxAdapter,
    convert_from_alerax,
    expand_alerax_inputs,
)
from maxtic_next.constraints.constraint import ConstraintSet
from maxtic_next.tree.tree import Tree

DATA_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data", "alerax_example")
SPECIES = os.path.join(DATA_DIR, "species.tree")
RUN = os.path.join(DATA_DIR, "run")
SUMMARIES = os.path.join(RUN, "reconciliations", "summaries")


def _read_tree(path):
    t = Tree()
    with open(path) as fh:
        t.read_newick(fh.readline())
    return t


class TestExpandInputs:
    """``expand_alerax_inputs`` 目录展开测试。"""

    def test_expand_root_dir(self):
        files = expand_alerax_inputs([RUN])
        assert len(files) == 2
        assert all(f.endswith("_transfers.txt") for f in files)

    def test_expand_summaries_dir(self):
        files = expand_alerax_inputs([SUMMARIES])
        assert len(files) == 2

    def test_expand_explicit_files(self):
        f = os.path.join(SUMMARIES, "FAM1_transfers.txt")
        assert expand_alerax_inputs([f]) == [f]


class TestConvertFromAleRax:
    """``convert_from_alerax`` 端到端测试。"""

    def test_cross_family_aggregation(self):
        species = _read_tree(SPECIES)
        edges = convert_from_alerax(species, [RUN]).informative_edges()
        # FAM1: 61->62=0.9, 62->65=0.7 ; FAM2: 59->61=0.6, 61->62=0.4
        assert edges["61,62"] == pytest.approx(1.3)  # 0.9 + 0.4
        assert edges["62,65"] == pytest.approx(0.7)
        assert edges["59,61"] == pytest.approx(0.6)

    def test_frequency_is_weight(self):
        species = _read_tree(SPECIES)
        cset = convert_from_alerax(species, [SUMMARIES])
        for c in cset.constraints:
            assert c.metadata["support"] == c.weight  # 频率即支持度

    def test_min_support_threshold(self):
        species = _read_tree(SPECIES)
        # min_support=0.5 过滤掉 59->61 (0.6>0.5 保留) 与 61->62(FAM2,0.4<0.5)
        # 注意聚合前按家族频率过滤：FAM2 的 0.4 被丢弃，61,62 仅剩 FAM1 的 0.9
        edges = convert_from_alerax(species, [RUN], min_support=0.5).informative_edges()
        assert edges["61,62"] == pytest.approx(0.9)
        assert edges["62,65"] == pytest.approx(0.7)
        assert edges["59,61"] == pytest.approx(0.6)

    def test_input_forms_equivalent(self):
        species = _read_tree(SPECIES)
        f1 = os.path.join(SUMMARIES, "FAM1_transfers.txt")
        f2 = os.path.join(SUMMARIES, "FAM2_transfers.txt")
        by_files = convert_from_alerax(species, [f1, f2]).informative_edges()
        by_dir = convert_from_alerax(species, [RUN]).informative_edges()
        assert by_files == by_dir

    def test_metadata(self):
        species = _read_tree(SPECIES)
        cset = convert_from_alerax(species, [RUN])
        fams = {c.metadata["family"] for c in cset.constraints}
        assert fams == {"FAM1", "FAM2"}

    def test_empty_inputs(self):
        species = _read_tree(SPECIES)
        cset = convert_from_alerax(species, [])
        assert isinstance(cset, ConstraintSet) and len(cset) == 0

    def test_mean_transfers_naming_supported(self, tmp_path):
        """回归：现行 AleRax 写 ``<fam>_meanTransfers.txt``
        （AleOptimizer.cpp summaries 输出），必须同样被展开与解析。"""
        src = os.path.join(SUMMARIES, "FAM1_transfers.txt")
        new_dir = tmp_path / "run" / "reconciliations" / "summaries"
        new_dir.mkdir(parents=True)
        (new_dir / "FAM9_meanTransfers.txt").write_text(
            open(src, encoding="utf-8").read(), encoding="utf-8"
        )
        run_root = str(tmp_path / "run")
        files = expand_alerax_inputs([run_root])
        assert len(files) == 1 and files[0].endswith("_meanTransfers.txt")
        species = _read_tree(SPECIES)
        by_old = convert_from_alerax(
            species, [os.path.join(SUMMARIES, "FAM1_transfers.txt")]
        ).informative_edges()
        by_new = convert_from_alerax(species, [run_root]).informative_edges()
        assert by_new == by_old

    def test_recognized_dir_without_transfer_files_raises(self, tmp_path):
        """识别为 AleRax 目录但展开为空时必须报错，绝不静默返回空约束集。"""
        species = _read_tree(SPECIES)
        empty_run = tmp_path / "run" / "reconciliations" / "summaries"
        empty_run.mkdir(parents=True)
        with pytest.raises(ValueError, match="AleRax"):
            convert_from_alerax(species, [str(tmp_path / "run")])

    def test_default_attributes(self):
        a = AleRaxAdapter()
        assert a.min_support == 0.05
