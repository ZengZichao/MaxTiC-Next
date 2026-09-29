# Copyright (C) 2026 MaxTiC-Next rewrite team.
# MaxTiC-Next is free software under the CeCILL 2.1 license (see LICENSE).

"""适配器注册表 + 自动检测 + CLI/API 接线测试（多软件兼容框架）。

覆盖：
1. ``REGISTRY`` 登记 5 个内置工具；
2. ``detect_format`` 对 ale/ranger/eccetera/artra/alerax 的识别与未知返回 None；
3. ``detect_formats`` 混合格式冲突抛错；
4. ``convert`` / ``convert_auto`` 分发一致；
5. ``api.rank(from_tool=..., constraints_out=...)`` 两阶段写出约束（覆盖 5 工具）；
6. CLI 解析 ``--from`` / ``--from-auto`` / ``--artra-transfer-kind``。
"""

import os

import pytest

from maxtic_next.constraints.adapters import registry
from maxtic_next.tree.tree import Tree
from maxtic_next import api, cli

DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")
RANGER = os.path.join(DATA, "ranger_dtl_example", "FAM1.dtl")
RANGER_SP = os.path.join(DATA, "ranger_dtl_example", "species.tree")
ECCE = os.path.join(DATA, "eccetera_example", "FAM1.recphyloxml")
ARTRA = os.path.join(DATA, "artra_example", "FAM1.txt")
ALERAX_RUN = os.path.join(DATA, "alerax_example", "run")
ALE = os.path.join(DATA, "ale_example", "rec.uml_rec")


def _tree(path):
    t = Tree()
    with open(path) as fh:
        t.read_newick(fh.readline())
    return t


class TestRegistry:
    """注册表结构与分发。"""

    def test_five_tools_registered(self):
        assert set(registry.REGISTRY.names()) == {"ale", "ranger", "eccetera", "artra", "alerax"}

    def test_unknown_tool_raises(self):
        with pytest.raises(KeyError):
            registry.REGISTRY.get("nope")

    def test_convert_dispatch_matches_direct(self):
        sp = _tree(RANGER_SP)
        via_reg = registry.convert("ranger", sp, [RANGER], min_family_size=0)
        assert via_reg.informative_edges()["61,62"] == pytest.approx(1.0)


class TestDetection:
    """自动检测。"""

    def test_detect_each_tool(self):
        assert registry.detect_format(RANGER) == "ranger"
        assert registry.detect_format(ECCE) == "eccetera"
        assert registry.detect_format(ARTRA) == "artra"
        assert registry.detect_format(ALERAX_RUN) == "alerax"
        assert registry.detect_format(ALE) == "ale"

    def test_detect_unknown_returns_none(self, tmp_path):
        f = str(tmp_path / "plain.tsv")
        open(f, "w").write("61 62 3.0\n")  # 通用文本约束，非特定工具
        assert registry.detect_format(f) is None

    def test_detect_formats_conflict_raises(self):
        with pytest.raises(ValueError):
            registry.detect_formats([RANGER, ECCE])

    def test_convert_auto(self):
        sp = _tree(RANGER_SP)
        cset = registry.convert_auto(sp, [ECCE], min_family_size=0)
        assert cset.informative_edges() == {"61,62": pytest.approx(1.0)}

    def test_convert_auto_unknown_raises(self, tmp_path):
        f = str(tmp_path / "plain.tsv")
        open(f, "w").write("61 62 3.0\n")
        with pytest.raises(ValueError):
            registry.convert_auto(_tree(RANGER_SP), [f])


class TestApiWiring:
    """``api.rank`` 的 from_tool 两阶段（constraints_out）接线。"""

    @pytest.mark.parametrize(
        "tool,inputs",
        [
            ("ranger", [RANGER]),
            ("eccetera", [ECCE]),
            ("artra", [ARTRA]),
            ("alerax", [ALERAX_RUN]),
            ("auto", [ECCE]),
        ],
    )
    def test_constraints_out_stage(self, tool, inputs, tmp_path):
        out = str(tmp_path / "cons.tsv")
        res = api.rank(
            RANGER_SP,
            inputs,
            from_tool=tool,
            constraints_out=out,
            ale_min_family_size=0,
            print_summary=False,
        )
        assert os.path.exists(out)
        content = open(out).read()
        assert "#family,donor,receptor,weight,distance" in content
        assert res.run_metadata.get("from_tool") in (tool, "eccetera")

    def test_from_ale_backcompat_maps_to_tool(self, tmp_path):
        out = str(tmp_path / "cons.tsv")
        api.rank(
            os.path.join(DATA, "ale_example", "species.tree"),
            [ALE],
            from_ale=True,
            constraints_out=out,
            ale_min_family_size=0,
            print_summary=False,
        )
        assert os.path.exists(out)


class TestCliParsing:
    """CLI 新开关解析。"""

    def test_from_choice(self):
        args = cli.build_parser().parse_args(["sp.tree", "f.txt", "--from", "ranger"])
        assert args.from_tool == "ranger"

    def test_from_auto_flag(self):
        args = cli.build_parser().parse_args(["sp.tree", "f.txt", "--from-auto"])
        assert args.from_auto is True

    def test_artra_transfer_kind(self):
        args = cli.build_parser().parse_args(
            ["sp.tree", "f.txt", "--from", "artra", "--artra-transfer-kind", "replacing"]
        )
        assert args.artra_transfer_kind == "replacing"

    def test_from_ale_still_works(self):
        args = cli.build_parser().parse_args(["sp.tree", "f.uml_rec", "--from-ale"])
        assert args.from_ale is True

    def test_output_style_default_short(self):
        args = cli.build_parser().parse_args(["sp.tree", "f.txt"])
        assert args.output_style == "short"

    def test_output_style_legacy(self):
        args = cli.build_parser().parse_args(["sp.tree", "f.txt", "--output-style", "legacy"])
        assert args.output_style == "legacy"
