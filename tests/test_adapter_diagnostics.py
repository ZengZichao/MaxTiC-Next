# Copyright (C) 2026 MaxTiC-Next rewrite team.
# MaxTiC-Next is free software under the CeCILL 2.1 license (see LICENSE).

"""上游格式适配器的口径与可观测性回归测试。

每条用例断言的都是用户**看得见的**事实：

* ALE 默认口径 ``trf``（ALE 官方 MaxTiC 集成规定），且口径 + 边数 + 总权重
  在 stderr 声明；
* ARTra / RANGER-DTLx 用自己声明的 ``Total number of optimal
  solutions: N`` 作分母；分母不可确立时权重退化为整数计数并**说出来**；
* ecceTERA 的数字物种节点 ID 可反解（新增夹具
  ``data/eccetera_numeric_example``）；完全命中不了时抛 ``LabelMismatchError``
  而不是静默 ``{}``；家族规模探测失败时跳过过滤；
* ALE 文件级缓存键含物种树**拓扑**指纹 + 版本号（新增夹具
  ``data/ale_topology_example`` 的两棵同标签集、不同拓扑的树）；
* 无距离列时 ``filter_by_distance`` 如实报告"该选项被忽略"；
  未知标签的距离是 ``None`` 而不是 0.0；
* XML 解析失败抛出带文件/行列/原始异常的错误（旧的静默路径至少响亮告警）；
* 所有丢弃被计数；混用两种 donor 层级约定时 ``--from-auto``
  会警告。
"""

import os

import pytest

from maxtic_next.constraints.adapters import registry
from maxtic_next.constraints.adapters._diagnostics import (
    CONVENTION_DONOR_ITSELF,
    CONVENTION_PARENT_OF_DONOR,
    LabelMismatchError,
    UpstreamParseError,
    WEIGHT_DECLARED_SUPPORT,
    WEIGHT_INTEGER_COUNT,
)
from maxtic_next.constraints.adapters._recon_report import (
    parse_reconciliation_report_detailed,
)
from maxtic_next.constraints.adapters._species import (
    build_valid_labels,
    distance_from,
    numeric_species_id_map,
    species_cache_identity,
    species_labels_digest,
    species_topology_digest,
)
from maxtic_next.constraints.adapters.ale import (
    ALEAdapter,
    DEFAULT_SOURCE,
    convert_from_ale,
)
from maxtic_next.constraints.adapters.eccetera import (
    convert_from_eccetera,
    parse_recphyloxml,
)
from maxtic_next.constraints.adapters.artra import convert_from_artra
from maxtic_next.constraints.adapters.ranger_dtl import convert_from_ranger_dtl
from maxtic_next.constraints.constraint import Constraint, ConstraintSet
from maxtic_next.tree.tree import Tree

DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")
ALE_TOPO = os.path.join(DATA, "ale_topology_example")
ECCE_NUM = os.path.join(DATA, "eccetera_numeric_example")
ARTRA_DECL = os.path.join(DATA, "artra_declared_example")
RANGER = os.path.join(DATA, "ranger_dtl_example")
ALE_DIR = os.path.join(DATA, "ale_example")


def _tree(path):
    t = Tree()
    with open(path, encoding="utf-8") as fh:
        t.read_newick(fh.readline())
    return t


def _write(path, text):
    with open(str(path), "w", encoding="utf-8") as fh:
        fh.write(text)
    return str(path)


# ======================================================================
#  —— ALE 默认 source = trf（+ stderr 口径声明）
# ======================================================================
class TestM20AleDefaultSource:
    def test_default_source_is_trf(self):
        sp = _tree(os.path.join(ALE_DIR, "species.tree"))
        assert DEFAULT_SOURCE == "trf"
        assert ALEAdapter(sp).source == "trf"

    def test_default_equals_explicit_trf_and_rec_still_selectable(self):
        sp = _tree(os.path.join(ALE_DIR, "species.tree"))
        rec_file = os.path.join(ALE_DIR, "rec.uml_rec")
        default = convert_from_ale(sp, [rec_file], min_family_size=0, quiet=True)
        trf = convert_from_ale(sp, [rec_file], min_family_size=0, source="trf", quiet=True)
        rec = convert_from_ale(sp, [rec_file], min_family_size=0, source="rec", quiet=True)
        assert default.informative_edges() == trf.informative_edges()
        # trf 取 parent(donor)->receptor，rec 取 donnor_search->receptor_search：
        # 两条路径的端点不同
        assert set(trf.informative_edges()) != set(rec.informative_edges())
        assert ("65", "62") in {(c.donor, c.receptor) for c in trf}
        assert ("65", "59") in {(c.donor, c.receptor) for c in rec}

    def test_stderr_announces_convention_edges_and_weight(self, capsys):
        sp = _tree(os.path.join(ALE_DIR, "species.tree"))
        convert_from_ale(sp, [os.path.join(ALE_DIR, "rec.uml_rec")], min_family_size=0)
        err = capsys.readouterr().err
        assert "source=trf" in err
        assert "供体端点约定=parent_of_donor" in err
        assert "总权重=" in err
        assert "[ale] 约束口径" in err

    def test_rec_notice_names_the_convention(self, capsys):
        sp = _tree(os.path.join(ALE_DIR, "species.tree"))
        convert_from_ale(
            sp, [os.path.join(ALE_DIR, "rec.uml_rec")], min_family_size=0, source="rec"
        )
        assert "source=rec" in capsys.readouterr().err


# ======================================================================
#  /  —— 声明样本数是分母；不可确立时是整数计数
# ======================================================================
class TestM21SampleDenominator:
    def test_declared_optimal_solution_count_is_denominator(self):
        """真实形状：声明 24 个最优解、只打印 1 个 => 分母 24、权重 1/24。"""
        sp = _tree(os.path.join(ARTRA_DECL, "species.tree"))
        cset = convert_from_artra(
            sp,
            [os.path.join(ARTRA_DECL, "FAM1.txt")],
            min_family_size=0,
            min_support=0.0,
            quiet=True,
        )
        diag = cset.diagnostics
        assert diag["blocks_seen"] == 1
        assert diag["declared_sample_counts"] == [24]
        assert diag["sample_denominator"] == 24
        assert diag["weight_semantics"] == WEIGHT_DECLARED_SUPPORT
        for w in cset.informative_edges().values():
            assert w == pytest.approx(1.0 / 24.0)
        assert diag["min_support_is_fractional"] is True

    def test_support_threshold_is_no_longer_a_silent_noop(self):
        """旧实现：块数=1 => 权重恒 1.0 => min_support 0.95 全保留（空操作）。
        现在 1/24 < 0.05 默认阈值即被丢弃，且丢弃数可见。"""
        sp = _tree(os.path.join(ARTRA_DECL, "species.tree"))
        kept = convert_from_artra(
            sp,
            [os.path.join(ARTRA_DECL, "FAM1.txt")],
            min_family_size=0,
            min_support=0.95,
            quiet=True,
        )
        assert len(kept) == 0
        assert kept.diagnostics["dropped_by_support"] == 2
        assert kept.diagnostics["transfers_seen"] == 2

    def test_integer_counts_when_denominator_unestablished(self, tmp_path, capsys):
        """多块且叶子集合不同（= 多棵不同基因树）：不得用块数当分母，
        改为整数计数 + 警告（阈值语义同时降级并说明）。"""
        sp = _tree(os.path.join(RANGER, "species.tree"))
        f = _write(
            tmp_path / "multi_trees.dtl",
            "\n".join(
                [
                    "Reconciliation for Gene Tree 1:",
                    "g1: Leaf Node",
                    "g2: Leaf Node",
                    "m1 = LCA[a,b]: Transfer, Mapping --> 61, Recipient --> 62",
                    "Reconciliation for Gene Tree 2:",
                    "h9: Leaf Node",
                    "h10: Leaf Node",
                    "m1 = LCA[x,y]: Transfer, Mapping --> 61, Recipient --> 62",
                ]
            )
            + "\n",
        )
        cset = convert_from_ranger_dtl(sp, [f], min_family_size=0)
        diag = cset.diagnostics
        assert diag["weight_semantics"] == WEIGHT_INTEGER_COUNT
        assert diag["sample_denominator"] is None
        assert diag["blocks_seen"] == 2
        # 权重是计数 2，而不是被块数归一化后的 1.0
        assert cset.informative_edges()["61,62"] == pytest.approx(2.0)
        assert diag["min_support_is_fractional"] is False
        err = capsys.readouterr().err
        assert "整数计数" in err
        assert "叶子集合" in err

    def test_replicate_blocks_still_supported_but_flagged(self, tmp_path):
        """叶子集合一致的多块（重复调和）：可用块数作分母，但必须记警告。"""
        sp = _tree(os.path.join(RANGER, "species.tree"))
        f = _write(
            tmp_path / "replicates.dtl",
            "\n".join(
                [
                    "Reconciliation for Gene Tree 1:",
                    "g1: Leaf Node",
                    "g2: Leaf Node",
                    "m1 = LCA[a,b]: Transfer, Mapping --> 61, Recipient --> 62",
                    "Reconciliation for Gene Tree 2:",
                    "g1: Leaf Node",
                    "g2: Leaf Node",
                    "m1 = LCA[a,b]: Transfer, Mapping --> 61, Recipient --> 62",
                ]
            )
            + "\n",
        )
        cset = convert_from_ranger_dtl(sp, [f], min_family_size=0, quiet=True)
        assert cset.diagnostics["weight_semantics"] == "support_over_replicate_blocks"
        assert cset.informative_edges()["61,62"] == pytest.approx(1.0)  # 2/2
        assert any("分母" in w for w in cset.diagnostics["warnings"])

    def test_parser_exposes_declared_counts(self):
        lines = [
            "Reconciliation for Gene Tree 1:",
            "g1: Leaf Node",
            "m1 = LCA[a,b]: Transfer, Mapping --> 61, Recipient --> 62",
            "Total number of optimal solutions: 24",
        ]
        transfers, blocks, fam, stats = parse_reconciliation_report_detailed(lines, {"61", "62"})
        assert (transfers, blocks, fam) == ({"61,62": 1}, 1, 1)
        assert stats["declared_sample_counts"] == [24]
        assert stats["weight_semantics"] == WEIGHT_DECLARED_SUPPORT


# ======================================================================
#  /  —— ecceTERA 数字节点 ID、命中率守卫、家族规模探测失败
# ======================================================================
class TestM22EcceTERALabels:
    def test_numeric_internal_species_ids_resolve(self):
        """真实 ecceTERA 输出：内部物种节点写 ``getId()`` —— 必须能反解。"""
        sp = _tree(os.path.join(ECCE_NUM, "species.tree"))
        cset = convert_from_eccetera(
            sp, [os.path.join(ECCE_NUM, "FAM1.recphyloxml")], min_family_size=0, quiet=True
        )
        edges = cset.informative_edges()
        assert edges == {"H61,H62": pytest.approx(1.0), "H65,H59": pytest.approx(1.0)}
        assert cset.diagnostics["numeric_id_resolutions"] == 4
        assert cset.diagnostics["dropped_by_label_miss"] == 0
        assert cset.diagnostics["endpoint_hit_rate"] == pytest.approx(1.0)

    def test_numeric_id_numbering_matches_upstream_rule(self):
        """编号规则复刻 ``MySpeciesTree.cpp:112-136``：排序叶子 0..L-1，
        其父节点自底向上逐层广度优先依次编号。"""
        sp = _tree(os.path.join(ECCE_NUM, "species.tree"))
        id_map = numeric_species_id_map(sp)
        assert id_map["0"] == "ANAVT"  # 排序后首个叶名
        assert id_map["5"] == "NOSP7"  # 6 片叶子占 0..5
        assert id_map["6"] == "H62"  # 第一个被访问到的内部节点
        assert id_map["10"] == "H69"  # 根最后编号

    def test_zero_hit_rate_raises_naming_both_sides(self, tmp_path):
        """端点全不命中：报错并给出两侧标签样本，绝不返回 ``{}``。"""
        sp = _tree(os.path.join(ECCE_NUM, "species.tree"))
        xml = (
            '<recPhylo><recGeneTree><phylogeny rooted="true">'
            "<clade><name>0</name>"
            '<eventsRec><branchingOut speciesLocation="ZZZ_DONOR">'
            "</branchingOut></eventsRec>"
            "<clade><name>1</name>"
            '<eventsRec><transferBack destinationSpecies="YYY_REC">'
            "</transferBack>"
            '<leaf speciesLocation="YYY_REC" geneName="g1"></leaf></eventsRec>'
            "</clade></clade></phylogeny></recGeneTree></recPhylo>"
        )
        f = _write(tmp_path / "FAM_badlabels.recphyloxml", xml)
        with pytest.raises(LabelMismatchError) as exc:
            convert_from_eccetera(sp, [f], min_family_size=0, quiet=True)
        msg = str(exc.value)
        assert "ZZZ_DONOR" in msg and "YYY_REC" in msg  # 上游侧样本
        assert "H61" in msg or "CYAP8" in msg  # 物种树侧样本
        assert exc.value.diagnostics["dropped_by_label_miss"] == 1

    def test_partial_hit_rate_warns_below_threshold(self, tmp_path, capsys):
        sp = _tree(os.path.join(ECCE_NUM, "species.tree"))
        xml = (
            '<recPhylo><recGeneTree><phylogeny rooted="true">'
            "<clade><name>0</name>"
            '<eventsRec><branchingOut speciesLocation="7"></branchingOut>'
            "</eventsRec>"
            "<clade><name>1</name>"
            '<eventsRec><transferBack destinationSpecies="6"></transferBack>'
            '<leaf speciesLocation="6" geneName="g1"></leaf></eventsRec>'
            "</clade></clade>"
            "<clade><name>2</name>"
            '<eventsRec><branchingOut speciesLocation="9"></branchingOut>'
            "</eventsRec>"
            "<clade><name>3</name>"
            '<eventsRec><transferBack destinationSpecies="SPU_NOPE">'
            "</transferBack>"
            '<leaf speciesLocation="0" geneName="g2"></leaf></eventsRec>'
            "</clade></clade>"
            "</phylogeny></recGeneTree></recPhylo>"
        )
        f = _write(tmp_path / "FAM_partial.recphyloxml", xml)
        cset = convert_from_eccetera(sp, [f], min_family_size=0, min_endpoint_hit_rate=0.9)
        assert cset.diagnostics["endpoint_hit_rate"] == pytest.approx(0.5)
        assert "50.0%" in capsys.readouterr().err
        # 关闭守卫（threshold<=0）时不报错也不告警
        quiet = convert_from_eccetera(
            sp, [f], min_family_size=0, min_endpoint_hit_rate=0.0, quiet=True
        )
        assert quiet.diagnostics["dropped_by_label_miss"] == 1

    def test_family_size_probe_failure_skips_filter(self, tmp_path):
        """无 ``<leaf>`` 可探测家族规模时：跳过过滤而不是清空。"""
        sp = _tree(os.path.join(ECCE_NUM, "species.tree"))
        xml = (
            '<recPhylo><recGeneTree><phylogeny rooted="true">'
            "<clade><name>0</name>"
            '<eventsRec><branchingOut speciesLocation="7"></branchingOut>'
            "</eventsRec>"
            "<clade><name>1</name>"
            '<eventsRec><transferBack destinationSpecies="6"></transferBack>'
            '<speciation speciesLocation="6"></speciation></eventsRec>'
            "</clade></clade></phylogeny></recGeneTree></recPhylo>"
        )
        f = _write(tmp_path / "FAM_noleaves.recphyloxml", xml)
        # 连 --ale-min-family-size 0 也不得把全部转移清掉
        cset = convert_from_eccetera(sp, [f], min_family_size=0, quiet=True)
        assert cset.informative_edges() == {"H61,H62": pytest.approx(1.0)}
        assert cset.diagnostics["family_size_probe_ok"] is False
        assert any("跳过" in w for w in cset.diagnostics["warnings"])

    def test_report_without_leaf_lines_keeps_transfers(self, tmp_path):
        """RANGER/ARTra 报告无 ``: Leaf Node`` 行时的同一口径。"""
        sp = _tree(os.path.join(RANGER, "species.tree"))
        f = _write(
            tmp_path / "noleaf.dtl",
            "Reconciliation for Gene Tree 1:\n"
            "m1 = LCA[a,b]: Transfer, Mapping --> 61, Recipient --> 62\n",
        )
        assert len(convert_from_ranger_dtl(sp, [f], min_family_size=0, quiet=True)) == 1
        # 旧行为：family_size 探测失败 => 0 <= 5 => 全部丢弃
        assert len(convert_from_ranger_dtl(sp, [f], quiet=True)) == 1


# ======================================================================
#  —— 缓存键必须含物种树拓扑指纹
# ======================================================================
class TestM23TopologyAwareCache:
    def test_cache_identity_distinguishes_same_labels_different_topology(self):
        a = _tree(os.path.join(ALE_TOPO, "species_a.tree"))
        b = _tree(os.path.join(ALE_TOPO, "species_b.tree"))
        # 标签集合完全相同（旧实现只看这个 => 缓存跨树复用）
        assert species_labels_digest(build_valid_labels(a)) == species_labels_digest(
            build_valid_labels(b)
        )
        # 拓扑指纹不同
        assert species_topology_digest(a) != species_topology_digest(b)
        assert species_cache_identity(a) != species_cache_identity(b)
        assert species_cache_identity(a).startswith("v2|")

    def test_ale_cache_dir_is_not_reused_across_topologies(self, tmp_path):
        sp_a = _tree(os.path.join(ALE_TOPO, "species_a.tree"))
        sp_b = _tree(os.path.join(ALE_TOPO, "species_b.tree"))
        rec = os.path.join(ALE_TOPO, "fam.uml_rec")
        cache = str(tmp_path / "ale_cache")
        a = convert_from_ale(sp_a, [rec], min_family_size=0, cache_dir=cache, quiet=True)
        b = convert_from_ale(sp_b, [rec], min_family_size=0, cache_dir=cache, quiet=True)
        assert list(a.informative_edges()) == ["65,62"]  # parent(61) = 65
        assert list(b.informative_edges()) == ["69,62"]  # parent(61) = 69
        assert len(os.listdir(cache)) == 2  # 两棵树各自的条目

    def test_legacy_cache_entries_are_not_read(self, tmp_path):
        """无拓扑、无版本号前缀的缓存条目必须完全不被命中。"""
        import hashlib
        import pickle

        from maxtic_next.constraints.adapters.ale import _process_family_worker

        sp_a = _tree(os.path.join(ALE_TOPO, "species_a.tree"))
        rec = os.path.join(ALE_TOPO, "fam.uml_rec")
        adapter = ALEAdapter(sp_a, min_family_size=0, source="trf")
        st = os.stat(rec)
        # 逐字复刻旧键：无 v2 前缀、只有标签集合摘要
        raw = (
            f"{os.path.abspath(rec)}|{st.st_mtime}|{st.st_size}|"
            f"{adapter.min_support}|{adapter.min_family_size}|"
            f"{species_labels_digest(adapter._get_parent_map().keys())}"
        )
        cache = tmp_path / "legacy_cache"
        cache.mkdir()
        with open(cache / (hashlib.sha256(raw.encode()).hexdigest() + ".pkl"), "wb") as fh:
            pickle.dump({"61,62": 999}, fh)  # 恶意/过期的载荷
        got = _process_family_worker(
            (
                rec,
                0,
                adapter._get_parent_map(),
                adapter.extant_species,
                str(cache),
                adapter.min_support,
            )
        )
        # 命中的不是旧条目：解析结果仍来自真实文件（计数 1，而非 999）
        assert got[1].get("65,62") == 1
        assert got[5]["transfers_seen"] == 1


# ======================================================================
#  /  —— 距离列缺失与未知标签距离
# ======================================================================
class TestM24AndM10Distance:
    def test_distance_from_unknown_label_is_none_not_zero(self):
        sp = _tree(os.path.join(RANGER, "species.tree"))
        assert distance_from(sp, "61", "62") == 2.0
        assert distance_from(sp, "61", "NOT_A_LABEL") is None
        assert distance_from(sp, "NOT_A_LABEL", "61") is None
        # 真距离 0（同一节点）与"未知"严格可分
        assert distance_from(sp, "61", "61") == 0.0

    def test_filter_by_distance_reports_ignored_option(self):
        cset = ConstraintSet(
            [
                Constraint("61", "62", 1.0, metadata={}),  # 无距离列
                Constraint("62", "65", 2.0, metadata={}),
            ]
        )
        report = cset.filter_by_distance(5.0)
        assert report["ignored"] is True
        assert report["dropped"] == 0
        assert report["without_distance_column"] == 2
        assert report["kept"] == 2
        assert len(cset.constraints) == 2  # 与原版一致：全部保留
        assert cset.diagnostics["distance_filter"]["ignored"] is True

    def test_filter_by_distance_with_distance_column_actually_filters(self):
        cset = ConstraintSet(
            [
                Constraint("61", "62", 1.0, metadata={"distance": 1.0}),
                Constraint("62", "65", 1.0, metadata={"distance": 9.0}),
                Constraint("59", "61", 1.0, metadata={"distance": None}),
            ]
        )
        report = cset.filter_by_distance(5.0)
        assert report["ignored"] is False
        assert report["dropped"] == 1
        assert report["without_distance_column"] == 1
        assert {c.to_edge_key() for c in cset} == {"62,65", "59,61"}

    def test_self_loop_dropped_count_is_visible(self):
        cset = ConstraintSet([Constraint("61", "61", 1.0, metadata={"distance": 0.0})])
        report = cset.filter_by_distance(0)
        assert report["dropped"] == 1
        assert report["self_loops_dropped"] == 1


# ======================================================================
#  —— XML 解析失败不再被吞成"无约束"
# ======================================================================
class TestM25XmlErrors:
    def test_convert_raises_parse_error_with_file_and_position(self, tmp_path):
        sp = _tree(os.path.join(ECCE_NUM, "species.tree"))
        bad = _write(tmp_path / "broken.recphyloxml", "<recPhylo><recGeneTree><clade>\n")
        with pytest.raises(UpstreamParseError) as exc:
            convert_from_eccetera(sp, [bad], min_family_size=0, quiet=True)
        msg = str(exc.value)
        assert "broken.recphyloxml" in msg
        assert "第" in msg and "行" in msg
        assert exc.value.line >= 1
        assert exc.value.cause is not None  # 原始异常随身携带

    def test_strict_parse_function_raises_legacy_warns(self, tmp_path, capsys):
        text = "<recPhylo><clade>"
        with pytest.raises(UpstreamParseError):
            parse_recphyloxml(text, {"61"}, strict=True, source="x.recphyloxml")
        assert parse_recphyloxml(text, {"61"}) == ({}, 0, 0)
        err = capsys.readouterr().err
        assert "解析失败" in err and "不等于" in err  # 旧兼容路径也必须响亮

    def test_doctype_rejected_is_loud_in_adapter_path(self, tmp_path):
        sp = _tree(os.path.join(ECCE_NUM, "species.tree"))
        bad = _write(
            tmp_path / "entity.recphyloxml",
            '<!DOCTYPE recPhylo [<!ENTITY x "y">]><recPhylo><recGeneTree></recGeneTree></recPhylo>',
        )
        with pytest.raises(UpstreamParseError):
            convert_from_eccetera(sp, [bad], min_family_size=0, quiet=True)


# ======================================================================
#  /  —— 每一次静默 continue 都被计数
# ======================================================================
class TestM11DropAccounting:
    def test_ranger_diagnostics_count_every_drop(self, tmp_path):
        sp = _tree(os.path.join(RANGER, "species.tree"))
        f = _write(
            tmp_path / "mixed.dtl",
            "\n".join(
                [
                    "Reconciliation for Gene Tree 1:",
                    "g1: Leaf Node",
                    "g2: Leaf Node",
                    "g3: Leaf Node",
                    "m1 = LCA[a,b]: Transfer, Mapping --> 61, Recipient --> 62",
                    "m2 = LCA[c,d]: Transfer, Mapping --> H117_0, Recipient --> 62",
                    "m3 = LCA[e,f]: Transfer, Mapping --> 62, Recipient --> QQQ",
                    "m4 = LCA[g,h]: Transfer, Mapping --> 65, Recipient --> 65",
                ]
            )
            + "\n",
        )
        cset = convert_from_ranger_dtl(sp, [f], min_family_size=0, quiet=True)
        diag = cset.diagnostics
        assert diag["transfers_seen"] == 4
        assert diag["transfers_resolved"] == 2
        assert diag["dropped_by_label_miss"] == 2
        assert diag["self_loop_transfers"] == 1
        assert diag["constraints_kept"] == 2
        assert diag["edges_kept"] == 2
        assert "H117_0" in diag["unresolved_endpoint_samples"]

    def test_stderr_summary_line_is_printed(self, tmp_path, capsys):
        sp = _tree(os.path.join(RANGER, "species.tree"))
        convert_from_ranger_dtl(sp, [os.path.join(RANGER, "FAM1.dtl")], min_family_size=0)
        err = capsys.readouterr().err
        for token in (
            "donor 端点约定=",
            "权重口径=",
            "转移=",
            "保留=",
            "标签不匹配=",
            "支持度=",
            "家族规模=",
            "自环=",
        ):
            assert token in err

    def test_alerax_counts_malformed_and_label_misses(self, tmp_path):
        sp = _tree(os.path.join(DATA, "alerax_example", "species.tree"))
        f = _write(
            tmp_path / "FAMX_transfers.txt", "61 62 0.9\n62 65 notanumber\nZZZ 59 0.5\nshort\n"
        )
        cset = registry.convert("alerax", sp, [f], quiet=True)
        diag = cset.diagnostics
        assert diag["transfers_seen"] == 3
        assert diag["malformed_lines"] == 2
        assert diag["dropped_by_label_miss"] == 1
        assert diag["constraints_kept"] == 1


# ======================================================================
#  / --from-auto —— 端点约定机器可读 + 混用告警
# ======================================================================
class TestM10Conventions:
    def test_registry_exposes_native_conventions(self):
        assert registry.endpoint_convention("ale") == CONVENTION_PARENT_OF_DONOR
        for tool in ("ranger", "artra", "eccetera", "alerax"):
            assert registry.endpoint_convention(tool) == CONVENTION_DONOR_ITSELF

    def test_each_adapter_reports_its_own_convention(self):
        sp = _tree(os.path.join(RANGER, "species.tree"))
        ale_sp = _tree(os.path.join(ALE_DIR, "species.tree"))
        ale = convert_from_ale(
            ale_sp, [os.path.join(ALE_DIR, "rec.uml_rec")], min_family_size=0, quiet=True
        )
        ran = convert_from_ranger_dtl(
            sp, [os.path.join(RANGER, "FAM1.dtl")], min_family_size=0, quiet=True
        )
        ecc = convert_from_eccetera(
            _tree(os.path.join(ECCE_NUM, "species.tree")),
            [os.path.join(ECCE_NUM, "FAM1.recphyloxml")],
            min_family_size=0,
            quiet=True,
        )
        assert ale.diagnostics["donor_endpoint_convention"] == CONVENTION_PARENT_OF_DONOR
        assert ran.diagnostics["donor_endpoint_convention"] == CONVENTION_DONOR_ITSELF
        assert ecc.diagnostics["donor_endpoint_convention"] == CONVENTION_DONOR_ITSELF
        # 每条约束也带口径元数据，聚合后可逐条追溯
        assert all(c.metadata["weight_semantics"] for c in ale)

    def test_merging_mixed_conventions_warns(self, capsys):
        sp = _tree(os.path.join(RANGER, "species.tree"))
        ale_sp = _tree(os.path.join(ALE_DIR, "species.tree"))
        ale = convert_from_ale(
            ale_sp, [os.path.join(ALE_DIR, "rec.uml_rec")], min_family_size=0, quiet=True
        )
        ran = convert_from_ranger_dtl(
            sp, [os.path.join(RANGER, "FAM1.dtl")], min_family_size=0, quiet=True
        )
        merged = registry.merge_constraint_sets([ale, ran])
        err = capsys.readouterr().err
        assert "不同的供体端点层级约定" in err
        assert merged.diagnostics["convention_conflict"] == [
            CONVENTION_DONOR_ITSELF,
            CONVENTION_PARENT_OF_DONOR,
        ]
        assert merged.diagnostics["donor_endpoint_convention"] == "mixed"
        assert len(merged.constraints) == len(ale.constraints) + len(ran.constraints)

    def test_same_convention_merge_is_silent(self, capsys):
        sp = _tree(os.path.join(RANGER, "species.tree"))
        a = convert_from_ranger_dtl(
            sp, [os.path.join(RANGER, "FAM1.dtl")], min_family_size=0, quiet=True
        )
        b = convert_from_ranger_dtl(
            sp, [os.path.join(RANGER, "FAM1.dtl")], min_family_size=0, quiet=True
        )
        merged = registry.merge_constraint_sets([a, b])
        assert "convention_conflict" not in merged.diagnostics
        assert capsys.readouterr().err == ""

    def test_convert_auto_allow_mixed_groups_and_warns(self, capsys):
        sp = _tree(os.path.join(RANGER, "species.tree"))
        inputs = [os.path.join(ALE_DIR, "rec.uml_rec"), os.path.join(RANGER, "FAM1.dtl")]
        # 默认：混合格式仍然抛错（保持既有 CLI 语义）
        with pytest.raises(ValueError):
            registry.convert_auto(sp, inputs)
        merged = registry.convert_auto(sp, inputs, allow_mixed=True, min_family_size=0, quiet=True)
        assert "不同的供体端点层级约定" in capsys.readouterr().err
        assert len(merged.constraints) >= 3
        assert {m["tool"] for m in merged.diagnostics["merged_from"]} == {"ale", "ranger"}
