"""自环转移（``donor == receptor``）与 ``--d`` 距离过滤的回归测试。

背景
----
供体与受体是同一谱系（``donor == receptor``）的"自环转移"在物种树上的拓扑距离
恒为 ``0``。若适配器把这个 0 写进 ``metadata["distance"]``，下游
``ConstraintSet.filter_by_distance(d)`` 的保留规则是 ``distance > d``，于是
**连默认阈值 ``d = 0``（``config.DEFAULT_MIN_TRANSFER_DIST = 0``）都会把它判掉**，
这类事件便从原版口径的 stdout 统计 ``to itself`` 里静默消失。

RANGER-DTLx / ecceTERA / ARTra / AleRax / ALE 五个适配器对自环一律写
``distance = None``（= 无距离信息，与"文本输入没有第 4/5 列"同义，一律保留）。
本文件锁定 ALE 侧的这条不变量，确保它
（``constraints/adapters/ale.py`` 的 ``ALEAdapter._merge_payload``，全仓库唯一构造
ALE ``Constraint`` 的地方，顺序与进程池并行两条路径共用）不再退回 ``0.0``。

为什么在适配器层 + ``api.rank`` 两层都测
----------------------------------------
* 适配器层（用例 1/2/3）直接盯住不变量本身：自环的 ``distance is None``。
* 端到端层（用例 4/5）用**默认** ``min_transfer_distance`` 跑一次真实排序，断言
  ``Result.uninformative["to_itself"]`` 非零且 stdout 摘要里有 ``1.0 to itself``，
  并且不出现"自环约束被丢弃"的 NOTE——这才是用户看得见的症状；用例 5 另外锁住
  两阶段产物（``-o`` 的 5 列逗号格式）第 5 列留空、读回后端点一致。
* 用例 3 是**判别力自检**：把同一条自环的距离人为改回 ``0.0``，确认它确实会被
  默认阈值删掉。没有这条对照，用例 1/2 的通过毫无信息量（不变量被破坏时无人报警）。

全部产物写入 ``tmp_path``，``api.rank`` 一律 ``force=True``，不改动仓库任何文件。
"""

from __future__ import annotations

import pytest

from maxtic_next import api
from maxtic_next.config import DEFAULT_MIN_TRANSFER_DIST
from maxtic_next.constraints.adapters.ale import ALEAdapter
from maxtic_next.tree.tree import Tree

# 与 tests/data/ale_example/species.tree（= examples/adapters/ale/species.tree）
# 同一棵树：内部标签 59/61/62/65/69，parent(61)=parent(62)=65、parent(65)=69。
SPECIES_NEWICK = "((CYAP8,CYAP0)59,((CYAA5,CYAP2)61,(NOSP7,ANAVT)62)65)69;\n"

# 一份含**自环转移**的 ALE ``.uml_rec``（默认口径 ``trf``，donor = 供体的父）：
#   * ``T@62->65`` → 65,65 —— parent(62) == 受体 65 ⇒ **自环**（供体 = 受体谱系）
#   * ``T@61->59`` → 65,59 —— 两个不同子树，普通信息性约束（拓扑距离 2.0）
# 头行声明 ``1 reconciled``，故权重口径为支持度（1/1 = 1.0 > min_support 0.05）；
# 基因树有 6 个叶子，故默认的 ``min_family_size=5``（判据 len(leaves) > 5）不会
# 把整个家族跳过 —— 本用例因此**不需要**任何放宽阈值的开关。
REC_NEWICK = (
    "1 reconciled\n"
    "FAM_SELFLOOP\n"
    "(((CYAP8:1.0,CYAP0:1.0)59:1.0,"
    "(CYAA5:1.0,CYAP2:1.0)61.T@61->59:1.0)65:1.0,"
    "(NOSP7:1.0,ANAVT:1.0)62.T@62->65:1.0)69:1.0;\n"
)


def _write(tmp_path, name, text):
    path = tmp_path / name
    path.write_text(text, encoding="utf-8")
    return str(path)


def _species_tree(tmp_path):
    path = _write(tmp_path, "species.tree", SPECIES_NEWICK)
    tree = Tree()
    with open(path, encoding="utf-8") as fh:
        tree.read_newick(fh.readline())
    return tree


def _rec_file(tmp_path, name="selfloop.uml_rec"):
    return _write(tmp_path, name, REC_NEWICK)


def _convert(tmp_path, parallel="process", **kwargs):
    kwargs.setdefault("quiet", True)
    adapter = ALEAdapter(_species_tree(tmp_path), **kwargs)
    return adapter.convert([_rec_file(tmp_path)])


def _by_key(cset):
    return {c.to_edge_key(): c for c in cset.constraints}


# ----------------------------------------------------------------------
# 1. 自环转移的 distance 必须是 None（写成 0.0 就会被过滤掉）
# ----------------------------------------------------------------------
def test_ale_self_loop_distance_is_none(tmp_path):
    cset = _convert(tmp_path)
    self_loop = _by_key(cset)["65,65"]
    assert self_loop.donor == self_loop.receptor

    # 这一行就是  的不变量：None = 无距离信息，而不是"真距离 0"。
    assert self_loop.metadata["distance"] is None
    # 对照：普通约束仍带真实拓扑距离（65 → 59 相差两代）。
    assert _by_key(cset)["65,59"].metadata["distance"] == 2.0
    # 自环被计数，但**不是**被丢弃。
    assert cset.diagnostics["self_loop_transfers"] == 1
    assert cset.diagnostics["dropped_by_label_miss"] == 0


# ----------------------------------------------------------------------
# 2. 默认 --min-transfer-distance 下自环保留（全仓库唯一构造点 ⇒ 顺序=并行）
# ----------------------------------------------------------------------
@pytest.mark.parametrize("parallel", ["process", "thread", "none"])
def test_ale_self_loop_survives_default_distance_filter(tmp_path, parallel):
    cset = _convert(tmp_path, parallel=parallel)
    assert len(cset) == 2

    report = cset.filter_by_distance(DEFAULT_MIN_TRANSFER_DIST)
    assert report["threshold"] == 0
    assert report["dropped"] == 0
    assert report["self_loops_dropped"] == 0
    assert report["without_distance_column"] == 1  # 仅自环无距离列，普通约束带 2.0
    assert report["kept"] == 2
    assert "65,65" in cset.informative_edges()
    assert cset.informative_edges()["65,65"] == pytest.approx(1.0)


# ----------------------------------------------------------------------
# 3. 判别力自检：距离写回 0.0 时，默认阈值确实会吃掉自环
# ----------------------------------------------------------------------
def test_zero_distance_self_loop_would_be_dropped_by_default_threshold(tmp_path):
    cset = _convert(tmp_path)
    self_loop = _by_key(cset)["65,65"]
    assert self_loop.metadata["distance"] is None

    self_loop.metadata["distance"] = 0.0  # 人为写入 0.0，复现被吃掉的症状
    report = cset.filter_by_distance(DEFAULT_MIN_TRANSFER_DIST)
    assert report["dropped"] == 1
    assert report["self_loops_dropped"] == 1
    assert "65,65" not in cset.informative_edges()  # ⇒ to itself 恒为 0


# ----------------------------------------------------------------------
# 4. 端到端：默认阈值下排序结果的 to itself 非零，且无"自环被丢弃"NOTE
# ----------------------------------------------------------------------
def test_rank_keeps_ale_self_loop_in_to_itself_statistic(tmp_path, capsys):
    tree = _write(tmp_path, "species.tree", SPECIES_NEWICK)
    prefix = str(tmp_path / "ale_selfloop")

    res = api.rank(
        tree,
        [_rec_file(tmp_path)],
        from_tool="ale",
        output_prefix=prefix,
        force=True,
        seed=42,
        html_report=False,
        print_summary=True,
        adapter_quiet=True,
    )

    # 默认 min_transfer_distance 未被改动 ⇒ 自环必须活着进到统计里
    assert res.uninformative["to_itself"] == pytest.approx(1.0)
    assert res.uninformative["total"] == pytest.approx(1.0)
    assert res.total_weight == pytest.approx(2.0)

    out = capsys.readouterr().out
    assert "1.0 to itself" in out
    assert "0.0 to itself" not in out
    assert not any("自环约束" in w for w in res.warnings), res.warnings


def test_self_loop_to_itself_is_independent_of_d_threshold(tmp_path):
    """``--d 0`` 与 ``--d 1`` 的 ``to itself`` 完全相同：自环穿过了距离过滤。

    ``65 → 59`` 的拓扑距离 2.0 在两种阈值下都保留，故两次运行的约束集一致。
    """
    tree = _write(tmp_path, "species.tree", SPECIES_NEWICK)
    weights = {}
    for d in (0, 1):
        res = api.rank(
            tree,
            [_rec_file(tmp_path, name=f"selfloop_{d}.uml_rec")],
            from_tool="ale",
            output_prefix=str(tmp_path / f"ale_d{d}"),
            force=True,
            seed=42,
            min_transfer_distance=d,
            print_summary=False,
            html_report=False,
            adapter_quiet=True,
        )
        weights[d] = res.uninformative["to_itself"]
        assert res.total_weight == pytest.approx(2.0)
    assert weights[0] == pytest.approx(1.0)
    assert weights[0] == pytest.approx(weights[1])


# ----------------------------------------------------------------------
# 5. 两阶段产物（-o 的 5 列 ALE 逗号格式）第 5 列留空，读回后症状一致
# ----------------------------------------------------------------------
def test_constraints_out_leaves_distance_column_empty(tmp_path):
    tree = _write(tmp_path, "species.tree", SPECIES_NEWICK)
    out = str(tmp_path / "unified.tsv")
    api.rank(
        tree,
        [_rec_file(tmp_path)],
        from_tool="ale",
        constraints_out=out,
        output_prefix=str(tmp_path / "ale_two_stage"),
        force=True,
        print_summary=False,
        html_report=False,
        adapter_quiet=True,
    )

    with open(out, encoding="utf-8") as fh:
        lines = [ln.strip() for ln in fh if ln.strip() and not ln.startswith("#")]
    assert "selfloop,65,65,1.0," in lines  # 第 5 列为空 = None
    assert "selfloop,65,59,1.0,2.0" in lines  # 普通约束仍写真实距离

    back = api.rank(
        tree,
        out,
        output_prefix=str(tmp_path / "ale_readback"),
        force=True,
        seed=42,
        print_summary=False,
        html_report=False,
    )
    assert back.uninformative["to_itself"] == pytest.approx(1.0)
