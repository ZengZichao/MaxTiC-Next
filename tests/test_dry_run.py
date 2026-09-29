"""静态预检单元测试。"""

import os

from maxtic_next.dry_run import dry_run_check, format_dry_run
from maxtic_next.io.parsing import read_newick_file

DATA_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")


def _tree():
    return read_newick_file(os.path.join(DATA_DIR, "minitree.tree"))


def test_dry_run_clean():
    """合法输入：预检通过，无错误。"""
    tree = _tree()
    rep = dry_run_check(tree, [os.path.join(DATA_DIR, "Cyano_CUTConstraints.tsv")])
    assert rep["ok"] is True
    assert rep["n_errors"] == 0
    assert isinstance(format_dry_run(rep), str)


def test_dry_run_missing_constraint_file():
    """约束文件不存在：报告致命错误。"""
    tree = _tree()
    rep = dry_run_check(tree, ["/no/such/file.tsv"])
    assert rep["ok"] is False
    assert rep["n_errors"] >= 1
    assert any("不存在" in it["message"] for it in rep["issues"])


def test_dry_run_distance_semantics_is_error_and_truthful():
    """d>0 但约束文件不含 phylogenetic distance 列：报告 error，且如实说明"被忽略"。

    ：旧文案称"按距离过滤将丢弃全部约束"，与实现（无距离列者一律
    保留，即该选项被忽略）恰好相反。
    """
    tree = _tree()
    # Cyano 文件为 3 列双格式（无 distance 列）
    rep = dry_run_check(
        tree, [os.path.join(DATA_DIR, "Cyano_CUTConstraints.tsv")], min_transfer_distance=1.0
    )
    hits = [it for it in rep["issues"] if "min-transfer-distance" in it["message"]]
    assert hits, rep["issues"]
    assert hits[0]["severity"] == "error"
    assert "忽略" in hits[0]["message"]
    assert "丢弃全部" not in hits[0]["message"]
    assert rep["ok"] is False


def test_dry_run_reports_non_finite_and_negative_weights_with_line():
    """：nan / inf / 负权重在预检即报 error，且信息含文件名与行号。"""
    import tempfile

    tree = _tree()
    path = os.path.join(tempfile.gettempdir(), "bad_weights.tsv")
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("61 62 1.0\n61 67 nan\n")
    try:
        rep = dry_run_check(tree, [path])
        assert rep["ok"] is False
        msgs = " ".join(it["message"] for it in rep["issues"])
        assert "nan" in msgs and "第 2 行" in msgs and "bad_weights.tsv" in msgs
    finally:
        os.remove(path)


def test_dry_run_ale_missing_file():
    """--from-ale 模式下 ``.uml_rec`` 文件不存在：报告致命错误。"""
    tree = _tree()
    rep = dry_run_check(tree, ["/no/such.rec"], from_ale=True)
    assert rep["ok"] is False
    assert any("uml_rec" in it["message"] for it in rep["issues"])


def test_dry_run_bad_endpoint():
    """约束端点不在物种树中：报告致命错误。"""
    import tempfile

    tree = _tree()
    path = os.path.join(tempfile.gettempdir(), "bad_endpoint.tsv")
    with open(path, "w", encoding="utf-8") as fh:
        # 约束文件为制表符分隔三列格式；999 不是合法标签
        fh.write("61\t999\t10.0\n")
    try:
        rep = dry_run_check(tree, [path])
        assert rep["ok"] is False
        assert any("不在物种树中" in it["message"] for it in rep["issues"])
    finally:
        os.remove(path)


def test_dry_run_polytomy_rejected():
    """物种树存在多歧分支（内部节点子节点数 != 2）应预检报错（二叉树强假设）。"""
    import tempfile
    from maxtic_next.io.parsing import read_newick_file

    tree_str = "((A,B,C)X,Y)Z;"
    path = os.path.join(tempfile.gettempdir(), "polytomy.tree")
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(tree_str)
    try:
        tree = read_newick_file(path)
        rep = dry_run_check(tree, [os.path.join(DATA_DIR, "Cyano_CUTConstraints.tsv")])
        assert rep["ok"] is False
        assert any("多歧" in it["message"] or "polytomy" in it["message"] for it in rep["issues"])
    finally:
        os.remove(path)


def test_dry_run_invalid_target_clade_rejected():
    """--target-clade 指定不存在的标签应在预检阶段报错（不等到正式运行 KeyError）。"""
    tree = _tree()
    rep = dry_run_check(
        tree, [os.path.join(DATA_DIR, "Cyano_CUTConstraints.tsv")], target_clade="999"
    )
    assert rep["ok"] is False
    assert any("target-clade" in it["message"] for it in rep["issues"])
