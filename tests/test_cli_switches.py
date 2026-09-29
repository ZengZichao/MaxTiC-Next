"""CLI / API 开关冒烟测试（Task #5）。

覆盖 ``--rd`` / ``--d`` / ``--ts`` / ``--random-type`` / ``--from-ale``：
断言不崩溃且输出合理、可复现。

说明：与原版 Python 2.7 ``MaxTiC.py`` 黄金值的逐字节比对需要在容器化
Python 2.7 环境进行（本机无 Py2.7），此处仅做可复现的合理性校验。
"""

import os

from maxtic_next.api import rank

DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")
TREE = os.path.join(DATA, "minitree.tree")
CONS = os.path.join(DATA, "Cyano_CUTConstraints.tsv")
ALE_TREE = os.path.join(DATA, "ale_example", "species.tree")
ALE_REC = os.path.join(DATA, "ale_example", "rec.uml_rec")


def _write_space_constraints(path):
    """写一个空格格式约束文件，含第 4 列 phylogenetic distance。

    61->67 的 distance=0.5 应被 ``--d 3`` 过滤丢弃；
    61->62 的 distance=5.0 与无距离列的 62->45 应保留。
    """
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("61 62 10.0 5.0\n")
        fh.write("61 67 20.0 0.5\n")
        fh.write("62 45 30.0\n")


def test_output_style_default_short_then_legacy(tmp_path):
    """端到端：默认产出简洁短名；output_style=legacy 回到原版长名；两者内容一致。"""
    p_short = str(tmp_path / "short")
    r_short = rank(TREE, CONS, seed=42, output_prefix=p_short, print_summary=False)
    assert r_short.informative_file.endswith(".mt.informative.tsv")
    assert r_short.conflicting_file.endswith(".mt.conflicts.tsv")
    assert r_short.partial_order_file.endswith(".mt.partial_order.tsv")
    assert os.path.exists(r_short.informative_file)

    p_legacy = str(tmp_path / "legacy")
    r_legacy = rank(
        TREE, CONS, seed=42, output_prefix=p_legacy, print_summary=False, output_style="legacy"
    )
    assert r_legacy.informative_file.endswith(
        "_MT_output_filtered_list_of_weighted_informative_constraints"
    )
    # 内容一致（仅文件名不同）
    assert open(r_short.informative_file).read() == open(r_legacy.informative_file).read()


def test_rd_pvalue_off_by_default():
    """默认 random_trees=0 时不计算 p 值，保持 P0 行为不变。"""
    res = rank(TREE, CONS, seed=42, print_summary=False)
    assert res.random_stats is None


def test_rd_pvalue_computed(capsys):
    """--rd 计算随机树 p 值，stdout 含原版 4 行打印，且同 seed 可复现。"""
    res = rank(TREE, CONS, seed=42, random_trees=15, print_summary=True)
    captured = capsys.readouterr().out

    assert res.random_stats is not None
    rs = res.random_stats
    assert rs["n"] == 15
    assert 0.0 <= rs["value_pvalue"] <= 1.0
    assert 0.0 <= rs["sim_pvalue"] <= 1.0
    assert "pvalue of the found order:" in captured
    assert "pvalue of the similarity with the input order:" in captured

    # 确定性：同 seed 复现相同 p 值
    res2 = rank(TREE, CONS, seed=42, random_trees=15, print_summary=False)
    assert res2.random_stats["value_pvalue"] == rs["value_pvalue"]
    assert res2.random_stats["sim_pvalue"] == rs["sim_pvalue"]


def test_min_transfer_distance_space(tmp_path):
    """--d 对空格格式第 4 列 distance 生效：d=3 丢弃 61->67。"""
    cons = str(tmp_path / "cons.txt")
    _write_space_constraints(cons)
    r0 = rank(TREE, cons, seed=42, min_transfer_distance=0.0, print_summary=False)
    r3 = rank(TREE, cons, seed=42, min_transfer_distance=3.0, print_summary=False)
    # 被过滤的 61->67（weight 20，distance 0.5 <= 3）不再进入边集合，
    # 故总权重下降；其余（61->62 distance 5、62->45 无距离）保留。
    assert r3.total_weight < r0.total_weight


def test_min_transfer_distance_drops_low_distance_constraint(tmp_path):
    """--d 端到端：d=3 时 distance=0.5 的 61->67 被过滤出输出，d=0 时保留。"""
    cons = str(tmp_path / "cons.txt")
    _write_space_constraints(cons)
    prefix0 = str(tmp_path / "r0")
    prefix3 = str(tmp_path / "r3")
    r0 = rank(
        TREE, cons, seed=42, min_transfer_distance=0.0, output_prefix=prefix0, print_summary=False
    )
    r3 = rank(
        TREE, cons, seed=42, min_transfer_distance=3.0, output_prefix=prefix3, print_summary=False
    )

    def _out_keys(res):
        # 信息性 / 冲突行格式为 "key weight"，key = "donor,receptor"
        return {line.split()[0] for line in res.informative_lines + res.conflicting_lines}

    # d=0：61->67（distance 0.5）进入输出
    assert "61,67" in _out_keys(r0)
    # d=3：被过滤，不再出现在任何输出边
    assert "61,67" not in _out_keys(r3)
    # 高 distance 的 61->62（distance 5.0）保留
    assert "61,62" in _out_keys(r3)


def test_threshold_constraints(tmp_path):
    """--ts 按权重阈值移除最低权重边，信息性边数量不增。"""
    cons = str(tmp_path / "cons.txt")
    _write_space_constraints(cons)
    r0 = rank(TREE, cons, seed=42, threshold_constraints=0.0, print_summary=False)
    r5 = rank(TREE, cons, seed=42, threshold_constraints=0.5, print_summary=False)
    assert len(r5.informative_lines) <= len(r0.informative_lines)


def test_random_type_deterministic():
    """--random-type 0/1/2 均不崩溃，且同 seed 可复现。"""
    for rt in (0, 1, 2):
        r = rank(TREE, CONS, seed=42, random_type=rt, print_summary=False)
        assert r.best_order  # 非空合法排序
    a = rank(TREE, CONS, seed=7, random_type=1, print_summary=False).best_order
    b = rank(TREE, CONS, seed=7, random_type=1, print_summary=False).best_order
    assert a == b


def test_from_ale_two_stage(tmp_path):
    """--from-ale 两阶段：仅转换约束写出 5 列 ALE 格式后提前返回。"""
    out = str(tmp_path / "out_cons.tsv")
    res = rank(ALE_TREE, ALE_REC, from_ale=True, constraints_out=out, print_summary=False)
    assert os.path.exists(out)
    with open(out, encoding="utf-8") as fh:
        content = fh.read().strip()
    assert content  # 非空
    # 5 列 ALE 逗号格式：family,donor,receptor,weight,distance
    first = content.splitlines()[0].split(",")
    assert len(first) == 5
    # 两阶段模式提前返回，不进入正式排序
    assert res.best_order == []
