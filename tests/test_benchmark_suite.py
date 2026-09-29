"""基准套件的计时口径、重复结构、规模计数与正确性断言。

覆盖 :mod:`maxtic_next.benchmarks.suite`：

* CPU 用 ``time.process_time``、墙钟单独报告（**行为可验证**：把 ``process_time``
  冻成常数后 CPU 样本必须为 0，而墙钟仍为正）；
* 输出 I/O 与 ``tracemalloc`` 都不在计时窗口内；
* 默认 1 次预热 + 3 次重复，报告中位数与 IQR，而不是单次测量；
* 合成树确实受 ``seed`` 驱动；约束**无放回**生成、口径如实上报；
* 数据集配对不再静默丢弃未匹配文件；``#`` / FRQ 行不算约束；
* 除"跑完了没报错"外，还对交付值做**独立重算**的正确性断言。
"""

import os
import tracemalloc

import pytest

from maxtic_next.benchmarks import suite
from maxtic_next.tree.tree import Tree

DATA_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")
TREE = os.path.join(DATA_DIR, "minitree.tree")
CONS = os.path.join(DATA_DIR, "Cyano_CUTConstraints.tsv")


# =========================================================================
# 重复结构与统计量
# =========================================================================
def test_run_benchmark_reports_repeat_structure():
    res = suite.run_benchmark(
        n_internal=6, n_edges=12, seed=1, time_budget=30.0, n_repeats=3, warmup=2
    )
    assert res["error"] is None, res["error"]
    assert res["n_repeats"] == 3 and res["warmup"] == 2
    assert len(res["cpu_time_samples"]) == 3
    assert len(res["wall_time_samples"]) == 3
    assert res["cpu_time_median"] == pytest.approx(sorted(res["cpu_time_samples"])[1])
    assert res["cpu_time_iqr"] is not None
    assert res["wall_time_iqr"] is not None
    assert res["cpu_time"] == res["cpu_time_median"]
    assert res["wall_time"] == res["wall_time_median"]
    # I/O 被单列，且内核时间 = 总 CPU - I/O（不为负）
    assert res["io_time"] is not None and res["io_time"] >= 0.0
    assert 0.0 <= res["kernel_cpu_time"] <= res["cpu_time"] + 1e-12
    assert res["peak_memory_mb"] >= 0.0


def test_run_suite_has_one_row_per_grid_point():
    rows = suite.run_suite([3, 5], [4, 8], seed=2, time_budget=30.0, n_repeats=2, warmup=0)
    assert len(rows) == 4
    for r in rows:
        assert r["error"] is None
        assert len(r["cpu_time_samples"]) == 2


# =========================================================================
# CPU 时间名副其实：process_time，而不是墙钟
# =========================================================================
def test_cpu_time_uses_process_time_not_wall_clock(monkeypatch):
    real = suite.time.process_time
    monkeypatch.setattr(suite.time, "process_time", lambda: 0.0)  # 冻结 CPU 钟
    res = suite.run_benchmark(n_internal=6, n_edges=10, seed=3, n_repeats=3, warmup=0)
    assert res["error"] is None
    for v in res["cpu_time_samples"]:
        assert v == pytest.approx(0.0, abs=1e-9)
    assert res["wall_time"] > 0.0  # 墙钟仍在走
    assert res["kernel_cpu_time"] == pytest.approx(0.0, abs=1e-9)
    monkeypatch.undo()
    assert suite.time.process_time is real


def test_tracemalloc_not_inside_timed_window(monkeypatch):
    tracing_at_rank_call = []
    original_make_ranker = suite._make_ranker

    def spy(tree, cset, seed, output_prefix, **kw):
        tracing_at_rank_call.append(tracemalloc.is_tracing())
        return original_make_ranker(tree, cset, seed, output_prefix, **kw)

    starts = {"n": 0}
    real_start = suite.tracemalloc.start

    def counting_start(*a, **kw):
        starts["n"] += 1
        return real_start(*a, **kw)

    monkeypatch.setattr(suite, "_make_ranker", spy)
    monkeypatch.setattr(suite.tracemalloc, "start", counting_start)
    res = suite.run_benchmark(n_internal=5, n_edges=8, seed=4, n_repeats=3, warmup=1)
    assert res["error"] is None
    # 3 次计时 + 1 次预热调用时都不得处于 tracemalloc 之下；
    # 只有独立的那趟内存测量在 tracing 中
    assert tracemalloc.is_tracing() is False
    assert any(v is False for v in tracing_at_rank_call)
    assert starts["n"] == 1, "峰值内存应由单独一趟测量"
    assert len(tracing_at_rank_call) == 3 + 1 + 1  # 预热 + 3 重复 + 内存趟
    assert res["peak_memory_mb"] > 0.0


# =========================================================================
# 正确性断言（不是只有计时）
# =========================================================================
def test_run_benchmark_includes_correctness_assertion():
    res = suite.run_benchmark(n_internal=7, n_edges=15, seed=5, n_repeats=1, warmup=0)
    assert res["correctness_ok"] is True
    checks = res["correctness"]
    assert checks["matches_conflicting_sum"] is True
    assert checks["not_worse_than_heuristics"] is True
    assert checks["finite_value"] is True
    # 独立重算的交付值必须不劣于任一启发式（择优/局部搜索没有把结果改坏）
    assert res["delivered_value"] <= min(res["best_value"], res["delivered_value"]) + 1e-9
    assert checks["informative_count"] > 0


def test_verify_delivered_value_catches_wrong_result():
    """基准必须能抓出"计时正常、但交付值不对"的运行。

    三个夹具都**内部自洽**（``conflicting_lines`` 与 ``best_order`` 的真实违反集
    一致）：``verify_delivered_value`` 是"从交付序独立重算并与上报值对账"，因此
    自相矛盾的夹具无法用于判定实现对错。
    """
    from types import SimpleNamespace

    # 交付序 ["a","c","b"] 违反 b→c（权重 1.0）：三点对账全部吻合
    good = SimpleNamespace(
        informative_lines=["a,b 2.0", "b,c 1.0"],
        conflicting_lines=["b,c 1.0"],
        best_order=["a", "c", "b"],
        values={"input": 3.0, "greedy": 1.0, "mixing": 2.0, "best": 1.0},
    )
    checks = suite.verify_delivered_value(good)
    assert checks["ok"] is True, checks
    assert checks["recomputed_value"] == pytest.approx(1.0)

    # 谎报：声称 value=0.0 且冲突文件为空，但交付序确实违反了 b→c（1.0）
    lying = SimpleNamespace(
        informative_lines=good.informative_lines,
        conflicting_lines=[],
        best_order=["a", "c", "b"],
        values={"input": 3.0, "greedy": 1.0, "mixing": 2.0, "best": 0.0},
    )
    checks = suite.verify_delivered_value(lying)
    assert checks["ok"] is False
    assert checks["matches_conflicting_sum"] is False
    assert checks["matches_reported_best"] is False
    assert checks["recomputed_value"] == pytest.approx(1.0)

    # 比启发式还差（择优 / 局部搜索把结果改坏）：交付序 3.0 > greedy 0.2
    worse = SimpleNamespace(
        informative_lines=good.informative_lines,
        conflicting_lines=["a,b 2.0", "b,c 1.0"],
        best_order=["c", "b", "a"],
        values={"input": 3.0, "greedy": 0.2, "mixing": 0.3, "best": 3.0},
    )
    checks = suite.verify_delivered_value(worse)
    assert checks["ok"] is False
    assert checks["not_worse_than_heuristics"] is False
    assert checks["matches_conflicting_sum"] is True


def test_verify_delivered_value_accepts_non_informative_violations():
    """被违反的约束**不止** informative 一类。

    "供体是受体后代"这类反向约束在任何合法序里都被违反，它进 ``conflicting_lines``
    与 ``value()``，但按分类**不**是信息性约束。旧实现只用 informative 重建边集，
    于是把完全正确的真实数据结果误判为 BAD。
    """
    from types import SimpleNamespace

    res = SimpleNamespace(
        informative_lines=["a,b 2.0"],  # 信息性：未被违反
        conflicting_lines=["a,b 2.0", "c,a 5.0"],  # c→a 为反向（平凡冲突）约束
        best_order=["b", "a", "c"],
        values={"input": 7.0, "greedy": 7.0, "mixing": 7.0, "best": 7.0},
    )
    checks = suite.verify_delivered_value(res)
    assert checks["ok"] is True, checks
    assert checks["recomputed_value"] == pytest.approx(7.0)
    assert checks["informative_count"] == 1


# =========================================================================
# 规模口径：seed 有效、无放回、注释行不算约束
# =========================================================================
def test_synthetic_tree_seed_is_effective():
    n1, l1 = suite.make_synthetic_tree(9, seed=1)
    n1b, _ = suite.make_synthetic_tree(9, seed=1)
    n2, l2 = suite.make_synthetic_tree(9, seed=2)
    assert n1 == n1b, "同一 seed 必须逐字符可复现"
    assert n1 != n2, "不同 seed 应给出不同的合成树"
    tree = Tree()
    tree.read_newick(n2)
    assert sorted(l1) == sorted(l2)  # 标签集合不变：仍是 1..9
    assert l1 != l2  # 但落位随 seed 变化
    assert len(tree.internal_node_labels()) == 9
    assert len(set(tree.get_leaves_names())) == 10  # 叶子数 = n + 1


def test_synthetic_constraints_without_replacement():
    _, labels = suite.make_synthetic_tree(6, seed=1)
    stats = {}
    cset = suite.make_synthetic_constraints(labels, 10, seed=1, stats_out=stats)
    pairs = [(c.donor, c.receptor) for c in cset.constraints]
    assert len(pairs) == len(set(pairs)), "不得有放回抽样后再靠聚合去重"
    assert stats["n_edges_generated"] == 10
    assert stats["n_edge_keys_effective"] == 10
    assert stats["n_pairs_available"] == 6 * 5


def test_synthetic_constraints_excludes_topology_forced_pairs():
    newick, labels = suite.make_synthetic_tree(7, seed=2)
    tree = Tree()
    tree.read_newick(newick)
    stats = {}
    cset = suite.make_synthetic_constraints(labels, 50, seed=2, tree=tree, stats_out=stats)
    anc = set()
    ids = {}
    for nid in tree.get_nodes():
        if not tree.is_leaf(nid):
            ids[tree.get_bootstrap(nid)] = nid
    for a in labels:
        for b in labels:
            if a != b and tree.is_ancestor(ids[a], ids[b]):
                anc.add((a, b))
    for c in cset.constraints:
        assert (c.donor, c.receptor) not in anc
        assert (c.receptor, c.donor) not in anc
    assert stats["n_pairs_available"] < 7 * 6


def test_synthetic_constraints_clamps_and_reports_honestly():
    _, labels = suite.make_synthetic_tree(3, seed=1)  # 只有 3*2=6 个有序对
    stats = {}
    cset = suite.make_synthetic_constraints(labels, 100, seed=1, stats_out=stats)
    assert len(cset.constraints) == stats["n_edges_generated"] == 6
    assert stats["n_edges_requested"] == 100
    assert stats["n_pairs_available"] == 6


def test_count_constraints_ignores_comments_and_frq(tmp_path):
    p = tmp_path / "c.tsv"
    p.write_text("# a comment line\n1\t2\t1.0\n\n2\t3\t0.5\nFRQ\tdonor\treceptor\n")
    assert suite._count_constraints(str(p)) == 2
    stats = suite._constraint_stats(str(p))
    assert stats["n_lines_nonempty"] == 4
    assert stats["n_constraints_parsed"] == 2
    assert stats["n_comment_or_header_lines"] == 2
    assert suite._count_constraints(str(tmp_path / "missing.tsv")) == 0


def test_count_constraints_on_real_test_data():
    assert suite._count_constraints(CONS) == 123


# =========================================================================
# 数据集配对不再静默丢弃
# =========================================================================
def test_find_dataset_pairs_keeps_unmatched_files(tmp_path):
    (tmp_path / "a.tree").write_text("((A:1,B:1)1:1,C:1)2;")
    (tmp_path / "a.tsv").write_text("1\t2\t1.0\n")
    (tmp_path / "b.tree").write_text("((A:1,B:1)1:1,C:1)2;")
    pairs = suite.find_dataset_pairs(str(tmp_path))
    # 旧实现：a 已同名配对 -> b.tree 被静默丢弃（返回 1 对）
    assert len(pairs) == 2, f"未匹配的树不得被丢弃，实得 {pairs}"
    trees = {os.path.basename(t) for t, _ in pairs}
    assert trees == {"a.tree", "b.tree"}


def test_find_dataset_pairs_one_to_one_for_equal_remainders(tmp_path):
    for name in ("x", "y"):
        (tmp_path / f"{name}.tree").write_text("((A:1,B:1)1:1,C:1)2;")
    (tmp_path / "p.tsv").write_text("1\t2\t1.0\n")
    (tmp_path / "q.tsv").write_text("1\t2\t1.0\n")
    pairs = suite.find_dataset_pairs(str(tmp_path))
    assert len(pairs) == 2
    assert {(os.path.basename(t), os.path.basename(c)) for t, c in pairs} == {
        ("x.tree", "p.tsv"),
        ("y.tree", "q.tsv"),
    }


def test_find_dataset_pairs_skips_maxtic_outputs_including_random_dist(tmp_path):
    (tmp_path / "m.tree").write_text("((A:1,B:1)1:1,C:1)2;")
    (tmp_path / "m.tsv").write_text("1\t2\t1.0\n")
    (tmp_path / "m.mt.random_dist.tsv").write_text("1.0\n2.0\n")
    (tmp_path / "m.mt.informative.tsv").write_text("1,2 1.0\n")
    (tmp_path / "m.html").write_text("<html/>")
    pairs = suite.find_dataset_pairs(str(tmp_path))
    assert pairs == [(str(tmp_path / "m.tree"), str(tmp_path / "m.tsv"))]
    assert suite._is_maxtic_output("anything.mt.random_dist.tsv") is True


# =========================================================================
# 真实基准：同样口径 + 报告列
# =========================================================================
def test_run_real_benchmark_uses_same_accounting(tmp_path):
    res = suite.run_real_benchmark(
        TREE, CONS, seed=42, n_repeats=2, warmup=1, output_prefix=str(tmp_path / "real")
    )
    assert res["error"] is None, res["error"]
    assert len(res["cpu_time_samples"]) == 2
    assert res["correctness_ok"] is True
    assert res["delivered_value"] == pytest.approx(755.88, abs=1e-6)
    assert res["n_constraints"] == 123 and res["n_distinct_pairs"] == 123
    assert 0.0 <= res["kernel_cpu_time"] <= res["cpu_time"] + 1e-12


def test_format_tables_expose_new_columns():
    res = suite.run_benchmark(n_internal=5, n_edges=8, seed=1, n_repeats=1, warmup=0)
    table = suite._format_table([res])
    for col in ("CPU(s)", "IQR", "wall(s)", "kernel(s)", "correct", "reps", "|E|real", "|E|req"):
        assert col in table, col
    assert suite._format_table([]).startswith("     n")
    real = suite.run_real_benchmark(TREE, CONS, seed=42, n_repeats=1, warmup=0)
    rtable = suite._format_real_table([real])
    for col in ("dataset", "CPU(s)", "wall(s)", "correct", "ok"):
        assert col in rtable, col


def test_format_table_flags_bad_correctness():
    row = {
        "n_internal": 3,
        "n_edges": 3,
        "n_edge_keys_effective": 3,
        "error": None,
        "exceeded": False,
        "cpu_time": 0.1,
        "cpu_time_iqr": 0.0,
        "wall_time": 0.1,
        "kernel_cpu_time": 0.09,
        "io_time": 0.01,
        "peak_memory_mb": 1.0,
        "delivered_value": 1.0,
        "n_repeats": 1,
        "correctness_ok": False,
    }
    assert "BAD" in suite._format_table([row])


def test_cli_rejects_bad_repeat_counts():
    with pytest.raises(SystemExit):
        suite.main(["--repeats", "0"])
    with pytest.raises(SystemExit):
        suite.main(["--warmup", "-1"])


def test_cli_synthetic_mode_exits_zero_and_prints_medians(capsys):
    rc = suite.main(
        ["--max-n", "6", "--max-e", "12", "--steps", "2", "--repeats", "1", "--warmup", "0"]
    )
    assert rc == 0
    out = capsys.readouterr().out
    assert "基准扫描" in out and "预热" in out and "中位数" in out
