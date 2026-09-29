"""真实数据集基准测试（ / P1-PERF）。

验证 ``benchmarks/suite.py`` 中新增的真实数据基准功能：

* :func:`run_real_benchmark`：单组真实数据运行
* :func:`run_real_suite`：目录批量运行
* :func:`find_dataset_pairs`：数据对自动发现
* CLI ``--real`` / ``--data-dir`` 标志
"""

import os

import pytest

from maxtic_next.benchmarks.suite import (
    find_dataset_pairs,
    run_real_benchmark,
    run_real_suite,
    _format_real_table,
    _count_constraints,
    _count_internal_nodes,
    main,
)

# 测试数据路径（tests/data/ 中的 minitree + Cyano 真实数据）
_TESTS_DIR = os.path.dirname(os.path.abspath(__file__))
TEST_DATA_DIR = os.path.join(_TESTS_DIR, "data")
TEST_TREE = os.path.join(TEST_DATA_DIR, "minitree.tree")
TEST_CONS = os.path.join(TEST_DATA_DIR, "Cyano_CUTConstraints.tsv")

# examples/ 目录（项目根目录下）
_PROJECT_ROOT = os.path.dirname(_TESTS_DIR)
EXAMPLES_DIR = os.path.join(_PROJECT_ROOT, "examples")


# ----------------------------------------------------------------------
# 辅助函数测试
# ----------------------------------------------------------------------
class TestHelperFunctions:
    """测试真实基准的辅助函数。"""

    def test_count_constraints(self):
        """约束行数统计正确。"""
        n = _count_constraints(TEST_CONS)
        assert n > 0, f"约束文件 {TEST_CONS} 行数应为正"
        # 手动统计验证
        with open(TEST_CONS) as fh:
            manual = sum(1 for line in fh if line.strip())
        assert n == manual

    def test_count_internal_nodes(self):
        """内部节点数统计正确。"""
        n = _count_internal_nodes(TEST_TREE)
        assert n > 0, f"树 {TEST_TREE} 内部节点数应为正"
        # minitree 有 13 个内部节点（42,41,59,43,45,61,56,65,40,46,62,67,69）
        assert n == 13, f"minitree 应有 13 个内部节点，实际 {n}"

    def test_count_constraints_nonexistent(self):
        """不存在的文件返回 0。"""
        assert _count_constraints("/nonexistent/path.tsv") == 0

    def test_count_internal_nodes_nonexistent(self):
        """不存在的树文件返回 0。"""
        assert _count_internal_nodes("/nonexistent/path.tree") == 0


# ----------------------------------------------------------------------
# find_dataset_pairs 测试
# ----------------------------------------------------------------------
class TestFindDatasetPairs:
    """测试数据对自动发现功能。"""

    def test_find_pairs_real_data_dir(self):
        """在 tests/data/ 目录中能找到 tree + constraints 对。"""
        pairs = find_dataset_pairs(TEST_DATA_DIR)
        assert len(pairs) >= 1, f"应至少找到 1 对数据，实际 {len(pairs)}"
        # 检查找到的对包含 minitree.tree
        tree_names = [os.path.basename(p[0]) for p in pairs]
        assert "minitree.tree" in tree_names, f"应找到 minitree.tree，实际找到 {tree_names}"

    def test_find_pairs_examples_dir(self):
        """在 examples/ 目录中能找到 tree + constraints 对。"""
        if not os.path.isdir(EXAMPLES_DIR):
            pytest.skip(f"examples/ 目录不存在: {EXAMPLES_DIR}")
        pairs = find_dataset_pairs(EXAMPLES_DIR)
        assert len(pairs) >= 1, f"应至少找到 1 对数据，实际 {len(pairs)}"

    def test_find_pairs_nonexistent_dir(self):
        """不存在的目录返回空列表。"""
        pairs = find_dataset_pairs("/nonexistent/directory/")
        assert pairs == []

    def test_find_pairs_name_matching(self, tmp_path):
        """同前缀的文件应按名称配对。"""
        # 创建同前缀的树和约束文件
        (tmp_path / "dataset1.tree").write_text("((A:1,B:1)1:1,C:1)2;")
        (tmp_path / "dataset1.tsv").write_text("1\t2\t1.0\n")
        (tmp_path / "dataset2.tree").write_text("((A:1,B:1)1:1,C:1)2;")
        (tmp_path / "dataset2.tsv").write_text("1\t2\t1.0\n")

        pairs = find_dataset_pairs(str(tmp_path))
        assert len(pairs) == 2
        stems = {os.path.splitext(os.path.basename(p[0]))[0] for p in pairs}
        assert stems == {"dataset1", "dataset2"}

    def test_find_pairs_single_pair_no_name_match(self, tmp_path):
        """目录内恰好 1 树 + 1 约束但名称不同，应自动配对。"""
        (tmp_path / "species.tree").write_text("((A:1,B:1)1:1,C:1)2;")
        (tmp_path / "constraints.tsv").write_text("1\t2\t1.0\n")

        pairs = find_dataset_pairs(str(tmp_path))
        assert len(pairs) == 1


# ----------------------------------------------------------------------
# run_real_benchmark 测试
# ----------------------------------------------------------------------
class TestRunRealBenchmark:
    """测试单组真实数据基准。"""

    def test_run_real_benchmark_success(self, tmp_path):
        """真实数据基准可运行、无异常、指标字段合法。"""
        out = str(tmp_path / "real_bench")
        res = run_real_benchmark(TEST_TREE, TEST_CONS, seed=42, output_prefix=out, time_budget=60.0)

        # 基本字段
        assert res["tree_path"] == TEST_TREE
        assert res["constraints_path"] == TEST_CONS
        assert res["n_internal"] == 13  # minitree 有 13 个内部节点
        assert res["n_constraints"] > 0
        assert res["error"] is None, f"基准出错: {res['error']}"

        # CPU 时间与内存
        assert res["cpu_time"] is not None and res["cpu_time"] >= 0.0
        assert res["peak_memory_mb"] is not None and res["peak_memory_mb"] >= 0.0

        # 排序结果
        assert res["best_value"] is not None and res["best_value"] >= 0.0
        assert res["best_source"] in ("greedy heuristic", "mixing heuristic")
        assert 0.0 <= res["similarity_to_input"] <= 1.0
        assert res["exceeded"] is False  # 小数据集不应超时

    def test_run_real_benchmark_default_output_prefix(self, tmp_path):
        """不指定 output_prefix 时使用临时目录，不报错。"""
        res = run_real_benchmark(TEST_TREE, TEST_CONS, seed=42, time_budget=60.0)
        assert res["error"] is None

    def test_run_real_benchmark_different_seeds(self, tmp_path):
        """不同种子均可运行且产出合法结果。"""
        for seed in [1, 42, 99]:
            out = str(tmp_path / f"real_seed{seed}")
            res = run_real_benchmark(
                TEST_TREE, TEST_CONS, seed=seed, output_prefix=out, time_budget=60.0
            )
            assert res["error"] is None, f"seed={seed} 出错: {res['error']}"
            assert res["best_value"] is not None

    def test_run_real_benchmark_deterministic(self, tmp_path):
        """同种子运行两次，best_value 与 best_source 完全一致。"""
        r1 = run_real_benchmark(
            TEST_TREE, TEST_CONS, seed=42, output_prefix=str(tmp_path / "d1"), time_budget=60.0
        )
        r2 = run_real_benchmark(
            TEST_TREE, TEST_CONS, seed=42, output_prefix=str(tmp_path / "d2"), time_budget=60.0
        )
        assert r1["best_value"] == r2["best_value"]
        assert r1["best_source"] == r2["best_source"]
        assert r1["similarity_to_input"] == r2["similarity_to_input"]

    def test_run_real_benchmark_nonexistent_files(self, tmp_path):
        """不存在的文件应捕获异常而非崩溃。"""
        res = run_real_benchmark("/nonexistent.tree", "/nonexistent.tsv", seed=42, time_budget=60.0)
        assert res["error"] is not None
        assert res["cpu_time"] is None
        assert res["best_value"] is None


# ----------------------------------------------------------------------
# run_real_suite 测试
# ----------------------------------------------------------------------
class TestRunRealSuite:
    """测试目录批量基准。"""

    def test_run_real_suite_tests_data(self, tmp_path, capsys):
        """在 tests/data/ 目录上运行批量基准。"""
        results = run_real_suite(TEST_DATA_DIR, seed=42, time_budget=60.0)
        assert len(results) >= 1
        for r in results:
            assert r["error"] is None, f"基准出错: {r['error']}"
            assert r["cpu_time"] is not None

    def test_run_real_suite_empty_dir(self, tmp_path, capsys):
        """空目录返回空列表并打印警告。"""
        results = run_real_suite(str(tmp_path), seed=42, time_budget=60.0)
        assert results == []

    def test_run_real_suite_nonexistent_dir(self, tmp_path, capsys):
        """不存在的目录返回空列表。"""
        results = run_real_suite("/nonexistent/dir/", seed=42, time_budget=60.0)
        assert results == []

    def test_run_real_suite_custom_dir(self, tmp_path, capsys):
        """自定义目录中放入数据对，批量基准可运行。"""
        # 创建测试数据
        tree_content = "((A:1,B:1)1:1,C:1)2;"
        cons_content = "1\t2\t1.0\n2\t1\t0.5\n"
        (tmp_path / "custom.tree").write_text(tree_content)
        (tmp_path / "custom.tsv").write_text(cons_content)

        results = run_real_suite(str(tmp_path), seed=42, time_budget=60.0)
        assert len(results) == 1
        assert results[0]["error"] is None
        assert results[0]["n_internal"] == 2
        assert results[0]["n_constraints"] == 2


# ----------------------------------------------------------------------
# 格式化与 CLI 测试
# ----------------------------------------------------------------------
class TestFormatAndCLI:
    """测试结果格式化与 CLI 入口。"""

    def test_format_real_table_with_results(self, tmp_path):
        """格式化函数能正确处理有结果的数据。"""
        out = str(tmp_path / "fmt")
        res = run_real_benchmark(TEST_TREE, TEST_CONS, seed=42, output_prefix=out, time_budget=60.0)
        table = _format_real_table([res])
        assert "minitree" in table
        assert "CPU(s)" in table
        assert "ok" in table

    def test_format_real_table_with_error(self):
        """格式化函数能正确处理错误结果。"""
        error_res = {
            "tree_path": "/bad.tree",
            "constraints_path": "/bad.tsv",
            "n_internal": 0,
            "n_constraints": 0,
            "cpu_time": None,
            "peak_memory_mb": None,
            "best_value": None,
            "best_source": None,
            "similarity_to_input": None,
            "exceeded": False,
            "error": "FileNotFoundError: bad",
        }
        table = _format_real_table([error_res])
        assert "ERROR" in table

    def test_format_real_table_empty(self):
        """空结果列表格式化不报错。"""
        table = _format_real_table([])
        assert "dataset" in table  # header 仍存在

    def test_cli_real_flag(self, capsys):
        """CLI --real 标志使用默认 examples/ 目录运行。"""
        if not os.path.isdir(EXAMPLES_DIR):
            pytest.skip("examples/ 目录不存在")
        ret = main(["--real", "--seed", "42", "--time-budget", "60"])
        assert ret == 0
        captured = capsys.readouterr().out
        assert "真实数据集基准" in captured or "运行" in captured

    def test_cli_real_with_data_dir(self, capsys):
        """CLI --real --data-dir 指定 tests/data/ 目录。"""
        ret = main(["--real", "--data-dir", TEST_DATA_DIR, "--seed", "42", "--time-budget", "60"])
        assert ret == 0
        captured = capsys.readouterr().out
        assert "真实数据集基准" in captured

    def test_cli_real_empty_data_dir(self, tmp_path, capsys):
        """CLI --real --data-dir 指向空目录时返回 0 且打印警告。"""
        empty_dir = str(tmp_path / "empty")
        os.makedirs(empty_dir, exist_ok=True)
        ret = main(["--real", "--data-dir", empty_dir, "--seed", "42", "--time-budget", "60"])
        assert ret == 0
        captured = capsys.readouterr().out
        assert "未找到" in captured or "警告" in captured

    def test_cli_default_synthetic_mode(self, capsys):
        """不带 --real 时仍走合成数据模式（回归保护）。"""
        ret = main(
            ["--max-n", "5", "--max-e", "8", "--steps", "1", "--seed", "1", "--time-budget", "30"]
        )
        assert ret == 0
        captured = capsys.readouterr().out
        assert "基准扫描" in captured
