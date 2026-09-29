"""性能基准单元测试。

验证合成数据生成与基准运行正确性；基准内部捕获一切异常，故断言聚焦于
"无异常返回 + 指标字段合法"。
"""

from maxtic_next.benchmarks.suite import (
    make_synthetic_constraints,
    make_synthetic_tree,
    run_benchmark,
    run_suite,
)
from maxtic_next.tree.tree import Tree


def test_make_synthetic_tree_shape():
    """合成树含恰好 n 个内部节点且可解析。"""
    newick, labels = make_synthetic_tree(7, seed=1)
    tree = Tree()
    tree.read_newick(newick)  # 解析失败会抛异常，本身即校验
    assert len(labels) == 7
    assert len(tree.internal_node_labels()) == 7


def test_make_synthetic_constraints_count():
    """合成约束数恰好为 n_edges。"""
    _, labels = make_synthetic_tree(5, seed=1)
    cset = make_synthetic_constraints(labels, 12, seed=1)
    assert len(cset.constraints) == 12


def test_run_benchmark_small():
    """小规模基准可运行、无异常、指标字段合法。"""
    res = run_benchmark(n_internal=7, n_edges=10, seed=1, time_budget=30.0)
    assert res["n_internal"] == 7
    assert res["n_edges"] == 10
    assert res["error"] is None
    assert res["cpu_time"] is not None and res["cpu_time"] >= 0.0
    assert res["peak_memory_mb"] is not None and res["peak_memory_mb"] >= 0.0


def test_run_suite_small():
    """基准网格返回 (len(n_list) * len(e_list)) 个结果且均成功。"""
    results = run_suite([3, 5], [4, 8], seed=2, time_budget=30.0)
    assert len(results) == 4
    for r in results:
        assert r["error"] is None
