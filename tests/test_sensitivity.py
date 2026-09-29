"""近优解收集与稳健性/敏感性摘要单元测试。

术语铁律：本测试显式断言摘要中**不出现**"置信区间 / 后验概率 /
后验分布"等用语——这些是 P2-MCMC 的专属词汇，局部搜索访问分布不是后验分布。
"""

from maxtic_next.robustness.sensitivity import NearOptimalCollector


def test_collector_dedup_and_capacity_trim():
    """add 自动去重（计数+1）并按 top_k 裁剪最差解。"""
    col = NearOptimalCollector(top_k=3)
    col.add(["a", "b", "c"], 5.0)
    col.add(["a", "b", "c"], 5.0)  # 重复 -> 访问计数 +1
    col.add(["c", "b", "a"], 3.0)
    col.add(["a", "c", "b"], 7.0)
    col.add(["b", "a", "c"], 9.0)  # 第 4 个，最差(9.0)应被裁剪
    assert len(col) == 3
    assert ("b", "a", "c") not in col._visited
    # 去重：("a","b","c") 累计访问 2 次
    assert col._visited[("a", "b", "c")][1] == 2


def test_summary_title_and_no_posterior_language():
    """摘要标题必须为"基于局部搜索访问解的稳健性/敏感性摘要"，且不含后验用语。"""
    col = NearOptimalCollector(top_k=10)
    col.add(["a", "b", "c"], 1.0)
    col.add(["a", "c", "b"], 2.0)
    s = col.summary()
    assert s["title"] == "基于局部搜索访问解的稳健性/敏感性摘要"
    text = str(s)
    assert "置信区间" not in text
    assert "后验概率" not in text
    assert "后验分布" not in text


def test_summary_statistics():
    """节点位置分布与成对次序频率统计正确。"""
    col = NearOptimalCollector(top_k=10)
    col.add(["a", "b", "c"], 1.0)
    col.add(["a", "c", "b"], 2.0)
    s = col.summary()
    # a 在两个序中位置都是 0
    assert s["node_position_distribution"]["a"][0] == 2
    # a 在 b 之前：两个序都满足 -> 频率 2
    assert s["pairwise_order_frequency"][("a", "b")] == 2
    # best_value 取最小值
    assert s["best_value"] == 1.0


def test_top_k_must_be_positive():
    """top_k < 1 必须抛 ValueError。"""
    try:
        NearOptimalCollector(top_k=0)
        assert False, "应当抛出 ValueError"
    except ValueError:
        pass
