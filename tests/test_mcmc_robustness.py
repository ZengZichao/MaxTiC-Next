"""回归测试：MCMC 的诊断口径与稳健性摘要。

覆盖范围（全部为"分析层"，报告 HTML/图表见 ``test_report_rendering.py``，
基准见 ``test_benchmark_suite.py``）：

* ``robustness/mcmc.py``：宣称克制、burn-in/thin 参数、``iters>0``/``thin>=1``
  的显式报错、ESS / 唯一样本 / 相邻重复率诊断、**独立参考 + 有判别力 + 冻结链必失败**
  的自检；并**实测断言默认温度下链只访问极少数状态**。
* ``robustness/sensitivity.py``：同时报去重序数与真实访问计数/频率、
  缺失方向记 0、``top_k`` 参数化、平凡对（方向恒定）与有信息量分歧分离。

GUI 层（``gui/node_scores.py`` 的逐节点目标值量、``gui/params.py`` 的参数域校验、
``gui/cancellation.py`` 的 StopToken / 取消回调注入 / 每任务流）已随 GUI 一起迁到
Studio 仓库的 ``tests/test_node_scores.py``、``tests/test_params_validation.py``、
``tests/test_cancellation.py``。
"""

import math

import pytest

from maxtic_next.config import MCMC_STATUS_NOTE
from maxtic_next.random_ import RandomWrapper
from maxtic_next.report import plot
from maxtic_next.robustness import mcmc as mcmc_mod
from maxtic_next.robustness import sensitivity as sens
from maxtic_next.robustness.mcmc import MCMCSampler, effective_sample_size, mcmc_diagnostics
from maxtic_next.robustness.sensitivity import NearOptimalCollector
from maxtic_next.tree.tree import Tree


def _balance_tree() -> Tree:
    tree = Tree()
    tree.read_newick(mcmc_mod._BALANCE_TEST_NEWICK)
    return tree


def _balance_edge():
    edge = dict(mcmc_mod._BALANCE_TEST_EDGE)
    return edge, sorted(edge.keys())


def _sampler(temperature=0.01, seed=7):
    edge, keys = _balance_edge()
    return MCMCSampler(_balance_tree(), edge, keys, RandomWrapper(seed), temperature=temperature)


# =========================================================================
# 措辞：不再有"严格 / 后验 / 独立样本"宣称
# =========================================================================
def test_mcmc_wording_is_honest():
    src = open(mcmc_mod.__file__, encoding="utf-8").read()
    # 模块 docstring 必须使用交付的官方口径，且不再自称"严格采样器"
    assert "Metropolis" in src and "linear extensions" in src
    assert "preliminary" in src and "not validated" in src
    assert "严格采样器" not in src
    # "独立样本" 仅应在解释"为何不成立"的段落中出现，不应在官方口径声明之前出现
    assert "独立样本" not in src.split("为什么不再写")[0]
    assert MCMC_STATUS_NOTE == MCMCSampler.STATUS_NOTE
    assert "posterior" not in MCMC_STATUS_NOTE.lower()
    assert "rigorous" not in MCMC_STATUS_NOTE.lower()


# =========================================================================
# 默认温度下链几乎是冻结的：诊断必须把这件事说出来
# =========================================================================
def test_default_temperature_reports_frozen_chain():
    sampler = _sampler(temperature=0.01)
    sampler.sample(2000, thin=1)
    diag = sampler.diagnostics()
    n_states = diag["state_space_size"]
    assert n_states == 80, "合成偏序实例应恰有 80 个线性扩展"
    # ：默认档几千步只访问极少数状态（示例数据 4000 步仅 3/80）
    assert diag["n_samples"] == 2000
    assert diag["n_unique_samples"] <= 20
    assert diag["state_coverage_fraction"] <= 0.25
    assert diag["unique_sample_fraction"] < 0.02
    assert diag["adjacent_duplicate_fraction"] > 0.7
    assert diag["independent_samples"] is False
    # ESS 必须远小于样本数（"看着 2000 个样本、其实只有一条解"）
    assert diag["effective_sample_size"] <= 20.0
    # 均值/最小值这类统计量在冻结链上没有统计意义，必须与诊断一并阅读
    assert diag["energy_min"] <= diag["energy_mean"]


def test_higher_temperature_mixes():
    cold = _sampler(temperature=0.01)
    cold.sample(4000, thin=1)
    hot = _sampler(temperature=2.0)
    hot.sample(4000, thin=1)
    cd, hd = cold.diagnostics(), hot.diagnostics()
    assert hd["state_coverage_fraction"] > 0.5
    assert hd["effective_sample_size"] > 50
    assert hd["n_unique_samples"] > 5 * cd["n_unique_samples"]
    assert hd["adjacent_duplicate_fraction"] < cd["adjacent_duplicate_fraction"]


# =========================================================================
# 参数校验：no ZeroDivisionError / no min([])
# =========================================================================
@pytest.mark.parametrize("bad", [0, -1, -100])
def test_sample_rejects_nonpositive_iters(bad):
    sampler = _sampler()
    with pytest.raises(ValueError) as exc:
        sampler.sample(bad)
    assert "iters" in str(exc.value)


@pytest.mark.parametrize("bad", [0, -3])
def test_sample_rejects_bad_thin(bad):
    sampler = _sampler()
    with pytest.raises(ValueError):
        sampler.sample(10, thin=bad)


def test_constructor_rejects_bad_burn_in_and_thin():
    edge, keys = _balance_edge()
    with pytest.raises(ValueError):
        MCMCSampler(_balance_tree(), edge, keys, RandomWrapper(1), temperature=1.0, thin=0)
    with pytest.raises(ValueError):
        MCMCSampler(_balance_tree(), edge, keys, RandomWrapper(1), temperature=1.0, burn_in=-1)
    with pytest.raises(ValueError):
        MCMCSampler(_balance_tree(), edge, keys, RandomWrapper(1), temperature=0.0)


def test_thinning_and_burn_in_shapes():
    sampler = _sampler(temperature=1.0)
    out = sampler.sample(100, thin=7, burn_in=50)
    assert len(out) == math.ceil(100 / 7)
    assert sampler.last_burn_in == 50
    for order in out:
        assert sampler.is_valid(order)
    assert len(sampler.last_energies) == len(out)
    assert len(sampler.last_chain_energies) == 100


def test_burn_in_params_are_stored_on_constructor():
    edge, keys = _balance_edge()
    s = MCMCSampler(
        _balance_tree(), edge, keys, RandomWrapper(3), temperature=1.0, burn_in=25, thin=5
    )
    assert s.burn_in == 25 and s.thin == 5
    out = s.sample(60)  # 用构造期默认 burn-in/thin
    assert s.last_burn_in == 25
    assert len(out) == 12


def test_diagnostics_of_empty_chain_does_not_crash():
    d = mcmc_diagnostics([], [])
    assert d["n_samples"] == 0 and d["energy_min"] is None
    assert d["effective_sample_size"] == 0.0


# =========================================================================
# ESS / 相邻重复率 / 唯一样本诊断
# =========================================================================
def test_effective_sample_size_helpers():
    n = 400
    iid = [float(i % 7) for i in range(n)]
    frozen = [1.0] * n
    assert effective_sample_size(frozen) == 1.0
    assert effective_sample_size(iid) > 50
    d = mcmc_diagnostics([[str(i % 3)] for i in range(n)], [float(i % 3) for i in range(n)])
    assert d["n_unique_samples"] == 3
    assert d["adjacent_duplicate_fraction"] == 0.0
    assert d["independent_samples"] is False
    # 有连续重复的链：相邻重复率被量化出来
    repeated = [[str(i // 4 % 3)] for i in range(n)]
    d2 = mcmc_diagnostics(repeated, [float(i // 4 % 3) for i in range(n)])
    assert 0.6 < d2["adjacent_duplicate_fraction"] < 0.9
    assert d2["max_chain_run_fraction"] >= 4 / n
    # 冻结链（全部相同）：ESS = 1
    frozen_diag = mcmc_diagnostics([["a", "b"]] * 50, [1.0] * 50)
    assert frozen_diag["effective_sample_size"] == 1.0
    assert frozen_diag["adjacent_duplicate_fraction"] == 1.0
    assert frozen_diag["n_unique_samples"] == 1


# =========================================================================
# 自检：独立参考 + 有判别力 + 冻结链不得蒙混
# =========================================================================
class _EnergyBlind(MCMCSampler):
    """坏链：只要提议合法就接受（完全忽略能量）。"""

    def _advance(self):
        cand, valid = self._propose_swap(self.current, self.rng)
        if valid:
            self.current = cand
            self.current_energy = self.value_computer.value(cand)
        return self.current, self.current_energy


class _Frozen(MCMCSampler):
    """坏链：从不移动（模拟默认低温下实际发生的事）。"""

    def _advance(self):
        self.rng.random()
        return self.current, self.current_energy


def _report_of(cls, **kw):
    edge, keys = _balance_edge()
    s = cls(_balance_tree(), edge, keys, RandomWrapper(11), temperature=0.01)
    return s.self_test_detailed_balance_report(**kw)


def test_self_test_passes_with_independent_reference():
    rep = _report_of(MCMCSampler, n_iter=20000, tolerance=0.05)
    assert rep["passed"] is True
    assert rep["necessary_not_sufficient"] is True
    assert "必要非充分" in rep["interpretation"]
    # 参考值不是用被检代码算的：独立枚举 + 位掩码 DP 交叉校验
    assert rep["dp_count_matches_enumeration"] is True
    assert rep["enumeration_agrees_with_is_valid"] is True
    assert rep["state_space_size"] == 80
    assert rep["state_coverage_fraction"] == 1.0
    assert rep["power_vs_energy_blind"] > 0.05


def test_self_test_rejects_energy_blind_chain():
    rep = _report_of(_EnergyBlind, n_iter=20000, tolerance=0.05)
    assert rep["passed"] is False
    assert any("独立参考" in r or "能量盲" in r or "> 允许" in r for r in rep["failure_reasons"])


def test_self_test_rejects_frozen_chain():
    rep = _report_of(_Frozen, n_iter=20000, tolerance=0.05)
    assert rep["passed"] is False
    assert rep["state_coverage_fraction"] < 0.1
    assert any("冻结" in r for r in rep["failure_reasons"])


def test_self_test_bool_wrapper_and_alias():
    s = _sampler(temperature=1.0)
    assert s.self_test_detailed_balance(n_iter=5000, tolerance=0.05) in (True, False)
    rep = s.self_test_report(n_iter=5000, tolerance=0.05)
    assert isinstance(rep, dict) and "passed" in rep


def test_reference_energy_helper_is_independent_of_value_computer():
    edge, keys = _balance_edge()
    order = ["R", "Z", "X", "W", "U", "Y", "V"]
    from maxtic_next.ranking.value import value as ref_value

    assert mcmc_mod._reference_energy(order, edge) == ref_value(order, edge, keys)


# =========================================================================
# ：真实访问计数 / 频率，缺失方向记 0；top_k 参数化
# =========================================================================
def test_summary_reports_true_access_counts_not_dedup_only():
    col = NearOptimalCollector(top_k=10)
    for _ in range(100):
        col.add(["x", "y", "z"], 1.0)
    col.add(["y", "x", "z"], 2.0)
    s = col.summary()
    assert s["n_orders_collected"] == 2  # 去重序数
    assert s["n_accesses_retained"] == 101  # 真实访问次数
    assert s["n_accesses_total"] == 101
    assert s["n_unique_orders_total"] == 2
    f = s["pairwise_order_frequency"]
    assert f[("x", "y")] == 100 and f[("y", "x")] == 1  # 反向不再省略，而是记真实值
    assert f[("z", "x")] == 0  # 从未出现的方向记 0
    p = s["pairwise_order_proportion"]
    assert abs(p[("x", "y")] - 100 / 101) < 1e-12
    assert abs(p[("x", "y")] + p[("y", "x")] - 1.0) < 1e-12
    # 节点位置分布补齐所有位置（未访问的位置为 0）
    assert s["node_position_distribution"]["x"] == {0: 100, 1: 1, 2: 0}
    assert s["is_resampling_robustness"] is False
    assert "重抽样" in s["basis"] and "访问" in s["basis"]


def test_summary_top_k_is_a_parameter():
    col = NearOptimalCollector(top_k=5)
    col.add(["a", "b", "c"], 1.0)
    col.add(["b", "a", "c"], 2.0)
    col.add(["c", "a", "b"], 3.0)
    full = col.summary()
    cut = col.summary(top_k=1)
    assert full["n_orders_collected"] == 3
    assert cut["n_orders_collected"] == 1
    assert cut["top_k_reported"] == 1
    assert cut["n_accesses_retained"] == 1
    assert full["top_k"] == 5
    with pytest.raises(ValueError):
        col.summary(top_k=0)


def test_default_top_k_constant_is_exposed():
    assert sens.DEFAULT_TOP_K == 50
    assert NearOptimalCollector().top_k == sens.DEFAULT_TOP_K
    assert NearOptimalCollector(top_k=500).top_k == 500


# =========================================================================
# ：方向恒定的平凡对不再淹没"有信息量"的成对表
# =========================================================================
def _toy_summary():
    col = NearOptimalCollector(top_k=10)
    # 6 个近优序：root 恒在最前（拓扑强制的平凡对），a/b 互换（真正的分歧）
    for _ in range(30):
        col.add(["root", "a", "b", "leaf"], 1.0)
    for _ in range(20):
        col.add(["root", "b", "a", "leaf"], 1.0)
    return col.summary()


def test_trivial_pairs_are_separated_from_informative_ones():
    s = _toy_summary()
    trivial = {tuple(p) for p in s["constant_direction_pairs"]}
    assert ("root", "leaf") in trivial and ("root", "a") in trivial
    info = s["informative_pairwise_order_frequency"]
    assert ("a", "b") in info and ("b", "a") in info
    assert ("root", "leaf") not in info and ("root", "a") not in info


def test_plot_pair_table_does_not_flood_with_trivial_pairs():
    html = plot.robustness_html(_toy_summary(), max_pairs=2)
    # top-2 必须是有信息量的分歧对（a/b 两个方向），而不是恒定的 root≺* 平凡对
    body = html.split("成对次序频率")[1]
    assert "&gt;" not in body.split("<tbody>")[1].split("</tr>")[0]
    first_rows = body.split("<tbody>")[1].split("</tbody>")[0]
    assert first_rows.count("<tr>") <= 2
    for row in first_rows.split("</tr>"):
        assert ">root<" not in row  # 平凡对被过滤
    assert "a" in first_rows and "b" in first_rows
    assert "平凡对" in html  # 有可见说明 + 折叠清单


def test_plot_shows_both_dedup_count_and_access_counts():
    html = plot.robustness_html(_toy_summary())
    assert "去重近优排序数" in html
    assert "真实访问次数" in html
    assert ">50<" in html  # 30 + 20 次访问


def test_esc_handles_integer_zero():
    assert plot._esc(0) == "0"
    assert plot._esc("") == ""
    assert plot._esc(None) == ""
    assert plot._esc("a<b") == "a&lt;b"
    html = plot.robustness_html(
        {
            "title": "t",
            "n_orders_collected": 1,
            "node_position_distribution": {"root": {0: 5, 1: 0}},
            "pairwise_order_frequency": {},
        }
    )
    assert "<td>0</td>" in html  # 众数位置 0 渲染成 "0"，不是空白
