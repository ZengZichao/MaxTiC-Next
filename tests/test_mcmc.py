# Copyright (C) 2026 MaxTiC-Next rewrite team.
# This file is part of MaxTiC-Next, a Python 3 rewrite of MaxTiC.
#
# MaxTiC is Copyright (C) Eric Tannier (Inria / CNRS / ENS Lyon).
# Original reference: MaxTiC: Fast ranking of a phylogenetic tree by
# Maximum Time Consistency with lateral gene transfers,
# Biorxiv doi.org/10.1101/127548
#
# MaxTiC-Next is free software: you can redistribute it and/or modify it
# under the terms of the CeCILL 2.1 license as published by the
# Commissariat a l'energie atomique et aux energies alternatives (CEA),
# the Centre National de la Recherche Scientifique (CNRS) and the
# Institut National de Recherche en Informatique et en Automatique (INRIA).
# A copy of the license is included in the LICENSE file.

"""P2-MCMC 严格采样器回归测试。

四类验证（对应 GAP 2 验收标准）：

(a) 合法线性扩展：``sample()`` 产出的每个排序都是偏序的合法线性扩展
    （满足全部树边约束、且为内部节点的一个排列）；
(b) 细致平衡自检：``self_test_detailed_balance()`` 在合成实例上经验验证
    采样器收敛到 Boltzmann 后验（经验频率 ≈ 理论平稳分布）；
(c) 确定性：相同种子 + 相同初始排序 -> 两次 ``sample()`` 结果逐元素一致；
(d) 低温聚焦最优：温度足够低时，采样能量最小值等于真实最优（=0），
    且均值逼近最优；高温时均值介于两极值之间（≈ 0.5）。

术语铁律：本测试仅校验 ``robustness/mcmc.py`` 的采样行为，不直接产出
"置信区间 / 后验分布" 报告用语（报告用语由该模块统一承载）。
"""

from typing import List

from maxtic_next.random_ import RandomWrapper
from maxtic_next.robustness.mcmc import MCMCSampler
from maxtic_next.tree.tree import Tree

# 受控合成树：((A,B)X,(C,D)Y)Z —— 内部节点 X,Y,Z；偏序 Z<X, Z<Y（X,Y 不可比）
# => 恰有两个线性扩展：[Z,X,Y] 与 [Z,Y,X]。
TREE_NEWICK = "((A:1.0,B:1.0)X:1.0,(C:1.0,D:1.0)Y:1.0)Z:1.0;"
# 约束边 "X,Y"（权重 1）：要求 X 早于 Y。
#   * [Z,X,Y]：X(1) < Y(2) -> 无违反，能量 0（= 最优）；
#   * [Z,Y,X]：X(2) > Y(1) -> 违反，能量 1。
EDGE = {"X,Y": 1.0}
EDGE_KEYS = ["X,Y"]


def _build(temperature: float, seed: int, initial_order: List[str] = None) -> MCMCSampler:
    """构造受控树上的采样器（统一封装 Newick 解析 + 随机驱动）。"""
    tree = Tree()
    tree.read_newick(TREE_NEWICK)
    rng = RandomWrapper(seed)
    return MCMCSampler(
        tree=tree,
        edge=EDGE,
        edge_keys=EDGE_KEYS,
        rng=rng,
        temperature=temperature,
        initial_order=initial_order,
    )


# ---------------------------------------------------------------------------
# (a) 合法线性扩展
# ---------------------------------------------------------------------------
def test_samples_are_valid_linear_extensions():
    """sample() 的每个输出都应是合法线性扩展（满足偏序且为全排列）。"""
    sampler = _build(temperature=0.01, seed=1)
    samples = sampler.sample(500)
    assert len(samples) == 500
    for order in samples:
        # 必须是内部节点集合的一个排列
        assert set(order) == set(sampler.internal_nodes), f"非排列: {order}"
        # 必须尊重全部树边偏序
        assert sampler.is_valid(order), f"非合法线性扩展: {order}"


def test_sample_of_size_one_still_valid():
    """即使只采 1 步，仍返回合法线性扩展。"""
    sampler = _build(temperature=0.01, seed=99)
    samples = sampler.sample(1)
    assert len(samples) == 1
    assert sampler.is_valid(samples[0])


# ---------------------------------------------------------------------------
# (b) 细致平衡自检
# ---------------------------------------------------------------------------
def test_self_test_detailed_balance_passes():
    """self_test_detailed_balance 在小规模合成实例上应返回 True。"""
    sampler = _build(temperature=0.01, seed=7)
    ok = sampler.self_test_detailed_balance(n_iter=20000, tolerance=0.05)
    assert ok is True


def test_self_test_detailed_balance_works_at_higher_temperature():
    """较高温度下自检仍应通过（平稳分布定义与温度无关，仅形状变化）。"""
    sampler = _build(temperature=1.0, seed=7)
    ok = sampler.self_test_detailed_balance(n_iter=40000, tolerance=0.05)
    assert ok is True


# ---------------------------------------------------------------------------
# (c) 确定性：固定种子 -> 逐元素一致
# ---------------------------------------------------------------------------
def test_determinism_with_fixed_seed():
    """相同种子 + 相同初始排序 -> 两次 sample() 结果完全一致。"""
    init = ["Z", "X", "Y"]
    s1 = _build(temperature=0.01, seed=1, initial_order=init).sample(300)
    s2 = _build(temperature=0.01, seed=1, initial_order=init).sample(300)
    assert s1 == s2


def test_determinism_with_default_initial_order():
    """初始排序取默认拓扑序（确定性）时，相同种子也应逐元素一致。"""
    s1 = _build(temperature=0.05, seed=123).sample(200)
    s2 = _build(temperature=0.05, seed=123).sample(200)
    assert s1 == s2


# ---------------------------------------------------------------------------
# (d) 低温聚焦最优 / 高温分散
# ---------------------------------------------------------------------------
def test_low_temperature_concentrates_on_optimum():
    """低温（T=0.001）时：最小值 == 最优(0)，且均值逼近最优。"""
    sampler = _build(temperature=0.001, seed=3)
    samples = sampler.sample(1000)
    energies = [sampler.value_computer.value(s) for s in samples]
    # 达到真实最优（能量 0）
    assert min(energies) == 0.0
    # 均值逼近最优（低温后验几乎全部质量落在最优态）
    mean_e = sum(energies) / len(energies)
    assert abs(mean_e - 0.0) < 0.05


def test_high_temperature_spreads_between_extremes():
    """高温（T=10）时：均值应介于两极值 (0 与 1) 之间（≈ 0.5）。"""
    sampler = _build(temperature=10.0, seed=3)
    samples = sampler.sample(1000)
    energies = [sampler.value_computer.value(s) for s in samples]
    mean_e = sum(energies) / len(energies)
    # 两线性扩展能量分别为 0 与 1，高温下平稳分布近似均匀
    assert 0.35 < mean_e < 0.65


def test_min_energy_equals_optimum_at_low_temperature():
    """独立验证：低温下采样集合的最小值等于该受控实例的理论最优（0）。"""
    # 理论最优：可追溯地计算两个线性扩展的能量，取最小
    optimum = min(
        MCMCSampler(
            tree=(lambda t: t.read_newick(TREE_NEWICK) or t)(Tree()),
            edge=EDGE,
            edge_keys=EDGE_KEYS,
            rng=RandomWrapper(0),
        ).value_computer.value(order)
        for order in (["Z", "X", "Y"], ["Z", "Y", "X"])
    )
    assert optimum == 0.0

    sampler = _build(temperature=0.0001, seed=5)
    samples = sampler.sample(500)
    energies = [sampler.value_computer.value(s) for s in samples]
    assert min(energies) == optimum
