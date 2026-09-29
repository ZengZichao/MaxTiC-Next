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

"""线性扩展上的 Metropolis–Hastings 采样——**初步实现**。

官方状态说明（与 ``config.MCMC_STATUS_NOTE`` 逐字一致，也是 stdout / HTML 报告
必须展示的口径）::

    Metropolis–Hastings over linear extensions of the species-tree partial
    order (preliminary; convergence diagnostics not validated)

------------------------------------------------------------------------------
：为什么不再写"严格 / 后验采样 / 独立样本"
------------------------------------------------------------------------------
接受率本身是标准 MH（对称提议 ⇒ MH 比率精确满足细致平衡方程），但**参数域与
"后验采样器"的宣称不符**，因此本模块的措辞降级为如实描述：

* 能量是约束权重之和（示例数据 ``total_weight = 2218.4``）。在旧的默认温度
  ``T = 0.01`` 下 uphill 接受率约 ``e^-50``…``e^-100``，实测 4000 步只访问
  80 个线性扩展中的 **3** 个（其中 2 个从未被访问）——链在实践时间尺度上
  **不可约**，其输出不是任何目标分布的样本，且强烈依赖 ``initial_order``。
* 马尔可夫链样本**自相关**，不能当独立样本处理：``thin=1`` 时实测相邻样本
  73% 完全相同。故本模块强制提供 burn-in / thinning，并把
  "唯一样本数 / 相邻重复率 / 有效样本量 (ESS) / 状态访问比例"作为
  **必须展示**的诊断量（见 :meth:`MCMCSampler.diagnostics`）。
* 内置自检 :meth:`MCMCSampler.self_test_report` 只给出**必要非充分**证据
  （见其文档），并且若链"冻结"（状态访问比例过低）会直接判为**不通过**。

因此：本模块产出的是"一条在偏序线性扩展上运行的可逆 MH 链 + 其混合诊断"，
不是"严格的后验采样"。任何"量化不确定性"的解释都必须先通过诊断量检查。

------------------------------------------------------------------------------
数学定义
------------------------------------------------------------------------------
考虑物种树内部节点上的**偏序（partial order）**：对每条树边 ``(父, 子)``（两端均为
内部节点），要求父节点在子节点之前发生。一个**线性扩展（linear extension）**即该偏序
的一个全序（满足所有树边约束的内部节点排序）。MaxTiC 的三种启发式（greedy / mixing /
局部搜索）产出的排序都是该偏序的合法线性扩展。

把目标函数 ``E(order) = value(order, edge, edge_keys)`` 视为能量（越小越优），在其上
定义 Boltzmann 型**目标权重**：

    π(order) ∝ exp(-E(order) / T)

提议移动：在当前线性扩展上均匀随机选取两个不同位置 ``i != j`` 并交换。交换是**对合**
（involution），故提议矩阵对称 ``Q(x→y) = Q(y→x)``；非法候选（不再是线性扩展）直接拒绝
（自环，不破坏细致平衡）；接受概率为标准 MH 比率

    A(order → order') = min( 1, exp( (E(order) - E(order')) / T ) )

在上述对称提议下细致平衡方程对 π 成立，因此**若链已收敛且已充分混合**，样本的经验分布
逼近 π。"已收敛且充分混合"必须由 :meth:`MCMCSampler.diagnostics` 的诊断量支持，不得假设。

能量计算复用 ``maxtic_next.ranking.value.ValueComputer``，与核心算法零漂移。
"""

import math
from typing import Callable, Dict, List, Optional, Sequence, Set, Tuple

from maxtic_next.config import MCMC_STATUS_NOTE
from maxtic_next.random_ import RandomWrapper, sample_two_distinct
from maxtic_next.ranking.value import ValueComputer
from maxtic_next.tree.tree import Tree

# 状态空间大小超过该值时不再精确枚举/位掩码计数（诊断量退化为下界）
_MAX_ENUMERABLE_STATES = 4096
_MAX_ENUMERABLE_NODES = 16

# 自检温度 = 能量跨度 / 该除数（经验值：既保证 π 明显非均匀，又保证有限长链
# 能覆盖全部"统计可见"状态；见 :func:`_choose_balance_temperature`）
BALANCE_SPREAD_DIVISOR = 3.5

# 自检判定容差 = max(给定 tolerance, BALANCE_MC_FACTOR × 批次均值法标准误)。
# 取 4 ≈ 双侧 6e-5 的宽松门限：链样本相关，纯固定容差会误判合法链。
BALANCE_MC_FACTOR = 4.0

__all__ = [
    "MCMCSampler",
    "mcmc_diagnostics",
    "effective_sample_size",
    "MCMC_STATUS_NOTE",
]


class MCMCSampler:
    """物种树偏序线性扩展上的可逆 Metropolis–Hastings 链（初步实现）。

    用法::

        sampler = MCMCSampler(tree, edge, edge_keys, rng, temperature=T,
                              initial_order=best_order,
                              burn_in=n_iter // 2, thin=max(1, n_iter // 100))
        samples = sampler.sample(n_iter)             # List[List[str]]（已抽稀）
        diag = sampler.diagnostics()                 # 混合诊断（必须展示）

    ``sample()`` 的返回序列是**马尔可夫链样本**，不是独立样本；任何均值/分位数
    解释都必须同时报告 :meth:`diagnostics` 的 ``effective_sample_size`` 与
    ``unique_sample_fraction``。
    """

    #: 报告 / stdout 统一使用的状态说明
    STATUS_NOTE = MCMC_STATUS_NOTE

    def __init__(
        self,
        tree: Tree,
        edge: Dict[str, float],
        edge_keys: List[str],
        rng: RandomWrapper,
        temperature: float = 0.01,
        initial_order: Optional[List[str]] = None,
        burn_in: int = 0,
        thin: int = 1,
    ) -> None:
        """初始化采样器。

        Args:
            tree: 物种树（内部节点标签在 bootstrap 字段），用于构造偏序。
            edge: 边权重字典（键 ``"donor,receptor"`` -> 权重），作为能量项。
            edge_keys: 参与能量计算的边键列表。
            rng: 统一随机性驱动（``RandomWrapper``），保证固定种子可复现。
            temperature: Boltzmann 温度 ``T``（必须 > 0）。温度远小于能量尺度时链
                实际不可约，请务必核对 :meth:`diagnostics`。
            initial_order: 初始线性扩展（需为合法排序）；``None`` 时内部拓扑排序生成。
            burn_in: 默认 burn-in 步数（``sample()`` 未显式传入时使用），须 >= 0。
            thin: 默认抽稀间隔，须 >= 1。

        Raises:
            ValueError: ``temperature <= 0``、``burn_in < 0`` 或 ``thin < 1``。
        """
        if temperature <= 0:
            raise ValueError(f"temperature 必须 > 0（Boltzmann 温度），收到 {temperature}。")
        _check_burn_in(burn_in)
        _check_thin(thin)
        self.tree = tree
        self.edge = edge
        self.edge_keys = edge_keys
        self.rng = rng
        self.temperature = float(temperature)
        self.burn_in = int(burn_in)
        self.thin = int(thin)
        self.value_computer = ValueComputer(edge, edge_keys)

        # 内部节点集合 + 标签 -> 节点 id 映射
        self.internal_nodes: List[str] = list(tree.internal_node_labels())
        self._label_to_id: Dict[str, int] = {}
        for nid in tree.get_nodes():
            if not tree.is_leaf(nid):
                self._label_to_id[tree.get_bootstrap(nid)] = nid

        # 预计算偏序（祖先 -> 后代 的所有内部节点有序对）
        self._precedence: Set[Tuple[str, str]] = _precedence_from_tree(tree)

        # 初始状态：使用给定合法排序，否则拓扑排序
        if initial_order is not None:
            self.current: List[str] = list(initial_order)
            if not self.is_valid(self.current):
                self.current = self._topological_order()
        else:
            self.current = self._topological_order()
        self.current_energy: float = self.value_computer.value(self.current)
        # 最近一次 sample() 的能量序列（与返回的排序对齐），供调用方统计
        self.last_energies: List[float] = []
        # 最近一次 sample() 抽稀**前**的整条链能量（用于 ESS / 相邻重复率）
        self.last_chain_energies: List[float] = []
        self.last_chain: List[Tuple[str, ...]] = []
        self.last_burn_in: int = 0
        # 最近一次 sample() 是否被 stop_check 提前打断
        self.last_stopped_early: bool = False

    # ------------------------------------------------------------------
    # 偏序 / 合法性
    # ------------------------------------------------------------------
    # 合法性判定直接使用 ``self._precedence`` 集合，不需要 ``_precedes(a, b)``
    # 这样的逐对辅助方法。
    def is_valid(self, order: Sequence[str]) -> bool:
        """判断 ``order`` 是否为合法的线性扩展（满足所有偏序约束）。

        Args:
            order: 内部节点标签排序。

        Returns:
            True 当且仅当 ``order`` 是内部节点集合的一个排列且尊重全部偏序关系。
        """
        if len(order) != len(self.internal_nodes):
            return False
        if set(order) != set(self.internal_nodes):
            return False
        index = {node: i for i, node in enumerate(order)}
        for a, b in self._precedence:
            if index[a] > index[b]:
                return False
        return True

    def state_space(self, max_states: int = _MAX_ENUMERABLE_STATES) -> Optional[List[List[str]]]:
        """精确枚举状态空间（全部线性扩展）；过大时返回 ``None``。

        枚举带 ``limit`` 早停：一旦确认"超过 ``max_states``"即返回 ``None``，不再
        为注定丢弃的巨大状态空间付指数级时间。
        """
        if len(self.internal_nodes) > _MAX_ENUMERABLE_NODES:
            return None
        ext = _enumerate_linear_extensions(
            self.internal_nodes, self._precedence, limit=max_states + 1
        )
        if len(ext) > max_states:
            return None
        return ext

    def _topological_order(self) -> List[str]:
        """通过 Kahn 式拓扑排序生成一个合法线性扩展（确定性：按默认节点顺序选点）。"""
        successors: Dict[str, List[str]] = {n: [] for n in self.internal_nodes}
        indeg: Dict[str, int] = {n: 0 for n in self.internal_nodes}
        for a, b in self._precedence:
            successors[a].append(b)
            indeg[b] += 1
        available = [n for n in self.internal_nodes if indeg[n] == 0]
        order: List[str] = []
        while available:
            nxt = available.pop(0)
            order.append(nxt)
            for s in successors[nxt]:
                indeg[s] -= 1
                if indeg[s] == 0:
                    available.append(s)
            available.sort(key=lambda x: self.internal_nodes.index(x))
        return order

    # ------------------------------------------------------------------
    # 提议移动（对称对合）
    # ------------------------------------------------------------------
    def _propose_swap(self, order: List[str], rng: RandomWrapper) -> Tuple[List[str], bool]:
        """均匀随机选取两个位置交换，返回 (候选排序, 是否合法线性扩展)。

        交换是对合：对候选再次交换同一对位置即回到原排序，故提议矩阵对称。
        """
        n = len(order)
        # n < 2 时无法交换，直接返回原序（合法但不改变状态）
        pair = sample_two_distinct(n, rng)
        if pair is None:
            return list(order), True
        i, j = pair
        if i > j:
            i, j = j, i
        candidate = list(order)
        candidate[i], candidate[j] = candidate[j], candidate[i]
        return candidate, self.is_valid(candidate)

    # ------------------------------------------------------------------
    # 采样
    # ------------------------------------------------------------------
    def sample(
        self,
        n_iter: Optional[int] = None,
        thin: Optional[int] = None,
        burn_in: Optional[int] = None,
        stop_check: Optional[Callable[[], bool]] = None,
    ) -> List[List[str]]:
        """运行 MH 链 :meth:`sample` 的 ``n_iter`` 步并返回抽稀后的排序列表。

        Args:
            n_iter: **记录阶段**的链步数（必须为 ``>= 1`` 的整数；``0`` 或负数
                直接抛 :class:`ValueError`，不再让调用方在
                ``ZeroDivisionError`` / ``min([])`` 上崩溃）。``None`` 时取 1000。
            thin: 抽稀间隔（``>= 1``）；``None`` 时取构造期的 ``self.thin``。
            burn_in: 记录前先推进并丢弃的步数（``>= 0``）；``None`` 时取构造值。
            stop_check: **协作式取消回调**（``() -> bool``）。在每次链步
                **之间**检查（burn-in 与记录阶段都检查，检查是循环体首条语句）；
                返回真值时立即停止推进，返回**已记录到当下**的样本（可能为空，
                由调用方决定如何报告）。默认 ``None``：与不传时逐字节一致。

        Returns:
            抽稀后的线性扩展列表 ``List[List[str]]``（长度 ``ceil(n_iter/thin)``，
            至少为 1；被提前取消时可短于此，甚至为空）。

        Raises:
            ValueError: ``n_iter`` 非正、``thin < 1`` 或 ``burn_in < 0``。
        """
        if n_iter is None:
            n_iter = 1000
        _check_iters(n_iter)
        if thin is None:
            thin = self.thin
        _check_thin(thin)
        if burn_in is None:
            burn_in = self.burn_in
        _check_burn_in(burn_in)
        self.last_burn_in = int(burn_in)
        self.last_stopped_early = False

        for _ in range(burn_in):
            if stop_check is not None and stop_check():
                # 取消：burn-in 也停下来（少烧链只是样本更少，不影响正确性）
                self.last_stopped_early = True
                break
            self._advance()

        chain: List[Tuple[str, ...]] = []
        chain_energies: List[float] = []
        for _ in range(n_iter):
            if stop_check is not None and stop_check():
                self.last_stopped_early = True
                break
            state, energy = self._advance()
            chain.append(tuple(state))
            chain_energies.append(energy)

        thinned = [list(chain[i]) for i in range(0, len(chain), thin)]
        self.last_chain = chain
        self.last_chain_energies = chain_energies
        self.last_energies = [self.value_computer.value(s) for s in thinned]
        return thinned

    def _advance(self) -> Tuple[List[str], float]:
        """推进一条链步（一次对称提议 + MH 接受），返回当前状态与能量。"""
        candidate, valid = self._propose_swap(self.current, self.rng)
        if valid:
            cand_energy = self.value_computer.value(candidate)
            # 最小化能量：接受率 = min(1, exp((E_old - E_new)/T))
            delta = (self.current_energy - cand_energy) / self.temperature
            accept = (delta >= 0.0) or (self.rng.random() < _safe_exp(delta))
            if accept:
                self.current = candidate
                self.current_energy = cand_energy
        return self.current, self.current_energy

    # ------------------------------------------------------------------
    # 诊断（混合情况必须与统计量一并展示）
    # ------------------------------------------------------------------
    def diagnostics(
        self,
        samples: Optional[Sequence[Sequence[str]]] = None,
        energies: Optional[Sequence[float]] = None,
    ) -> Dict:
        """返回最近一次采样的混合诊断量（供 stdout / 稳健性摘要 / HTML 展示）。

        Args:
            samples: 待诊断的排序序列；``None`` 时用 ``sample()`` 缓存的链。
            energies: 与 ``samples`` 对齐的能量；``None`` 时内部重算。

        Returns:
            见 :func:`mcmc_diagnostics`；额外含 ``temperature`` / ``burn_in`` /
            ``thin`` / ``status_note``。当状态空间可精确枚举时，附
            ``state_space_size``、``states_visited`` 与 ``state_coverage_fraction``
            ——链"冻结"时 ``state_coverage_fraction`` 会显著小于 1。
        """
        if samples is None:
            samples = self.last_chain
            energies = self.last_chain_energies
        if energies is None:
            energies = [self.value_computer.value(list(s)) for s in samples]
        diag = mcmc_diagnostics(samples, energies)
        diag["temperature"] = self.temperature
        diag["burn_in"] = self.last_burn_in
        diag["status_note"] = MCMC_STATUS_NOTE
        space = self.state_space()
        if space is not None:
            all_states = {tuple(s) for s in space}
            visited = {tuple(s) for s in samples}
            diag["state_space_size"] = len(all_states)
            diag["states_visited"] = len(visited & all_states)
            diag["state_coverage_fraction"] = diag["states_visited"] / float(len(all_states))
        else:
            diag["state_space_size"] = None
            diag["state_coverage_fraction"] = None
        return diag

    def energy_histogram(self, samples: Sequence[Sequence[str]]) -> Dict[float, int]:
        """统计样本能量直方图（用于诊断目标分布形状）。"""
        hist: Dict[float, int] = {}
        for s in samples:
            e = self.value_computer.value(list(s))
            hist[e] = hist.get(e, 0) + 1
        return hist

    # ------------------------------------------------------------------
    # 自检：独立参考分布下的必要非充分验证
    # ------------------------------------------------------------------
    def self_test_detailed_balance(self, n_iter: int = 20000, tolerance: float = 0.05) -> bool:
        """在小规模合成实例上做一次**有判别力**的经验自检（必要非充分）。

        等价于 ``self_test_detailed_balance_report(...)["passed"]``。语义必须被
        理解为**必要而非充分**条件（见 :meth:`self_test_detailed_balance_report`）。
        """
        return bool(
            self.self_test_detailed_balance_report(n_iter=n_iter, tolerance=tolerance)["passed"]
        )

    def self_test_detailed_balance_report(
        self, n_iter: int = 20000, tolerance: float = 0.05, temperature: Optional[float] = None
    ) -> Dict:
        """用**独立参考实现**经验验证 MH 细致平衡，并返回可展示的诊断字典。

        与旧实现（2 个状态、彼此一次交换即达、且参考 π 与链共用同一
        ``_safe_exp`` / ``ValueComputer``，故在 T=20 下"接受一切合法提议"的坏链
        也能通过）不同，本自检：

        1. **独立参考**：自行回溯枚举合成偏序（7 个内部节点 → 80 个线性扩展）的
           全部线性扩展，用**独立编写**的能量函数（不经过 ``ValueComputer`` /
           ``value()``）与**裸** ``math.exp``（不经过 ``_safe_exp``）计算目标权重；
           并用另一套**位掩码动态规划**独立计数交叉校验枚举结果。
        2. **有判别力的温度**：在温度网格上自动挑选使
           ``TV(π, 状态间均匀分布) > tolerance`` 的温度——即"只看合法性、不看能量"
           的坏链必然被拒绝的工作点；若不存在这样的温度，自检直接判**不通过**
           （无判别力的自检没有意义）。
        3. **冻结链检查**：报告 ``state_coverage_fraction``（被访问状态数 / 全部状态数）
           与 ``visible_states_coverage_fraction``（被访问的"统计可见"状态数，即
           ``π >= 1/n_iter`` 的状态中被访问的比例）；覆盖不足直接判**不通过**，
           因此"链从未离开初始态"不可能蒙混过关。
        4. **多重统计量**：能量分组质量、逐约束的成对次序边缘率、平均能量的
           经验值与参考值偏差均须在容差内。

        Args:
            n_iter: 记录阶段的链长（burn-in 取其 1/4，不计入）。
            tolerance: 各统计量允许的最大绝对偏差（同时也是判别力门限）。
            temperature: 覆盖自动选择的自检温度（``None`` = 自动）。

        Returns:
            dict，含 ``passed``、``necessary_not_sufficient``（恒为 ``True`` 的
            标签）、``interpretation``（中文说明）、``failure_reasons``、
            ``temperature``、``state_space_size``、``states_visited``、
            ``state_coverage_fraction``（**链访问了多少比例的状态，冻结链在此露出**）、
            ``visible_states`` / ``visible_states_coverage_fraction``、
            ``power_vs_energy_blind``（与"能量盲链"的 TV，衡量自检判别力）、
            ``deviations``（每个统计量的经验值 / 独立参考值 / 偏差 / 允许上限）、
            ``max_energy_group_deviation``、``max_pair_deviation``、
            ``mean_energy_reference`` / ``mean_energy_experienced``、
            ``chain_diagnostics``（同一份链的 :meth:`diagnostics`）。

        注意：``passed=True`` 只说明"在这个小规模实例、这个温度、这条链长下，
        链的经验分布与独立参考一致且未冻结"。**它不能**证明采样器在真实数据、
        真实温度（尤其远小于能量尺度的低温）下已收敛——那必须依靠
        :meth:`diagnostics` 的 ``unique_sample_fraction`` /
        ``effective_sample_size`` / ``state_coverage_fraction``。
        """
        tree = Tree()
        tree.read_newick(_BALANCE_TEST_NEWICK)
        edge = dict(_BALANCE_TEST_EDGE)
        keys = sorted(edge.keys())
        nodes = list(tree.internal_node_labels())
        precedence = _precedence_from_tree(tree)

        # (1) 独立参考：回溯枚举 + 独立能量函数 + 裸 math.exp
        extensions = [tuple(e) for e in _enumerate_linear_extensions(nodes, precedence)]
        n_states = len(extensions)
        report: Dict = {
            "n_iter": n_iter,
            "tolerance": tolerance,
            "state_space_size": n_states,
            "necessary_not_sufficient": True,
        }
        if n_states == 0:
            report.update({"passed": False, "interpretation": "合成实例枚举为空：自检无法执行。"})
            return report
        # 交叉校验：位掩码 DP 独立计数必须与回溯枚举一致
        dp_count = _count_linear_extensions(nodes, precedence)
        report["dp_count_matches_enumeration"] = dp_count == n_states
        # 交叉校验：采样器的合法性判定必须恰好吃掉全部枚举态
        sampler = MCMCSampler(
            tree,
            edge,
            keys,
            RandomWrapper(20240817),
            temperature=1.0,
            initial_order=list(extensions[0]),
        )
        report["enumeration_agrees_with_is_valid"] = all(
            sampler.is_valid(list(ext)) for ext in extensions
        )

        ref_energy = {tuple(ext): _reference_energy(ext, edge) for ext in extensions}
        temp = temperature if temperature else _choose_balance_temperature(ref_energy, tolerance)
        report["temperature"] = temp
        weights = {s: math.exp(-e / temp) for s, e in ref_energy.items()}
        Z = sum(weights.values())
        target = {s: w / Z for s, w in weights.items()}

        # (2) 判别力：与"只看合法性、不看能量"（状态间均匀）的 TV 必须超过容差
        uniform = {tuple(s): 1.0 / n_states for s in extensions}
        tv_vs_blind = 0.5 * sum(abs(target[s] - uniform[s]) for s in uniform)
        report["power_vs_energy_blind"] = tv_vs_blind
        if not temperature and tv_vs_blind <= tolerance:
            report.update(
                {
                    "passed": False,
                    "interpretation": (
                        f"自检在该温度下没有判别力（与能量盲链的 TV="
                        f"{tv_vs_blind:.4f} <= 容差 {tolerance}）：即使采样器完全忽略"
                        "能量也能通过，故判定不通过。"
                    ),
                }
            )
            return report

        # (3) 跑链：带 burn-in，避免"从最优态出发被冻住"蒙混过关。
        # 用 type(self) 而非硬编码 MCMCSampler —— 自检检验的必须是"被调用的那个类"。
        rng = RandomWrapper(12345)
        chain_sampler = type(self)(
            tree,
            edge,
            keys,
            rng,
            temperature=temp,
            initial_order=list(extensions[0]),
            burn_in=max(1, n_iter // 4),
        )
        samples = chain_sampler.sample(n_iter, thin=1)
        counts: Dict[Tuple[str, ...], int] = {}
        for s in samples:
            k = tuple(s)
            counts[k] = counts.get(k, 0) + 1
        total = len(samples)

        visited = set(counts)
        unknown_states = [k for k in visited if k not in ref_energy]
        visible = {s for s in extensions if target[s] >= 1.0 / total}
        cov_all = len(visited) / float(n_states)
        cov_vis = (len(visited & visible) / float(len(visible))) if visible else 1.0
        report.update(
            {
                "states_visited": len(visited),
                "chain_states_outside_state_space": len(unknown_states),
                "state_coverage_fraction": cov_all,
                "visible_states": len(visible),
                "visible_states_visited": len(visited & visible),
                "visible_states_coverage_fraction": cov_vis,
            }
        )

        # (4) 统计量偏差：能量分组质量 + 逐约束成对边缘率 + 平均能量。
        # 判据为"|经验 - 参考| <= max(tolerance, mc_factor * 蒙特卡洛标准误)"：
        # 马尔可夫链样本相关，固定容差要么把合法链误判为失败、要么失去意义，
        # 故用批次均值法（batch means）估计各统计量的抽样误差。
        mc_factor = BALANCE_MC_FACTOR
        series_chain = [tuple(s) for s in samples]
        deviations: Dict[str, Dict[str, float]] = {}

        def _check(name: str, series: List[float], ref: float) -> float:
            emp = sum(series) / float(len(series))
            dev = abs(emp - ref)
            limit = max(tolerance, mc_factor * _batch_standard_error(series))
            deviations[name] = {
                "empirical": emp,
                "reference": ref,
                "deviation": dev,
                "limit": limit,
                "passed": dev <= limit,
            }
            return dev

        ref_groups: Dict[float, float] = {}
        for state in extensions:
            ref_groups[ref_energy[state]] = ref_groups.get(ref_energy[state], 0.0) + target[state]
        group_dev = 0.0
        for e, p in sorted(ref_groups.items()):
            indicator = [1.0 if ref_energy.get(k) == e else 0.0 for k in series_chain]
            group_dev = max(group_dev, _check(f"energy_group[{e}]", indicator, p))
        report["max_energy_group_deviation"] = group_dev
        report["energy_group_reference"] = {k: round(v, 6) for k, v in sorted(ref_groups.items())}
        report["energy_group_experienced"] = {
            k: deviations[f"energy_group[{k}]"]["empirical"] for k in sorted(ref_groups)
        }

        pair_dev = 0.0
        for key in keys:
            a, _, b = key.partition(",")
            ref_p = sum(target[s] for s in extensions if _index_of(s, a) < _index_of(s, b))
            indicator = [1.0 if _index_of(k, a) < _index_of(k, b) else 0.0 for k in series_chain]
            pair_dev = max(pair_dev, _check(f"pair[{key}]", indicator, ref_p))
        report["max_pair_deviation"] = pair_dev

        mean_series = [ref_energy.get(k, float("nan")) for k in series_chain]
        mean_ref = sum(ref_energy[s] * target[s] for s in extensions)
        mean_dev = _check("mean_energy", mean_series, mean_ref)
        report["mean_energy_reference"] = mean_ref
        report["mean_energy_experienced"] = deviations["mean_energy"]["empirical"]
        report["deviations"] = deviations
        report["mc_error_factor"] = mc_factor
        del mean_dev

        # 混诊断（同一份链）：冻结链在此也会露出马脚
        report["chain_diagnostics"] = chain_sampler.diagnostics()

        reasons: List[str] = []
        if not report["dp_count_matches_enumeration"]:
            reasons.append("回溯枚举与位掩码 DP 计数不一致（状态空间参考不可信）")
        if not report["enumeration_agrees_with_is_valid"]:
            reasons.append("采样器的 is_valid 拒绝了部分合法线性扩展")
        if unknown_states:
            reasons.append(
                f"链产出了 {len(unknown_states)} 个不在枚举状态空间中的排序（提议/合法性实现有误）"
            )
        failed = [k for k, v in deviations.items() if not v["passed"]]
        for name in failed:
            entry = deviations[name]
            reasons.append(
                f"{name}：经验 {entry['empirical']:.4f} vs 独立参考 {entry['reference']:.4f}"
                f"（偏差 {entry['deviation']:.4f} > 允许 {entry['limit']:.4f}）"
            )
        if cov_vis < 1.0:
            reasons.append(
                f"链未访问全部统计可见状态（{len(visited & visible)}/{len(visible)}"
                f"，全部状态覆盖 {cov_all:.3f}）：链已冻结或不可约"
            )
        if tv_vs_blind <= tolerance:
            reasons.append(f"自检无判别力（与能量盲链 TV={tv_vs_blind:.4f} <= {tolerance}）")

        report["failure_reasons"] = reasons
        report["passed"] = not reasons
        report["interpretation"] = (
            "通过：经验分布与该实例上独立算出的参考权重一致，且链覆盖了全部统计可见状态。"
            if report["passed"]
            else "不通过：" + "；".join(reasons)
        )
        if report["passed"]:
            report["interpretation"] += (
                "（必要非充分：本自检只覆盖 7 节点合成实例与自动选定的温度，"
                "不能推断真实数据/低温长链已收敛；须同时查看 diagnostics。）"
            )
        return report

    # 向后兼容别名（旧名保留，语义同上）
    def self_test_report(self, n_iter: int = 20000, tolerance: float = 0.05) -> Dict:
        """:meth:`self_test_detailed_balance_report` 的简写别名。"""
        return self.self_test_detailed_balance_report(n_iter=n_iter, tolerance=tolerance)


# ----------------------------------------------------------------------
# 参数校验（不得 ZeroDivisionError / min([])）
# ----------------------------------------------------------------------
def _check_iters(n_iter) -> None:
    """``iters`` 必须是正整数。"""
    if not isinstance(n_iter, int) or isinstance(n_iter, bool):
        raise TypeError(f"iters（链步数）必须是整数，收到 {type(n_iter).__name__}。")
    if n_iter <= 0:
        raise ValueError(
            f"iters（链步数）必须是正整数（>= 1），收到 {n_iter}。"
            "空链没有样本，任何 mean/min/ESS 都不成立。"
        )


def _check_thin(thin) -> None:
    """``thin`` 必须是 >= 1 的整数。"""
    if not isinstance(thin, int) or isinstance(thin, bool):
        raise TypeError(f"thin（抽稀间隔）必须是整数，收到 {type(thin).__name__}。")
    if thin < 1:
        raise ValueError(f"thin（抽稀间隔）必须 >= 1，收到 {thin}。")


def _check_burn_in(burn_in) -> None:
    """``burn_in`` 必须是 >= 0 的整数。"""
    if not isinstance(burn_in, int) or isinstance(burn_in, bool):
        raise TypeError(f"burn_in 必须是整数，收到 {type(burn_in).__name__}。")
    if burn_in < 0:
        raise ValueError(f"burn_in 必须 >= 0，收到 {burn_in}。")


def _index_of(order: Sequence[str], label: str) -> int:
    """``label`` 在 ``order`` 中的位置；不存在返回 -1（不参与比较）。"""
    for i, x in enumerate(order):
        if x == label:
            return i
    return -1


def _batch_standard_error(series: Sequence[float], n_batches: int = 10) -> float:
    """批次均值法（batch means）估计链样本**均值**的标准误。

    把链切成 ``n_batches`` 段，用段均值的样本方差除以段数开方：
    ``se = sd(batch_means) / sqrt(n_batches)``。相比 iid 公式
    ``sd/sqrt(n)``，它自动吸收了自相关带来的方差膨胀（
    马尔可夫链样本不是独立样本）。样本过少时返回 ``inf``（不放行任何判据）。
    """
    n = len(series)
    if n < 2 * n_batches:
        return float("inf")
    m = n // n_batches
    means: List[float] = []
    for b in range(n_batches):
        chunk = series[b * m : (b + 1) * m]
        if chunk:
            means.append(sum(chunk) / len(chunk))
    if len(means) < 2:
        return float("inf")
    mu = sum(means) / len(means)
    var = sum((x - mu) ** 2 for x in means) / (len(means) - 1)
    return math.sqrt(var / len(means))


# ----------------------------------------------------------------------
# 模块级工具
# ----------------------------------------------------------------------
def mcmc_diagnostics(
    samples: Sequence[Sequence[str]], energies: Optional[Sequence[float]] = None, max_lag: int = 500
) -> Dict:
    """MH 链混合诊断：唯一样本数 / 相邻重复率 / 有效样本量 (ESS)。

    Args:
        samples: 链上记录的排序（已抽稀或未抽稀均可）。
        energies: 与 ``samples`` 对齐的能量序列；``None`` 时全为 0（只统计重复）。
        max_lag: ESS 自相关估计的最大滞后。

    Returns:
        dict，含 ``n_samples``、``n_unique_samples``、``unique_sample_fraction``、
        ``adjacent_duplicate_fraction``（相邻记录完全相同的比例，1.0 = 链完全冻结）、
        ``max_chain_run_fraction``（最长同态连续段占样本数比例）、
        ``effective_sample_size``（对能量序列的初始正段自相关 ESS）、
        ``ess_fraction``、``energy_mean``、``energy_min``、``energy_sd``、
        ``independent_samples``（恒为 ``False`` 的诚实标签）。
    """
    n = len(samples)
    out: Dict = {"n_samples": n, "independent_samples": False}
    if n == 0:
        out.update(
            {
                "n_unique_samples": 0,
                "unique_sample_fraction": 0.0,
                "adjacent_duplicate_fraction": 0.0,
                "max_chain_run_fraction": 0.0,
                "effective_sample_size": 0.0,
                "ess_fraction": 0.0,
                "energy_mean": None,
                "energy_min": None,
                "energy_sd": None,
            }
        )
        return out
    keys = [tuple(s) for s in samples]
    n_unique = len(set(keys))
    dup = sum(1 for i in range(1, n) if keys[i] == keys[i - 1]) / float(n - 1) if n > 1 else 0.0
    longest = run = 1
    for i in range(1, n):
        run = run + 1 if keys[i] == keys[i - 1] else 1
        longest = max(longest, run)
    out["n_unique_samples"] = n_unique
    out["unique_sample_fraction"] = n_unique / float(n)
    out["adjacent_duplicate_fraction"] = dup
    out["max_chain_run_fraction"] = longest / float(n)

    if energies is None:
        energies = [0.0] * n
    energies = [float(e) for e in energies]
    mean_e = sum(energies) / n
    out["energy_mean"] = mean_e
    out["energy_min"] = min(energies)
    out["energy_sd"] = math.sqrt(sum((e - mean_e) ** 2 for e in energies) / n)
    ess = effective_sample_size(energies, max_lag=max_lag)
    out["effective_sample_size"] = ess
    out["ess_fraction"] = ess / float(n)
    return out


def effective_sample_size(values: Sequence[float], max_lag: int = 500) -> float:
    """由初始正段自相关估计有效样本量 ``ESS = n / (1 + 2 Σ ρ_k)``（Geyer 型截断）。

    独立样本 ESS = n；强自相关（如低温冻结链）ESS → 1，从而把"看着很多样本、
    实际只有一条解"这一事实量化出来。
    """
    n = len(values)
    if n < 2:
        return float(n)
    mean = sum(values) / n
    devs = [v - mean for v in values]
    var0 = sum(d * d for d in devs) / n
    if var0 <= 0.0:
        return 1.0  # 常数链：没有任何有效信息
    tau = 1.0
    limit = min(max_lag, n - 1)
    for k in range(1, limit + 1):
        cov = sum(devs[i] * devs[i + k] for i in range(n - k)) / n
        rho = cov / var0
        if rho <= 0.0:
            break  # 初始正段估计：首个非正自相关处截断
        tau += 2.0 * rho
    return max(1.0, min(float(n), n / tau))


def _choose_balance_temperature(
    ref_energy: Dict[Tuple[str, ...], float], tolerance: float
) -> float:
    """在温度上挑选"有判别力且链不至于冻结"的自检工作点。

    判别力 = ``TV(π, 状态间均匀分布) > tolerance``（否则"只看合法性、不看能量"的
    坏链也能通过，见：旧自检在 T=20 时正是如此）；
    可混合性 = 最大单态权重 ``< 0.5``（否则一条有限长的链必然漏掉大量统计可见态）。

    初值取 ``T = 能量跨度 / BALANCE_SPREAD_DIVISOR``（经验上此温度既使 π 明显非均匀、
    又让 80 态实例在 2 万步内完全覆盖），然后按判据向"更冷 / 更热"方向各试若干步。
    """
    energies = list(ref_energy.values())
    n_states = max(2, len(energies))
    lo, hi = min(energies), max(energies)
    spread = (hi - lo) or max(1.0, abs(hi) or 1.0)
    temp = spread / BALANCE_SPREAD_DIVISOR
    for _ in range(12):
        weights = [math.exp(-(e - lo) / temp) for e in energies]
        Z = sum(weights) or 1.0
        p = [w / Z for w in weights]
        tv = 0.5 * sum(abs(pi - 1.0 / n_states) for pi in p)
        if tv <= tolerance:
            temp /= 2.0  # π 太接近均匀：更冷才有判别力
        elif max(p) >= 0.5:
            temp *= 2.0  # π 太集中：更热才可能覆盖状态空间
        else:
            return temp
    return temp


def _safe_exp(x: float) -> float:
    """数值安全的 exp：对过大负值返回 0.0，避免下溢异常。

    仅用于**接受率**（``exp(delta)`` 的上界是 1，无需精确下尾）。自检的参考分布
    刻意使用裸 ``math.exp``，以免与被检代码共享同一实现。
    """
    if x < -700.0:
        return 0.0
    return math.exp(x)


# ---- 自检用合成实例：7 个内部节点、恰 80 个线性扩展 -------------------
# 偏序：R < Z,W；Z < X,Y；W < U,V  =>  6!/(3*3) = 80 个线性扩展
_BALANCE_TEST_NEWICK = (
    "(((A:1.0,B:1.0)X:1.0,(C:1.0,D:1.0)Y:1.0)Z:1.0,"
    "((E:1.0,F:1.0)U:1.0,(G:1.0,H:1.0)V:1.0)W:1.0)R:1.0;"
)
# 只涉及**不可比**节点对 => 全部是信息性约束（不被树拓扑强制）
_BALANCE_TEST_EDGE = {"X,U": 1.0, "Y,V": 1.0, "Z,W": 1.0, "X,V": 0.5}


def _reference_energy(order: Sequence[str], edge: Dict[str, float]) -> float:
    """**独立**编写的目标函数（不经过 ``value()`` / ``ValueComputer``）。

    与主实现的算法路径刻意不同：用 ``dict.items()`` 遍历 + ``str.partition`` 拆键，
    以便自检的参考值不是"用被检代码算出来的"。
    """
    pos: Dict[str, int] = {}
    for i, lbl in enumerate(order):
        pos.setdefault(lbl, i)
    total = 0.0
    for key, w in edge.items():
        a, _, b = key.partition(",")
        ia, ib = pos.get(a, -1), pos.get(b, -1)
        if ia >= 0 and ib >= 0 and ia > ib:
            total += float(w)
    return total


def _precedence_from_tree(tree: Tree) -> Set[Tuple[str, str]]:
    """返回物种树内部节点偏序（祖先 -> 后代）的所有有序对。"""
    nodes = list(tree.internal_node_labels())
    label_to_id: Dict[str, int] = {}
    for nid in tree.get_nodes():
        if not tree.is_leaf(nid):
            label_to_id[tree.get_bootstrap(nid)] = nid
    prec: Set[Tuple[str, str]] = set()
    for a in nodes:
        for b in nodes:
            if a != b and tree.is_ancestor(label_to_id[a], label_to_id[b]):
                prec.add((a, b))
    return prec


def _enumerate_linear_extensions(
    nodes: List[str], precedence: Set[Tuple[str, str]], limit: Optional[int] = None
) -> List[List[str]]:
    """回溯枚举偏序的全部线性扩展（仅用于小规模自检/诊断）。

    Args:
        nodes: 节点标签列表。
        precedence: 偏序有序对集合。
        limit: 可选上限：一旦找到 ``limit`` 个即**停止扩展并返回**（长度 == ``limit``
            表示"至少这么多、可能更多"）。诊断路径用它避免为"反正要丢弃"的巨大
            状态空间付指数级时间（规模失控）。
    """
    result: List[List[str]] = []
    placed: Set[str] = set()

    def backtrack(current: List[str]) -> bool:
        if len(current) == len(nodes):
            result.append(list(current))
            return limit is not None and len(result) >= limit
        for n in nodes:
            if n in placed:
                continue
            # n 的所有前驱必须已放置
            blocked = any((p, n) in precedence and p not in placed for p in nodes)
            if blocked:
                continue
            placed.add(n)
            current.append(n)
            if backtrack(current):
                current.pop()
                placed.discard(n)
                return True
            current.pop()
            placed.discard(n)
        return False

    backtrack([])
    return result


def _count_linear_extensions(
    nodes: Sequence[str], precedence: Set[Tuple[str, str]]
) -> Optional[int]:
    """**另一套独立算法**（位掩码 DP）计数线性扩展，用于交叉校验枚举。

    节点数超过 ``_MAX_ENUMERABLE_NODES`` 时返回 ``None``（DP 复杂度 2^n·n）。
    """
    n = len(nodes)
    if n > _MAX_ENUMERABLE_NODES:
        return None
    idx = {lbl: i for i, lbl in enumerate(nodes)}
    preds = [0] * n
    for a, b in precedence:
        if a in idx and b in idx:
            preds[idx[b]] |= 1 << idx[a]
    full = (1 << n) - 1
    dp = [0] * (full + 1)
    dp[0] = 1
    for mask in range(full + 1):
        ways = dp[mask]
        if not ways:
            continue
        for k in range(n):
            if mask & (1 << k):
                continue
            if preds[k] & ~mask:  # 仍有前驱未放置
                continue
            dp[mask | (1 << k)] += ways
    return dp[full]
