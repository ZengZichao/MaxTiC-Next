"""近优解集合与稳健性/敏感性摘要。

本模块在**局部搜索**（``maxtic_next.ranking.local_search``）过程中收集被访问的近优排序，
并统计：

* 每节点在排序中的位置分布（``node_position_distribution``）；
* 成对节点的次序频率（``pairwise_order_frequency``，即 a 在 b 之前的次数）。

输出**仅**称为"基于局部搜索访问解的稳健性/敏感性摘要"（术语铁律）：
本模块**严禁**使用"置信区间 / 后验概率 / 后验分布"等用语——这些仅允许出现在
``robustness/mcmc.py`` 的 MCMC 报告里。局部搜索是一种 Metropolis 邻域采样，其访问
分布**不是**后验分布。

实现说明：``add`` 仅做 O(n) 的去重与暂存（``order`` 转 ``tuple`` 作键），统计学量在
``summary()`` 时一次性计算（O(K·n²)，K 为保留的近优序数量），以保持热路径开销可控。

------------------------------------------------------------------------------
：口径同时给出"去重占比"和"真实访问计数/频率"
------------------------------------------------------------------------------
局部搜索只在**邻域**里走，``add`` 又只保留目标值最小的 top-K 个序，因此本摘要衡量的是
"在近优集合这一条件下的次序敏感性"，**不是**对数据重抽样的稳健度，更不是概率区间。
旧 ``summary()`` 只遍历 ``_visited.keys()``，把 ``add()`` 维护的访问计数 ``cnt`` 整个丢掉，
于是"某序被访问 100 次、另一序 1 次"被报告成两个等权样本（还把缺失方向省略而不是记 0）。
现在两套口径都显式给出：

* ``n_orders_collected``：**去重**后保留的近优序数；
* ``n_accesses_retained`` / ``n_accesses_total``：真实访问次数（保留集合内 / 累计）；
* ``pairwise_order_frequency``：访问次数加权的计数，且**正反两个方向都出现**，
  未观测方向记 ``0`` 而非省略；
* ``pairwise_order_proportion``：上述计数除以 ``n_accesses_retained``（正反相加 = 1）；
* ``constant_direction_pairs``：在保留集合中方向恒定为 100% 的有序对（含被树种系拓扑
  强制的祖先-后代平凡对）；``informative_pairwise_order_frequency`` 为**剔除这些平凡对**
  后的计数，供报告展示真正有信息量的近优分歧。
"""

from typing import Dict, List, Optional, Set, Tuple

from maxtic_next.config import DEFAULT_NEAR_OPTIMAL_TOP_K

#: 默认保留的近优序数量上限（单一事实来源在 :mod:`maxtic_next.config`，
#: 现由 ``api.rank(top_k=...)`` / CLI ``--near-optimal-top-k`` / ``Ranker.run(top_k=)``
#: 全程可覆盖 —— 见 ``Ranker._run_impl`` 的 ``NearOptimalCollector(top_k=top_k)``）。
DEFAULT_TOP_K = DEFAULT_NEAR_OPTIMAL_TOP_K

#: ``add`` 累计"已见过的不同序"集合的大小上限（超过则只报下界，避免无界内存）。
_MAX_TRACKED_UNIQUE = 200000


class NearOptimalCollector:
    """在局部搜索过程中收集 top-K 去重近优排序，并汇总稳健性/敏感性统计。

    用法（在 ``optimisation_locale`` 中每接受一个新排序时调用 ``add(order, value)``）：
    ``collector.add(list(order), current)``。
    """

    def __init__(self, top_k: int = DEFAULT_TOP_K) -> None:
        """初始化收集器。

        Args:
            top_k: 保留的近优排序数量上限（按目标值升序保留最优的 K 个，值越小越优）。
        """
        if top_k < 1:
            raise ValueError(f"top_k 必须 >= 1，收到 {top_k}。")
        self.top_k = top_k
        # 去重近优序：tuple(order) -> (value, 访问计数)
        self._visited: Dict[Tuple[str, ...], Tuple[float, int]] = {}
        self.best_value: Optional[float] = None
        # 真实访问量（不随容量裁剪而减少）——
        self.total_accesses: int = 0
        self._ever_seen: Set[Tuple[str, ...]] = set()
        self._ever_seen_truncated: bool = False

    # ------------------------------------------------------------------
    # 收集接口
    # ------------------------------------------------------------------
    def add(self, order: List[str], value: float) -> None:
        """记录一个被局部搜索访问到的排序及其目标值（自动去重、容量裁剪）。

        Args:
            order: 节点排序（内部节点 bootstrap 标签列表）。
            value: 该排序的目标值（被违反约束权重和，越小越优）。
        """
        key = tuple(order)
        self.total_accesses += 1
        if len(self._ever_seen) < _MAX_TRACKED_UNIQUE:
            self._ever_seen.add(key)
        else:
            self._ever_seen_truncated = True
        if key in self._visited:
            _, cnt = self._visited[key]
            self._visited[key] = (value, cnt + 1)
        else:
            self._visited[key] = (value, 1)
        # 容量裁剪：保留目标值最小的 top_k 个（值越大越差，优先淘汰）。
        # 稳定淘汰（L3）：仅淘汰"最差（value 最大）"的序；若多个序并列最差，
        # 按插入顺序 FIFO 淘汰最早进入的那个。当前最优（最小 value）的序永不被淘汰，
        # 从而近优统计稳定、可复现（不会因字典遍历顺序波动而误删最优同值序）。
        if len(self._visited) > self.top_k:
            # `_visited` 的键是**排序元组**（旧标注误写成 str，导致 105/108/110 三处报错）。
            worst_key: Optional[Tuple[str, ...]] = None
            worst_val: float = float("inf")
            worst_pos: int = -1
            for pos, (key, (val, _cnt)) in enumerate(self._visited.items()):
                if worst_key is None or val > worst_val:
                    worst_key, worst_val, worst_pos = key, val, pos
                elif val == worst_val and pos < worst_pos:
                    # 等值最差：FIFO —— 淘汰更早插入（pos 更小）的序
                    worst_key, worst_pos = key, pos
            if worst_key is not None:
                del self._visited[worst_key]
        if self.best_value is None or value < self.best_value:
            self.best_value = value

    def __len__(self) -> int:
        """返回当前保留的（去重）近优序数量。"""
        return len(self._visited)

    # ------------------------------------------------------------------
    # 摘要（仅称"基于局部搜索访问解的稳健性/敏感性摘要"）
    # ------------------------------------------------------------------
    def summary(self, top_k: Optional[int] = None) -> Dict:
        """汇总稳健性/敏感性统计，返回 dict（**同时**给出去重占比与真实访问计数）。

        返回的 dict 以 ``title`` 字段显式标注为"基于局部搜索访问解的稳健性/敏感性摘要"，
        不含任何后验/置信区间用语（术语铁律）。

        Args:
            top_k: 本次**报告**只使用保留集合中目标值最小的前 ``top_k`` 个序
                （``None`` = 用构造期的 ``self.top_k``，即全部保留序）。
                取更小的值可让摘要聚焦头部近优解；不得大于已保留的序数时才生效。

        Returns:
            含以下键的 dict：
            * ``title``：摘要标题（术语铁律）；
            * ``n_orders_collected``：参与统计的（**去重**）近优序数量；
            * ``n_orders_retained``：容量裁剪后收集器实际保留的去重序数；
            * ``n_accesses_retained``：上述去重序的**真实访问次数之和**（重数保留）；
            * ``n_accesses_total``：``add()`` 被调用的累计次数（不受容量裁剪影响）；
            * ``n_unique_orders_total``：累计见过的不同序数量（超上限时为下界，
              并由 ``n_unique_orders_total_is_lower_bound`` 标注）；
            * ``top_k`` / ``top_k_reported``：保留上限 / 本次报告实际使用的上限；
            * ``best_value``：最优目标值（当前保留序中的最小 value）；
            * ``node_position_distribution``：节点 -> {位置: 访问次数}，
              **每个合法位置都出现**，未被访问的位置记 ``0`` 而非省略；
            * ``node_position_frequency``：同上但为访问次数占比；
            * ``pairwise_order_frequency``：有序对 ``(a, b)`` -> "a 在 b 之前"的
              **访问次数**（按重数加权），正反两个方向**都**给出键，未观测方向记 ``0``；
            * ``pairwise_order_proportion``：同口径的比例（``(a,b) + (b,a) == 1``）；
            * ``constant_direction_pairs``：比例恒为 1 的有序对列表
              ``[[a, b], ...]``（含树种系强制的祖先-后代平凡对）；
            * ``informative_pairwise_order_frequency``：剔除上述恒定对之后的计数
              （真正有信息量的近优分歧，报告表应优先使用它）；
            * ``is_resampling_robustness``：恒为 ``False`` —— 本摘要不是数据重抽样稳健度；
            * ``basis``：口径说明文本。
        """
        retained = list(self._visited.items())  # [(order, (value, cnt))]
        n_orders_retained = len(retained)
        if top_k is None:
            top_k_eff = self.top_k
        else:
            if top_k < 1:
                raise ValueError(f"top_k 必须 >= 1，收到 {top_k}。")
            top_k_eff = top_k
        # 报告口径：按目标值升序（值越小越优）取前 top_k_eff 个去重序
        reportable = sorted(retained, key=lambda kv: (kv[1][0],))[:top_k_eff]

        node_position: Dict[str, Dict[int, int]] = {}
        pair_freq: Dict[Tuple[str, str], int] = {}
        n_accesses = 0
        n_positions = max([len(o) for o, _ in reportable], default=0)

        # 1) 收集节点集合 / 真实访问次数
        for order, (_val, cnt) in reportable:
            n_accesses += cnt
            for node in order:
                node_position.setdefault(node, {})
        # 2) 补齐"节点 -> 每一个合法位置"的 0 基底（缺失位置记 0 而非省略）
        for node in node_position:
            for pos in range(n_positions):
                node_position[node][pos] = 0
        # 3) 补齐"成对 -> 正反两个方向"的 0 基底，再按访问次数累加
        all_pairs: Set[Tuple[str, str]] = set()
        for order, (_val, _cnt) in reportable:
            for i in range(len(order)):
                for j in range(i + 1, len(order)):
                    all_pairs.add((order[i], order[j]))
                    all_pairs.add((order[j], order[i]))
        for pair in all_pairs:
            pair_freq[pair] = 0
        for order, (_val, cnt) in reportable:
            for pos, node in enumerate(order):
                node_position[node][pos] += cnt
            for i in range(len(order)):
                for j in range(i + 1, len(order)):
                    pair_freq[(order[i], order[j])] += cnt

        denom = float(n_accesses) if n_accesses else 0.0
        pair_prop: Dict[Tuple[str, str], float] = {
            k: (v / denom if denom else 0.0) for k, v in pair_freq.items()
        }
        # "方向恒定"的对：在近优集合中 100% 同向（含树种系强制的祖先-后代平凡对）。
        # 正反两个朝向都算平凡，故一并从"有信息量"的表里剔除。
        unanimous: Set[Tuple[str, str]] = {
            k for k, p in pair_prop.items() if denom and p >= 1.0 - 1e-12
        }
        trivial_keys: Set[Tuple[str, str]] = set(unanimous) | {(b, a) for (a, b) in unanimous}
        constant_pairs: List[List[str]] = [[a, b] for (a, b) in sorted(unanimous)]
        informative: Dict[Tuple[str, str], int] = {
            k: v for k, v in pair_freq.items() if k not in trivial_keys
        }

        node_position_freq: Dict[str, Dict[int, float]] = {
            node: {pos: (c / denom if denom else 0.0) for pos, c in dist.items()}
            for node, dist in node_position.items()
        }

        return {
            "title": "基于局部搜索访问解的稳健性/敏感性摘要",
            "n_orders_collected": len(reportable),
            "n_orders_retained": n_orders_retained,
            "n_accesses_retained": n_accesses,
            "n_accesses_total": self.total_accesses,
            "n_unique_orders_total": len(self._ever_seen),
            "n_unique_orders_total_is_lower_bound": self._ever_seen_truncated,
            "top_k": self.top_k,
            "top_k_reported": top_k_eff,
            "best_value": self.best_value,
            "node_position_distribution": node_position,
            "node_position_frequency": node_position_freq,
            "pairwise_order_frequency": pair_freq,
            "pairwise_order_proportion": pair_prop,
            "constant_direction_pairs": constant_pairs,
            "informative_pairwise_order_frequency": informative,
            "is_resampling_robustness": False,
            "basis": (
                "口径：统计在'目标值最小的前 {k} 个去重近优排序'上进行，按每个序被局部"
                "搜索访问的真实次数加权；不是对输入数据重抽样的稳健度，也不是概率区间。"
            ).format(k=top_k_eff),
        }
