"""目标函数 ``value`` 与增量计算接口。

等价原版 ``MaxTiC.value(order)``：给定一个节点排序 ``order``，统计所有"被违反"的约束边
（即 ``donor`` 在 ``receptor`` 之后的边）的权重和。这是反馈弧集（FAS）权重的直接度量，
值越小排序越优。
"""

from typing import Dict, List, Optional


def value(order: List[str], edge: Dict[str, float], edge_keys: List[str]) -> float:
    """计算排序 ``order`` 下被违反的约束权重和。等价原版 ``value()``。

    Args:
        order: 节点排序（内部节点 bootstrap 标签列表）。
        edge: 边权重字典（键 ``"donor,receptor"`` -> 权重）。
        edge_keys: 参与计算的边键列表（原版全局 ``edge_keys``）。

    Returns:
        被违反约束的权重总和。
    """
    index = {}
    for i in range(len(order)):
        index[order[i]] = i
    result = 0.0
    for e in edge_keys:
        sommets = e.split(",")
        if sommets[0] in index and sommets[1] in index and index[sommets[0]] > index[sommets[1]]:
            result = result + edge[e]
    return result


class ValueComputer:
    """``value`` 的面向对象封装，供局部搜索复用。

    原版中 ``value()`` 依赖全局 ``edge`` / ``edge_keys``；此处显式持有。
    局部搜索在每次"区间旋转交换"后通过 ``value()`` 计算目标值，与全量
    ``value`` 严格一致（零漂移）。若后续需要增量重算，可在此扩展，
    但当前局部搜索的交换是"区间旋转"而非"单点交换"，并不直接复用
    单点增量接口。
    """

    def __init__(
        self, edge: Dict[str, float], edge_keys: List[str], total: Optional[float] = None
    ) -> None:
        self.edge = edge
        self.edge_keys = edge_keys
        self.total = total if total is not None else sum(edge.values())

    def value(self, order: List[str]) -> float:
        """计算给定排序的目标值（等价全量 ``value``）。"""
        return value(order, self.edge, self.edge_keys)


class IncrementalValueComputer:
    """增量目标函数计算器。

    在局部搜索的"区间旋转交换"邻域中，仅需计算被移动元素与被跨越元素之间
    的约束变化，无需全量遍历所有边。

    交换语义：将位置 ``a`` 的元素移动到位置 ``b``（a < b），位置 a+1..b 的元素
    左移一位。唯一改变相对顺序的是被移动元素 x 与位置 a+1..b 的元素 y 之间的关系：

    * 边 (x, y)：移动前 x 在 y 前（不违反），移动后 x 在 y 后（违反）→ delta += w
    * 边 (y, x)：移动前 y 在 x 后（违反），移动后 y 在 x 前（不违反）→ delta -= w

    复杂度：每次 ``trial_value`` / ``commit_move`` 为 O(b - a)，远优于全量 O(|E|)。

    为防止浮点漂移，每 ``verify_interval`` 次提交后自动用全量 ``value()``
    校正一次（默认 1000 次）。校正不消耗随机数，不影响搜索路径。
    """

    def __init__(
        self,
        edge: Dict[str, float],
        edge_keys: List[str],
        initial_order: List[str],
        verify_interval: int = 1000,
    ) -> None:
        """初始化增量计算器。

        Args:
            edge: 边权重字典。
            edge_keys: 参与计算的边键列表。
            initial_order: 初始排序。
            verify_interval: 浮点校正间隔（提交次数）。
        """
        self.edge = edge
        self.edge_keys = edge_keys
        self._verify_interval = verify_interval
        self._commit_count = 0

        # 当前排序与索引映射
        self._order: List[str] = list(initial_order)
        self._index: Dict[str, int] = {label: i for i, label in enumerate(initial_order)}
        # 全量计算初始值（保证与 value() 一致）
        self._current_value: float = value(initial_order, edge, edge_keys)

    @property
    def current_value(self) -> float:
        """当前排序的目标值。"""
        return self._current_value

    @property
    def order(self) -> List[str]:
        """当前排序的只读视图（返回副本）。"""
        return list(self._order)

    def _compute_delta(self, a: int, b: int) -> float:
        """计算将位置 a 的元素移到位置 b 的目标值增量。

        Args:
            a: 源位置（a < b）。
            b: 目标位置。

        Returns:
            目标值变化量（正值表示变差，负值表示变优）。
        """
        x = self._order[a]
        delta = 0.0
        # 遍历位置 a+1..b 的元素（被跨越的元素）
        for p in range(a + 1, b + 1):
            y = self._order[p]
            # 边 (x, y)：移动前 x 在 y 前不违反，移动后 x 在 y 后违反
            key_xy = x + "," + y
            if key_xy in self.edge:
                delta += self.edge[key_xy]
            # 边 (y, x)：移动前 y 在 x 后违反，移动后 y 在 x 前不违反
            key_yx = y + "," + x
            if key_yx in self.edge:
                delta -= self.edge[key_yx]
        return delta

    def trial_value(self, a: int, b: int) -> float:
        """计算假设移动后的目标值（不修改内部状态）。

        Args:
            a: 源位置（a < b）。
            b: 目标位置。

        Returns:
            移动后的目标值。
        """
        return self._current_value + self._compute_delta(a, b)

    def commit_move(self, a: int, b: int) -> float:
        """提交移动：更新排序、索引、目标值。

        Args:
            a: 源位置（a < b）。
            b: 目标位置。

        Returns:
            移动后的目标值。
        """
        delta = self._compute_delta(a, b)
        # 执行区间旋转：元素 a 移到 b，a+1..b 左移一位
        x = self._order[a]
        for p in range(a, b):
            self._order[p] = self._order[p + 1]
        self._order[b] = x
        # 更新索引映射（仅 a..b 范围内的位置改变）
        for p in range(a, b + 1):
            self._index[self._order[p]] = p
        # 更新目标值
        self._current_value += delta
        # 浮点校正
        self._commit_count += 1
        if self._commit_count % self._verify_interval == 0:
            self._current_value = value(self._order, self.edge, self.edge_keys)
        return self._current_value

    def full_recompute(self) -> float:
        """全量重算当前排序的目标值（用于校正浮点漂移或手动验证）。"""
        self._current_value = value(self._order, self.edge, self.edge_keys)
        return self._current_value
