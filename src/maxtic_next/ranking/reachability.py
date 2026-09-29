"""动态可达性矩阵（bitset 加速冲突检测）。

在贪婪启发式中，``path(graph, a, b)`` 每次做 DFS 判断可达性，复杂度 O(|V|)。
本模块维护增量更新的可达性矩阵，使 ``can_reach(a, b)`` 降为 O(1) 查表，
``add_edge(u, v)`` 增量更新为 O(|V|)（集合并运算，常数远小于 DFS 遍历）。

当 numpy 可用时，使用布尔数组（bitset）表示集合，``add_edge`` 的并运算可向量化；
否则纯 Python ``set`` 回退，接口一致。

注意：本矩阵维护的是**增广图**（树种系边 + 被选中的约束边）的可达关系，
与 ``tree/cache.py`` 的 ``PathCache``（仅树拓扑祖先关系）语义不同。
"""

from typing import Dict, List, Optional, Tuple

try:  # numpy 可选：缺失则纯 Python 回退
    import numpy as np

    _HAVE_NUMPY = True
except Exception:  # noqa: BLE001
    np = None  # type: ignore[assignment]  # 覆盖上面的模块名：可选依赖缺席时按 None 处理
    _HAVE_NUMPY = False


class ReachabilityMatrix:
    """动态可达性矩阵：增量维护增广图的可达关系，O(1) 查表替代 O(V) DFS。

    维护两组对称结构：
    * ``descendants[i]``：从节点 i 可达的所有节点集合（含 i 自身）；
    * ``ancestors[i]``：能到达节点 i 的所有节点集合（含 i 自身）。

    添加边 (u, v) 时，对每个能到达 u 的节点 x，将 v 的可达集合并入 x 的可达集合；
    对每个 v 可达的节点 y，将 u 的祖先集合并入 y 的祖先集合。
    """

    def __init__(
        self,
        labels: List[str],
        tree_edges: List[Tuple[str, str]],
        use_bitset: Optional[bool] = None,
    ) -> None:
        """初始化可达性矩阵。

        Args:
            labels: 所有节点标签列表（内部节点 bootstrap 标签 + 叶子名）。
            tree_edges: 初始树种系边列表，每条为 ``(parent_label, child_label)``。
            use_bitset: 是否使用 numpy bitset。``None``（默认）表示 numpy 可用时自动启用。
        """
        self._labels = list(labels)
        self._label_to_idx: Dict[str, int] = {label: i for i, label in enumerate(self._labels)}
        self._n = len(self._labels)
        self._use_bitset = _HAVE_NUMPY if use_bitset is None else (use_bitset and _HAVE_NUMPY)

        if self._use_bitset:
            # numpy bitset: descendants[i] = 布尔数组，bit j 为 True 表示 j 从 i 可达
            self._desc_bits: List["np.ndarray"] = [
                np.zeros(self._n, dtype=bool) for _ in range(self._n)
            ]
            self._anc_bits: List["np.ndarray"] = [
                np.zeros(self._n, dtype=bool) for _ in range(self._n)
            ]
            for i in range(self._n):
                self._desc_bits[i][i] = True
                self._anc_bits[i][i] = True
        else:
            # Python set 回退
            self._desc_sets: List[set] = [set([i]) for i in range(self._n)]
            self._anc_sets: List[set] = [set([i]) for i in range(self._n)]

        # 用树种系边初始化
        for u, v in tree_edges:
            self._add_edge_internal(u, v)

    def _add_edge_internal(self, u: str, v: str) -> None:
        """添加有向边 u -> v 并增量更新可达性矩阵。"""
        ui = self._label_to_idx[u]
        vi = self._label_to_idx[v]

        if self._use_bitset:
            # 复制避免别名（后续可能修改 desc_bits[vi] / anc_bits[ui]）
            anc_u_mask = self._anc_bits[ui].copy()
            desc_v_mask = self._desc_bits[vi].copy()
            anc_u_idx = np.nonzero(anc_u_mask)[0]
            desc_v_idx = np.nonzero(desc_v_mask)[0]
            # 对每个能到达 u 的 x，将 v 的后代集合并入 x 的后代集合
            for x in anc_u_idx:
                self._desc_bits[x] |= desc_v_mask
            # 对每个 v 可达的 y，将 u 的祖先集合并入 y 的祖先集合
            for y in desc_v_idx:
                self._anc_bits[y] |= anc_u_mask
        else:
            # 复制避免迭代中修改原始集合
            anc_u = set(self._anc_sets[ui])
            desc_v = set(self._desc_sets[vi])
            for x in anc_u:
                self._desc_sets[x] |= desc_v
            for y in desc_v:
                self._anc_sets[y] |= anc_u

    def can_reach(self, a: str, b: str) -> bool:
        """判断 b 是否从 a 可达（a == b 时返回 True，等价 ``path(graph, a, b)``）。

        Args:
            a: 起点标签。
            b: 终点标签。

        Returns:
            若 b 从 a 可达（含 a == b），返回 True。
        """
        ai = self._label_to_idx[a]
        bi = self._label_to_idx[b]
        if self._use_bitset:
            return bool(self._desc_bits[ai][bi])
        else:
            return bi in self._desc_sets[ai]

    def add_edge(self, u: str, v: str) -> None:
        """添加有向边 u -> v。

        调用方应先检查 ``can_reach(v, u)`` 以避免引入环。
        """
        self._add_edge_internal(u, v)

    @property
    def use_bitset(self) -> bool:
        """是否正在使用 numpy bitset 表示。"""
        return self._use_bitset
