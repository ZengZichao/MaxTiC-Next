"""祖先 / 后代集合缓存。

**现状说明（，务必先读）**：本类是**独立的工具实现，未被排序热路径调用**。
``--incremental`` 启用的是 ``ranking/value.IncrementalValueComputer``，
``--checkpoint`` 启用的是 ``ranking/checkpoint.CheckpointManager``，二者都**不会**
启用 ``PathCache``；排序算法里真正的加速来自 ``ranking/reachability.py`` 的
``ReachabilityMatrix``，它在 ``Ranker.run()`` 中无条件启用。
注意：``PathCache`` **不在排序热路径上**，也不受 ``--incremental`` / ``--checkpoint``
控制；那两个开关驱动的是 ``ranking/reachability.py`` 与检查点逻辑。

祖先 / 后代集合（bitset）与 ``path`` 统一在 ``PathCache`` 缓存，避免重复遍历。
除 ``ancestor_set`` / ``path_on_tree`` 外还提供：

* ``descendant_set(node_id)``：返回以 ``node_id`` 为根的子树全部节点 id 集合；
* 可选 **numpy bitset** 表示（``ancestor_bitset`` / ``descendant_bitset``）：当
  ``numpy`` 可用时，用布尔数组表示集合，支持 O(1) 集合交/并/判祖先，便于大规模
  冲突检测向量化；不可用则**纯 Python 回退**到 ``frozenset``，接口一致。

注意：本缓存针对**树拓扑**的祖先关系（``is_ancestor``），用于性能加速；排序算法中
的 ``path(graph, a, b)``（见 ``maxtic_next.ranking.path``）操作的是**增广图**（树种系边
+ 被选中的约束边），由该模块直接对图做可达性判断，二者语义不同，请勿混淆。
"""

from typing import TYPE_CHECKING, Dict, FrozenSet, Optional

try:  # numpy 可选：缺失则纯 Python 回退
    import numpy as np

    _HAVE_NUMPY = True
except Exception:  # noqa: BLE001
    np = None  # type: ignore[assignment]  # 覆盖上面的模块名：可选依赖缺席时按 None 处理
    _HAVE_NUMPY = False

if TYPE_CHECKING:
    from maxtic_next.tree.tree import Tree


class PathCache:
    """树拓扑祖先 / 后代集合缓存（可选 numpy bitset 回退）。"""

    def __init__(self, tree: "Tree", use_bitset: Optional[bool] = None) -> None:
        """初始化缓存。

        Args:
            tree: 物种树。
            use_bitset: 是否使用 numpy bitset 表示。``None``（默认）表示"numpy 可用
                时自动启用，否则纯 Python 回退"。
        """
        self.tree = tree
        self._ancestor_sets: Dict[int, FrozenSet[int]] = {}
        self._descendant_sets: Dict[int, FrozenSet[int]] = {}
        self._use_bitset = _HAVE_NUMPY if use_bitset is None else (use_bitset and _HAVE_NUMPY)
        self._n_nodes = max(tree.get_nodes()) + 1 if tree.get_nodes() else 0
        # bitset 缓存：node_id -> np.ndarray(bool)
        self._ancestor_bits: Dict[int, "np.ndarray"] = {}
        self._descendant_bits: Dict[int, "np.ndarray"] = {}

    # ------------------------------------------------------------------
    # 纯 Python（frozenset）接口
    # ------------------------------------------------------------------
    def ancestor_set(self, node_id: int) -> FrozenSet[int]:
        """返回 ``node_id`` 的全部祖先节点 id 集合（不含自身）。"""
        if node_id not in self._ancestor_sets:
            anc: set = set()
            cur = node_id
            while not self.tree.is_root(cur):
                cur = self.tree.get_parent(cur)
                anc.add(cur)
            self._ancestor_sets[node_id] = frozenset(anc)
        return self._ancestor_sets[node_id]

    def descendant_set(self, node_id: int) -> FrozenSet[int]:
        """返回以 ``node_id`` 为根的子树全部节点 id 集合（含自身）。"""
        if node_id not in self._descendant_sets:
            desc: set = set()

            def _walk(n: int) -> None:
                desc.add(n)
                for c in self.tree.get_children(n):
                    _walk(c)

            _walk(node_id)
            self._descendant_sets[node_id] = frozenset(desc)
        return self._descendant_sets[node_id]

    def path_on_tree(self, a: int, b: int) -> bool:
        """判断在树拓扑上 ``a`` 是否为 ``b`` 的祖先（等价 ``is_ancestor(a, b)``）。"""
        return a in self.ancestor_set(b)

    # ------------------------------------------------------------------
    # numpy bitset 接口（缺失时回退到 frozenset 判定）
    # ------------------------------------------------------------------
    def ancestor_bitset(self, node_id: int) -> "np.ndarray":
        """返回 ``node_id`` 祖先的 numpy 布尔数组（bit i 为 True 表示 i 是祖先）。

        仅当 numpy 可用时有效；否则抛出 ``RuntimeError``（调用方应使用
        ``ancestor_set`` 回退）。
        """
        if not self._use_bitset:
            raise RuntimeError("numpy 不可用或未启用 bitset，请使用 ancestor_set()。")
        if node_id not in self._ancestor_bits:
            bits = np.zeros(self._n_nodes, dtype=bool)
            for a in self.ancestor_set(node_id):
                bits[a] = True
            self._ancestor_bits[node_id] = bits
        return self._ancestor_bits[node_id]

    def descendant_bitset(self, node_id: int) -> "np.ndarray":
        """返回以 ``node_id`` 为根的子树节点的 numpy 布尔数组（bit i 为 True 表示 i 在子树内）。"""
        if not self._use_bitset:
            raise RuntimeError("numpy 不可用或未启用 bitset，请使用 descendant_set()。")
        if node_id not in self._descendant_bits:
            bits = np.zeros(self._n_nodes, dtype=bool)
            for d in self.descendant_set(node_id):
                bits[d] = True
            self._descendant_bits[node_id] = bits
        return self._descendant_bits[node_id]

    @property
    def use_bitset(self) -> bool:
        """是否正在使用 numpy bitset 表示。"""
        return self._use_bitset
