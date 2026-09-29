"""树模块公开接口。

注：``tree/cache.py`` 的 ``PathCache`` **不在排序热路径上**，也不受 ``--incremental`` /
``--checkpoint`` 控制。真正的加速来自 ``ranking/reachability.py`` 的
``ReachabilityMatrix``，它在 ``Ranker.run()`` 中**无条件启用**。
本模块仍导出 ``PathCache`` 以保持公开 import 面稳定，其适用边界见该模块自身的说明。
"""

from maxtic_next.tree.node import TreeNode
from maxtic_next.tree.tree import Tree
from maxtic_next.tree.cache import PathCache

__all__ = ["TreeNode", "Tree", "PathCache"]
