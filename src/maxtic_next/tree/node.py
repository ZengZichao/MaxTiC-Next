"""树节点数据结构。

等价原版 ``script_tree.py`` 的节点槽位：``[name, parent, children, length, isdup, species, bootstrap, bppnumber]``。
为清晰与可扩展性，使用 ``dataclass`` 表达同样的字段；``id`` 为节点在 ``Tree`` 中的整数标识，
内部节点的"节点标识"为 ``bootstrap`` 字符串（与原版 ``getBootstrap`` 对齐），叶子为 ``name``。
"""

from dataclasses import dataclass, field
from typing import List, Optional


@dataclass
class TreeNode:
    """等价于 ``script_tree`` 节点槽位的轻量节点。

    字段顺序与原版 ``[name, parent, children, length, isdup, species, bootstrap, bppnumber]`` 一一对应。
    """

    id: int = 0  # 节点在 Tree 中的整数 id（原版 dict 的 key）
    name: str = ""  # 槽位 0：叶子名（内部节点此字段为占位 "N{id}"）
    parent: int = -1  # 槽位 1：父节点 id，-1 表示根
    children: List[int] = field(default_factory=list)  # 槽位 2：子节点 id 列表
    length: float = 0.0  # 槽位 3：分支长度
    isdup: str = ""  # 槽位 4：是否为复制节点（"D" 表示 duplication）
    species: str = ""  # 槽位 5：物种注释（NHX S=）
    bootstrap: str = ""  # 槽位 6：bootstrap / 内部节点标签
    bppnumber: int = -1  # 槽位 7：bpp 编号（原版内部使用）
    nd: Optional[str] = None  # NHX ND= 注释（可选）

    def is_leaf(self) -> bool:
        """是否为叶子节点（无子节点）。"""
        return len(self.children) == 0

    def is_root(self) -> bool:
        """是否为根节点（父为 -1）。"""
        return self.parent == -1
