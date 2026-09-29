"""边构建与 ``edgeweights`` 原语。

等价原版：

* ``edgeweights(element, elements)``：对 ``element`` 的每条入边（``degre_entrant``），
  若其起点 ``s`` 落在候选集合 ``elements`` 中，则累加边 ``(s, element)`` 的权重。
* ``EdgeBuilder``：把解析后的约束（已按距离过滤）聚合为 ``edge`` 字典，并处理
  ``RANDOM_TYPE``（0/1/2）与 ``None`` 端点替换，等价原版约束解析主循环。
"""

import math
from typing import Collection, Dict, List

from maxtic_next.config import MAX_NUMBER, RANDOM_TYPE_CHOICES
from maxtic_next.constraints.constraint import ConstraintSet
from maxtic_next.random_ import RandomWrapper
from maxtic_next.tree.tree import Tree


def edgeweights(
    element: str,
    elements: Collection[str],
    degre_entrant: Dict[str, List[str]],
    edge: Dict[str, float],
) -> float:
    """等价原版 ``edgeweights(element, elements)``。

    **求和顺序 = ``degre_entrant[element]`` 的迭代顺序**，与 ``elements`` 的容器类型无关，
    因此把 ``elements`` 从列表换成集合不会改变浮点累加次序，输出逐字节一致
    （这正是 -1 敢改的原因；由 ``tests/test_output_equivalence.py`` 看守）。

    ``elements`` 接受任意容器，但**调用方应传 ``set``**：原版风格调用（每格重新切一片列表）
    会把成员判定退化成 O(len) 线性扫描，是默认排序路径的热点（量级可在仓库内
    ``benchmarks/run_benchmark.py`` 复现）。
    """
    result = 0.0
    for s in degre_entrant[element]:
        if s in elements:
            key = s + "," + element
            if key in edge:
                result = result + edge[key]
    return result


class EdgeBuilder:
    """从 ``ConstraintSet`` 构建算法所需的 ``edge`` 字典。"""

    def __init__(
        self, tree: Tree, cset: ConstraintSet, rng: RandomWrapper, random_type: int = 0
    ) -> None:
        self.tree = tree
        self.cset = cset
        self.rng = rng
        self.random_type = random_type

    def build(self) -> Dict[str, float]:
        """聚合约束为 ``edge`` 字典，等价原版约束解析主循环。

        处理 ``RANDOM_TYPE``（0=原样，1=保留节点随机化方向，2=完全随机化节点）与
        ``None`` 端点替换为根标签。权重累加到 ``edge[key]``，遇 ``MAX_NUMBER`` 哨兵则停止累加。

        ：``random_type`` 只接受 ``{0,1,2}``；旧实现对其它取值静默按 0 处理，
        使"打错随机化类型"表现为一份**未随机化**的假对照。
        """
        if self.random_type not in RANDOM_TYPE_CHOICES:
            raise ValueError(
                f"random_type 只能是 {list(RANDOM_TYPE_CHOICES)} 之一，"
                f"收到 {self.random_type!r}（其它取值过去被静默当作 0，即不做随机化）。"
            )
        edge: Dict[str, float] = {}
        internal_nodes = [str(x) for x in self.tree.internal_node_labels()]
        root_label = self.tree.get_bootstrap(self.tree.get_root())
        if self.random_type == 2 and len(set(internal_nodes)) < 2:
            raise ValueError(
                "random_type=2（完全随机化节点）至少需要物种树含 2 个内部节点，"
                f"当前只有 {len(set(internal_nodes))} 个，无法抽出两个不同节点。"
            )

        for c in self.cset.constraints:
            if not math.isfinite(c.weight) or c.weight < 0:
                raise ValueError(
                    f"约束 {c.donor!r}->{c.receptor!r} 的权重非法：{c.weight!r}"
                    "（必须是有限实数且 >= 0）。"
                )
            first = c.donor
            second = c.receptor
            # RANDOM_TYPE 处理（等价原版）
            if self.random_type == 1:
                if self.rng.random() < 0.5:
                    tmp = first
                    first = second
                    second = tmp
            elif self.random_type == 2:
                first = self.rng.choice(internal_nodes)
                second = first
                while second == first:
                    second = self.rng.choice(internal_nodes)
            if first == "None":
                first = root_label
            if second == "None":
                second = root_label
            key = first + "," + second
            if key not in edge:
                edge[key] = 0.0
            if not edge[key] >= MAX_NUMBER:
                edge[key] = edge[key] + c.weight
        return edge
