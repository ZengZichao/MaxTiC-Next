"""约束数据结构。

等价原版 ``MaxTiC.py`` 中 ``edge[key] = weight`` 的单条加权有向约束
（``donor -> receptor``，即"donor 应早于 receptor"）。
每条约束携带 ``metadata``，至少含 ``family`` / ``support`` / ``distance``，
其中 ``distance`` 用于 ``--min-transfer-distance`` 按 phylogenetic distance 过滤。

``metadata["distance"]`` 的三态：正数=真实拓扑距离；
``0.0``=供体与受体同节点（自环，会被任何 ``d >= 0`` 阈值删掉）；
``None``=**距离未知**（文本输入无第 4 列，或上游端点标签不在物种树中），
按"无距离列"处理而保留。
"""

from dataclasses import dataclass, field
from typing import Dict, List, Optional


@dataclass
class Constraint:
    """一条加权有向时间约束：donor 应早于 receptor 发生。"""

    donor: str  # 供体（内部节点 bootstrap 标签）
    receptor: str  # 受体（内部节点 bootstrap 标签）
    weight: float = 1.0  # 约束权重，缺失默认 1.0
    metadata: Dict = field(default_factory=dict)  # 含 family / support / distance

    def to_edge_key(self) -> str:
        """返回边键 ``"donor,receptor"``，用于 ``edge`` 字典的键与输出文件。"""
        return f"{self.donor},{self.receptor}"


class ConstraintSet:
    """约束集合：承载解析、过滤、聚合为 ``edge`` 字典所需的所有接口。

    ``diagnostics``：上游适配器写入的
    机器可读诊断（口径标签 + 各类丢弃计数 + stderr 摘要来源）。文本约束路径同样
    会由 ``filter_by_distance`` 填入 ``distance_filter`` 段，供 ``ranker`` / ``api``
    如实上报"该选项被忽略"。
    """

    def __init__(
        self, constraints: Optional[List[Constraint]] = None, diagnostics: Optional[Dict] = None
    ) -> None:
        self.constraints: List[Constraint] = list(constraints) if constraints else []
        self.diagnostics: Dict = dict(diagnostics) if diagnostics else {}

    def add(self, constraint: Constraint) -> None:
        """追加一条约束。"""
        self.constraints.append(constraint)

    def __len__(self) -> int:
        return len(self.constraints)

    def __iter__(self):
        return iter(self.constraints)

    def filter_by_distance(self, d: float) -> Dict:
        """按 phylogenetic distance 过滤（对应原版 ``d=MIN_TRANSFER_DIST``）。

        保留规则：``distance`` 缺失（无距离列 / 端点无法在物种树上解析）或
        ``distance > d`` 的约束；丢弃 ``distance <= d`` 的约束。等价于原版解析时的
        ``len(words) <= 3 or float(words[3]) > MIN_TRANSFER_DIST``。

        ：本方法的真实行为是"无距离列 ⇒ 一条都不删"，而旧预检文案声称
        "会丢弃全部约束"。故此处把实际发生的情况记入 ``self.diagnostics``
        （``distance_filter`` 段）并返回同一结构，调用方必须据此报告
        ``ignored=True``（error 级"该选项被忽略"），而不是反过来猜。

        Args:
            d: 最小转移距离阈值。

        Returns:
            本次过滤的诊断字典。
        """
        kept: List[Constraint] = []
        dropped = 0
        missing = 0
        self_loop_dropped = 0
        for c in self.constraints:
            dist = c.metadata.get("distance")
            if dist is None:
                missing += 1
                kept.append(c)
            elif dist > d:
                kept.append(c)
            else:
                dropped += 1
                if c.donor == c.receptor:
                    self_loop_dropped += 1
        self.constraints = kept
        report: Dict = {
            "threshold": d,
            "before": len(kept) + dropped,
            "kept": len(kept),
            "dropped": dropped,
            "without_distance_column": missing,
            "self_loops_dropped": self_loop_dropped,
            # 阈值 > 0 且一条都没删：该选项对本次输入完全没有作用
            "ignored": bool(d > 0 and dropped == 0),
        }
        self.diagnostics["distance_filter"] = report
        return report

    def informative_edges(self) -> Dict[str, float]:
        """聚合为 ``edge`` 字典（键 ``"donor,receptor"`` -> 权重和）。"""
        edge: Dict[str, float] = {}
        for c in self.constraints:
            key = c.to_edge_key()
            if key not in edge:
                edge[key] = 0.0
            edge[key] += c.weight
        return edge
