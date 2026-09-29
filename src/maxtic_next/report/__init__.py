"""交互式 HTML 报告模块。

导出：
* ``plot``：图表构建子模块（直方图 / 稳健性摘要 / 排序树可视化）；
* ``html``：报告构建与渲染子模块（Jinja2 优先，字符串/SVG 降级）。

术语铁律：本模块仅使用"稳健性/敏感性摘要"，不出现"置信区间 /
后验概率 / 后验分布"等用语。
"""

from maxtic_next.report import plot
from maxtic_next.report import html

__all__ = ["plot", "html"]
