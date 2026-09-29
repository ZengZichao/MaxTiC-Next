"""稳健性 / 敏感性分析子包。

* ``sensitivity``：基于局部搜索访问解的近优序收集与稳健性/敏感性摘要；
* ``mcmc``：可逆 MCMC 严格采样器。

术语铁律提醒：``sensitivity`` 与 HTML 报告**仅**可称"基于局部搜索访问解的稳健性/敏感性摘要"；
"后验 / 平稳分布 / 细致平衡"等概率分布用语**仅**允许出现在 ``mcmc`` 模块。
"""

from maxtic_next.robustness.mcmc import MCMCSampler
from maxtic_next.robustness.sensitivity import NearOptimalCollector

__all__ = ["NearOptimalCollector", "MCMCSampler"]
