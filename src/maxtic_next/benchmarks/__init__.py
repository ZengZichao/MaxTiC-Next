"""性能基准子包。

提供合成大规模数据集（大 n 内部节点数 / 大 |E| 约束数）的生成与基准测试，测量 CPU
时间与峰值内存，并标出可处理上限。纯标准库实现，不依赖真实 ALE 或任何第三方包。
"""

from maxtic_next.benchmarks.suite import run_benchmark, run_suite, make_synthetic_tree

__all__ = ["run_benchmark", "run_suite", "make_synthetic_tree"]
