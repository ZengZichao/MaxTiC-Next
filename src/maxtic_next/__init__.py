"""MaxTiC-Next —— MaxTiC 的 Python 3 重写版本。

本包把原版 ``MaxTiC.py``（Eric Tannier, Inria）的核心算法重构为模块化、可测试、
可复现的 Python 3 代码。

**关于"等价"的准确表述**：核心启发式（贪婪 / mixing / 局部搜索 / 目标函数）为忠实
移植，并在固定 ``--seed`` 下与 ``tests/baselines/`` 中**由本版本自产**的快照逐字节
一致（回归基线，防漂移）。相对原版的**已记录偏离**（输入守卫、统计口径、
哨兵阈值、局部搜索后回写、置换检验 p 值校正等）见
``maxtic_next.ranking.ranker`` 模块文档与 ``tests/test_reference_equivalence.py``；
这些偏离使输出**不再**与原版逐字节相同，但每一项都对应原版的一处缺陷。

原版作者：Eric Tannier
引用：MaxTiC: Fast ranking of a phylogenetic tree by Maximum Time Consistency
with lateral gene transfers, Biorxiv doi.org/10.1101/127548
许可证：CeCILL 2.1（继承自原版）
"""

__version__ = "0.1.0"
__author__ = "Zichao Zeng (曾子超, ORCID 0000-0001-6553-970X) — Python 3 rewrite; Eric Tannier — original MaxTiC"
__license__ = "CECILL-2.1"
__original_doi__ = "doi.org/10.1101/127548"

__all__ = ["__version__", "rank", "build_constraints", "Tree", "Ranker"]

from maxtic_next.tree.tree import Tree
from maxtic_next.ranking.ranker import Ranker
from maxtic_next.api import rank, build_constraints
