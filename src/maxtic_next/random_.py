"""统一的随机性驱动封装。

不可妥协约束：**所有随机性都由单个 ``random.Random(seed)``
实例驱动**，禁止直接使用全局 ``random`` 模块。默认种子为 ``42``，保证无参运行可复现。

原版 ``MaxTiC.py`` 在 ``mix()`` 平局与局部搜索中直接调用全局 ``random.random()``；
此处统一改由 ``RandomWrapper`` 驱动，从而可在固定 ``--seed`` 下复现确定性路径。
"""

import random as _random
from typing import Optional, Tuple


class RandomWrapper:
    """对 ``random.Random(seed)`` 的轻量封装，统一驱动所有随机性。

    方法命名刻意贴近原版调用习惯（``rng.random()``、``rng.randint(n)`` 等），
    其中 ``randint(n)`` 等价于原版 ``int(random.random() * n)``，保证与原版取值一致。
    """

    def __init__(self, seed: int = 42) -> None:
        """初始化随机包装器。

        Args:
            seed: 随机种子，默认 42。
        """
        self.seed: int = seed
        self._rng: _random.Random = _random.Random(seed)

    def random(self) -> float:
        """返回 [0, 1) 区间的随机浮点数，等价于原版 ``random.random()``。"""
        return self._rng.random()

    def randint(self, n: int) -> int:
        """返回 ``[0, n)`` 区间的随机整数，等价于原版 ``int(random.random() * n)``。"""
        return int(self._rng.random() * n)

    def choice(self, seq):
        """从非空序列中随机选取一项，等价于 ``seq[int(random.random() * len(seq))]``。"""
        return seq[int(self._rng.random() * len(seq))]

    def shuffle(self, seq) -> None:
        """原地打乱序列，使用底层 ``Random`` 实例。"""
        self._rng.shuffle(seq)

    def getstate(self):
        """返回底层随机数生成器的内部状态（用于调试/复现）。"""
        return self._rng.getstate()

    def setstate(self, state) -> None:
        """恢复底层随机数生成器的内部状态。"""
        self._rng.setstate(state)


def sample_two_distinct(n: int, rng: RandomWrapper) -> Optional[Tuple[int, int]]:
    """从 ``[0, n)`` 中均匀随机选取两个**不同**的整数，返回 ``(i, j)`` 且 ``i != j``。

    当 ``n < 2`` 时无法选出两个不同整数，返回 ``None``，避免调用方陷入死循环。

    本函数是项目中所有"取两个不同随机元素"场景的统一安全入口：
    ``while j == i: j = rng.randint(n)`` 这类写法在 ``n == 1`` 时恒为死循环。

    Args:
        n: 区间上界（不含），需 ``n >= 1``。
        rng: 统一随机性驱动。

    Returns:
        ``(i, j)`` 元组（``i != j``），或 ``n < 2`` 时返回 ``None``。
    """
    if n < 2:
        return None
    i = rng.randint(n)
    j = rng.randint(n)
    while j == i:
        j = rng.randint(n)
    return i, j
