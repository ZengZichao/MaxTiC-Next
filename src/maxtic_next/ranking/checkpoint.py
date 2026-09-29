"""局部搜索检查点/续跑。

在局部搜索过程中定期保存状态快照，使长时间搜索可在中断后续跑。
检查点为 pickle 文件，保存：

* 当前排序 ``current_order``
* 历史最优排序 ``best_order`` 与最优值 ``best_value``
* 当前目标值 ``current_value``
* 已用搜索时间 ``elapsed``
* RNG 内部状态 ``rng_state``（保证续跑可复现）
* 近优解收集器状态 ``collector_state``（可选）

注意：检查点保存/加载不消耗随机数，不影响搜索路径的复现性。
"""

import os
import pickle
import time
from typing import Any, Dict, List, Optional

from maxtic_next.random_ import RandomWrapper


class LocalSearchCheckpoint:
    """局部搜索检查点的序列化容器。"""

    def __init__(
        self,
        current_order: List[str],
        best_order: List[str],
        best_value: float,
        current_value: float,
        elapsed: float,
        rng_state: Any,
        collector_state: Optional[Dict] = None,
        iteration: int = 0,
    ) -> None:
        self.current_order = list(current_order)
        self.best_order = list(best_order)
        self.best_value = best_value
        self.current_value = current_value
        self.elapsed = elapsed
        self.rng_state = rng_state
        self.collector_state = collector_state
        self.iteration = iteration
        self.saved_at = time.time()


def save_checkpoint(path: str, checkpoint: LocalSearchCheckpoint) -> None:
    """将检查点保存到文件。

    Args:
        path: 检查点文件路径。
        checkpoint: 检查点对象。
    """
    dir_name = os.path.dirname(path)
    if dir_name:
        os.makedirs(dir_name, exist_ok=True)
    # 原子写入：先写临时文件，再重命名
    tmp_path = path + ".tmp"
    with open(tmp_path, "wb") as fh:
        pickle.dump(checkpoint, fh, protocol=pickle.HIGHEST_PROTOCOL)
    os.replace(tmp_path, path)


def load_checkpoint(path: str) -> Optional[LocalSearchCheckpoint]:
    """从文件加载检查点。

    Args:
        path: 检查点文件路径。

    Returns:
        检查点对象；文件不存在或损坏时返回 None。
    """
    if not os.path.exists(path):
        return None
    try:
        with open(path, "rb") as fh:
            return pickle.load(fh)
    except Exception:
        return None


def delete_checkpoint(path: str) -> None:
    """删除检查点文件（搜索完成后清理）。"""
    try:
        if os.path.exists(path):
            os.remove(path)
    except OSError:
        pass


class CheckpointManager:
    """局部搜索检查点管理器：自动定期保存。

    用法：
        mgr = CheckpointManager(path, interval=60)
        # 在搜索循环中
        if mgr.should_save(elapsed):
            mgr.save(...)
        # 搜索完成后
        mgr.cleanup()
    """

    def __init__(self, path: Optional[str], interval: float = 60.0) -> None:
        """初始化检查点管理器。

        Args:
            path: 检查点文件路径；``None`` 表示禁用检查点。
            interval: 自动保存间隔（秒），默认 60 秒。
        """
        self.path = path
        self.interval = interval
        self._last_save = 0.0

    @property
    def enabled(self) -> bool:
        """是否启用检查点。"""
        return self.path is not None

    def should_save(self, elapsed: float) -> bool:
        """判断是否应该保存检查点。"""
        if not self.enabled:
            return False
        return elapsed - self._last_save >= self.interval

    def save(
        self,
        current_order: List[str],
        best_order: List[str],
        best_value: float,
        current_value: float,
        elapsed: float,
        rng: RandomWrapper,
        collector_state: Optional[Dict] = None,
        iteration: int = 0,
    ) -> None:
        """保存检查点。"""
        path = self.path
        if path is None:  # 与 `enabled` 同义，但让类型收窄可见
            return
        ckpt = LocalSearchCheckpoint(
            current_order=current_order,
            best_order=best_order,
            best_value=best_value,
            current_value=current_value,
            elapsed=elapsed,
            rng_state=rng.getstate(),
            collector_state=collector_state,
            iteration=iteration,
        )
        save_checkpoint(path, ckpt)
        self._last_save = elapsed

    def load(self) -> Optional[LocalSearchCheckpoint]:
        """加载检查点。"""
        path = self.path
        if path is None:
            return None
        return load_checkpoint(path)

    def cleanup(self) -> None:
        """搜索完成后删除检查点文件。"""
        path = self.path
        if path is not None:
            delete_checkpoint(path)
