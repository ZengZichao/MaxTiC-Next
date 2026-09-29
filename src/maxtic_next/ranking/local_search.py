"""局部搜索（Metropolis 采样）。

等价原版 ``MaxTiC.optimisation_locale(order, duration)``：以当前排序为起点，在固定时长内
反复做"区间旋转交换"（取 i≠j，将位置 a..b 旋转），用 Metropolis 准则接受/拒绝；
维护历史最优 ``best_order``。

扩展：
* **增量 scoring**（``incremental=True``）：用 ``IncrementalValueComputer`` 在每次
  区间旋转后仅计算增量，避免全量遍历所有边（O(|E|) → O(b-a)）。每 1000 次提交
  自动校正浮点漂移，不消耗随机数。
* **检查点/续跑**（``checkpoint_path`` 不为 None）：定期保存搜索状态（排序、最优解、
  RNG 状态），中断后可从检查点续跑。检查点保存/加载不消耗随机数，不影响复现性。
* **协作式取消**（``stop_check=callable``）：Studio 的「取消」按钮经
  ``maxtic_studio.cancellation.StopToken`` → ``api.rank(stop_check=...)`` → ``Ranker.run`` →
  本函数，在**迭代边界**真正打断搜索（返回截至当下的历史最优，并在 ``stats_out``
  记 ``cancelled=True``）。默认 ``None`` 时与原版行为完全一致。
* **可复现性边界**：停止条件是**墙钟时长**，故 ``--local-search > 0`` 时
  即使固定 ``--seed``，迭代次数也依机器负载而变，**不保证跨机器逐字节复现**；
  只有默认路径（``--ls 0``）与显式设定 ``max_iterations`` 的路径是确定的。
  ``Ranker.run()`` 在开启局部搜索时会打印同义警告，并把实际迭代次数记入
  ``Result.run_metadata``。

> 注：本模块仅实现**核心搜索算法**（与原版逐字节等价的那部分）。
> 基于局部搜索访问解的**稳健性/敏感性摘要**由 ``NearOptimalCollector`` 负责，
> 本模块不产出"置信区间/后验概率"等用语（术语铁律）。``Ranker`` 已为此预留钩子。
"""

import math
import time
from typing import Callable, Dict, List, Optional

from maxtic_next.config import MAX_NUMBER
from maxtic_next.random_ import RandomWrapper, sample_two_distinct
from maxtic_next.ranking.checkpoint import CheckpointManager
from maxtic_next.ranking.value import IncrementalValueComputer, ValueComputer, value
from maxtic_next.robustness.sensitivity import NearOptimalCollector


def optimisation_locale(
    order: List[str],
    edge: Dict[str, float],
    edge_keys: List[str],
    rng: RandomWrapper,
    temperature: float,
    duration: float,
    order_input: Optional[List[str]] = None,
    similarity: Optional[Callable[[List[str], List[str]], float]] = None,
    collector: Optional[NearOptimalCollector] = None,
    value_computer: Optional[ValueComputer] = None,
    incremental: bool = False,
    checkpoint_path: Optional[str] = None,
    checkpoint_interval: float = 60.0,
    max_iterations: int = 0,
    stats_out: Optional[Dict] = None,
    stop_check: Optional[Callable[[], bool]] = None,
) -> List[str]:
    """等价原版 ``optimisation_locale``，返回搜索得到的最佳排序。

    Args:
        order: 起始排序。
        edge: 边权重字典。
        edge_keys: 参与计算的边键列表。
        rng: 随机包装器（受 ``--seed`` 驱动）。
        temperature: Metropolis 温度。
        duration: 搜索时长（秒）。
        order_input: 输入树给出的排序（用于采样相似度收集，供稳健性摘要使用）。
        similarity: 相似度函数（``similarity(a, b) -> float``），可选。保留参数以
            维持与既有调用方/测试的兼容；当前实现不再于热路径上计算相似度。
        collector: 近优解收集器（``NearOptimalCollector``），可选。传入时在每个被
            接受的排序上记录，用于生成稳健性/敏感性摘要。该调用不消耗
            随机数，不影响原版搜索的复现性。
        value_computer: 价值计算器（``ValueComputer``），可选。若提供且 ``incremental=False``，
            局部搜索使用其 ``value()`` 计算目标值（与全量 ``value`` 严格等价，零漂移）。
            默认 ``None`` 时使用模块级 ``value()``，行为与原版完全一致。
        incremental: 是否使用增量 scoring（``IncrementalValueComputer``）。启用后每次
            区间旋转仅计算增量 O(b-a)，远优于全量 O(|E|)。默认 ``False``（与原版一致）。
        checkpoint_path: 检查点文件路径；``None``（默认）表示禁用检查点。
        checkpoint_interval: 检查点自动保存间隔（秒），默认 60 秒。
        max_iterations: 迭代次数上限（``0`` = 不限，仅受 ``duration`` 约束）。
            **可复现性钩子**：``duration`` 是墙钟时长，迭代次数依机器
            负载而定，故 ``--ls > 0`` 时固定 ``--seed`` **不**保证跨机器复现；
            传入正的 ``max_iterations`` 可得到"时长或迭代数先到为准"的确定路径。
        stats_out: 可选字典，回填实际迭代次数 / 耗时 / 最优值，供上层记入元数据。
        stop_check: **协作式取消回调**（``() -> bool``）。它作为搜索主循环
            体的**第一条语句**被调用：返回真值时循环立即结束，函数**照常返回当前的
            历史最优排序**（不是异常，也不丢弃已找到的解），并在 ``stats_out`` 里记
            ``"cancelled": True`` 与实际迭代次数。默认 ``None``（不检查），故不传该
            参数时行为与原来**逐字节一致**。调用方约定见
            ``maxtic_studio.cancellation``（独立仓库，``StopToken.should_stop``）。

    Returns:
        搜索得到的历史最优排序（被取消时为"截至取消点"的历史最优）。
    """
    # temperature 必须 > 0（与 MCMCSampler.__init__ 一致）
    if temperature <= 0:
        raise ValueError("temperature 必须 > 0；如需纯贪心请把 local_search 设为 0。")
    if duration < 0:
        raise ValueError(f"duration（局部搜索时长）必须 >= 0，收到 {duration}。")
    if max_iterations < 0:
        raise ValueError(f"max_iterations 必须 >= 0，收到 {max_iterations}。")

    # ---- 检查点管理器 ----
    ckpt_mgr = CheckpointManager(checkpoint_path, checkpoint_interval)
    resumed_elapsed = 0.0
    resumed_iterations = 0
    loaded_ckpt = ckpt_mgr.load() if ckpt_mgr.enabled else None

    # ---- 增量计算器 ----
    inc: Optional[IncrementalValueComputer] = None
    if incremental:
        inc = IncrementalValueComputer(edge, edge_keys, order)

    # ---- 初始化搜索状态 ----
    if loaded_ckpt is not None:
        # 从检查点恢复
        order = list(loaded_ckpt.current_order)
        best_order = list(loaded_ckpt.best_order)
        current = loaded_ckpt.current_value
        # ：缺键 / None 时不可直接作为浮点使用（旧写法
        # ``best = state.get("best_value")`` 会把 None 带进后续浮点比较）。
        best = loaded_ckpt.best_value
        if best is None:
            best = (
                value(best_order, edge, edge_keys)
                if inc is None
                else IncrementalValueComputer(edge, edge_keys, best_order).current_value
            )
        resumed_elapsed = loaded_ckpt.elapsed or 0.0
        resumed_iterations = getattr(loaded_ckpt, "iteration", 0) or 0
        rng.setstate(loaded_ckpt.rng_state)
        # 恢复增量计算器（若启用）
        if inc is not None:
            inc = IncrementalValueComputer(edge, edge_keys, order)
            current = inc.current_value
        # 恢复收集器状态
        if collector is not None and loaded_ckpt.collector_state is not None:
            collector._visited = loaded_ckpt.collector_state.get("_visited", {})
            saved_best = loaded_ckpt.collector_state.get("best_value")
            if saved_best is not None:
                collector.best_value = saved_best
    else:
        current = value(order, edge, edge_keys) if inc is None else inc.current_value
        best = current
        best_order = list(order)

    # 记录初始排序（仅暂存，不消耗随机数，零漂移）
    if collector is not None and loaded_ckpt is None:
        collector.add(list(order), current)

    iterations = resumed_iterations
    cancelled = False
    time0 = time.time()
    while time.time() - time0 + resumed_elapsed < duration:
        # 协作式取消：停止检查是循环体的**第一条**语句 —— 取消请求到
        # 真正停下的延迟不超过"已完成的最后一次迭代"，不会再多做一次邻域采样。
        # 返回**已找到的历史最优**而非抛异常：取消不应毁掉已算出的部分结果。
        if stop_check is not None and stop_check():
            cancelled = True
            break
        if max_iterations and iterations >= max_iterations:
            break
        iterations += 1
        pair = sample_two_distinct(len(order), rng)
        if pair is None:
            break  # n < 2，无法交换
        i, j = pair
        a = min(i, j)
        b = max(i, j)

        # ---- 目标值计算 ----
        if inc is not None:
            # 增量 scoring：仅计算增量
            v = inc.trial_value(a, b)
        else:
            # 全量计算
            essai = list(order)
            temp = essai[a]
            essai[a:b] = essai[a + 1 : b + 1]
            essai[b] = temp
            if value_computer is not None:
                v = value_computer.value(essai)
            else:
                v = value(essai, edge, edge_keys)

        if v < MAX_NUMBER:
            if v <= current:
                metropolis_ratio: float = 1
            else:
                metropolis_ratio = math.exp((current - v) / temperature)
            coin = rng.random()
            if coin < metropolis_ratio:
                # 接受移动
                if inc is not None:
                    # ：使用 commit_move 的返回值（全量校正后的真值），
                    # 不再沿用校正前的 trial_value，避免增量状态与 current 脱同步。
                    current = inc.commit_move(a, b)
                    order = inc.order
                else:
                    order = essai
                    current = v
                # 记录被接受的排序（零漂移：仅暂存，不消耗随机数）
                if collector is not None:
                    collector.add(list(order), current)
                # 历史最优更新。注：此判断放在"接受"分支之内（与原版放在循环体
                # 之内、接受分支之外不同），二者在标准 Metropolis 准则下等价：
                # v <= current 必被接受且 best <= current，故被拒绝的移动不可能
                # 满足 v < best（旧注释把这一点说成"原版隐患"并不成立）。
                if v < best:
                    best = v
                    best_order = list(order)

        # ---- 检查点保存 ----
        elapsed = time.time() - time0 + resumed_elapsed
        if ckpt_mgr.should_save(elapsed):
            collector_state = None
            if collector is not None:
                collector_state = {
                    "_visited": collector._visited,
                    "best_value": collector.best_value,
                }
            ckpt_mgr.save(
                current_order=order,
                best_order=best_order,
                best_value=best,
                current_value=current,
                elapsed=elapsed,
                rng=rng,
                collector_state=collector_state,
                iteration=iterations,
            )

    # ---- 清理检查点 ----
    ckpt_mgr.cleanup()
    if stats_out is not None:
        stats_out["iterations"] = iterations
        stats_out["elapsed_sec"] = time.time() - time0 + resumed_elapsed
        stats_out["best_value"] = best
        stats_out["max_iterations"] = max_iterations
        # ：取消是否真的在下层生效，必须由上层可见地记下来
        stats_out["cancelled"] = cancelled
        stats_out["stop_requested"] = stop_check is not None
    return best_order
