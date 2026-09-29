"""主类 ``Ranker``：编排"校验 → 过滤 → 三启发式 → 选优 → 局部搜索 → 输出三文件 + stdout"。

本模块的主体流程移植自原版 ``MaxTiC.py`` 主流程（解析、建边、分类、三启发式、选优、
输出），并保持了原版的演化顺序与 ``dict.keys()`` 视图处理习惯。

**相对原版"逐字节等价"的已记录偏离**（；等价性测试须按此偏离比对）：

* ``Ranker.run()`` 增加**无条件前置守卫**：物种树必须是
  二叉树（本项目**不支持**非二叉物种树）、内部标签唯一、叶名唯一、约束端点存在、
  约束权重有限且非负、``--random-type ∈ {0,1,2}``、``0 <= --ts <= 1`` 等；
* ``run()`` **幂等**：内部在物种树副本与约束集副本上工作，不改写调用方
  的 ``Tree`` / ``ConstraintSet``，同一实例两次 ``run()`` 输出逐字节相同；
* ``uninformative`` 百分比分母为 ``total_weight``（；原版把无信息权重双重计入）；
* 偏序产物的哨兵边判定用 ``MAX_NUMBER``（；原版写死 ``100000``）；
* 信息性文件在"删除 0 权重边"**之后**写出（；与原版执行顺序一致）；
* ``--threshold-constraints`` 的统计口径在过滤**之后**，且哨兵边不可删除
  ；
* 局部搜索改进后回写 ``Result.values["local_search"]`` / ``["best"]`` 与 ``best_source``
  并打印原版两行；
* ``--random-trees`` 的 p 值用 ``(k+1)/(n+1)`` 校正，且检验对象为**最终交付序**
  。
* ``run(top_k=...)``：**近优解收集容量**不再写死 50 —— API ``api.rank(top_k=)`` 与
  CLI ``--near-optimal-top-k`` 可覆盖，生效值记入 ``run_metadata["params"]`` 与
  ``run_metadata["near_optimal_top_k"]``。
* ``run(stop_check=...)``：**协作式取消**透传给局部搜索与 MCMC 链，在迭代边界停止并
  交付"截至取消点"的最优序，是否真的打断由 ``run_metadata["cancelled"]`` 回传
  （；默认 ``None`` 时所有输出与基线逐字节一致）。

**行为等价的实现层改动**（不改变任何输出字节）：树上遍历 ``Tree.get_leaves`` /
``Tree._write`` / ``opt`` / ``maximum_distance`` / ``random_order`` 由递归改为**显式栈
迭代**，使 1000+ 内部节点的极端不平衡（caterpillar）物种树不再
``RecursionError``；正常树上的产出与随机数消耗序列与递归版逐位一致。
"""

import sys
import time
from dataclasses import dataclass, field
from functools import cmp_to_key
import math
from typing import Callable, Dict, List, Optional

from maxtic_next.config import (
    DEFAULT_MCMC_BURN_IN_FRACTION,
    DEFAULT_MCMC_ITERS,
    DEFAULT_MCMC_THIN_DIVISOR,
    DEFAULT_NEAR_OPTIMAL_TOP_K,
    DEFAULT_OUTPUT_STYLE,
    DEFAULT_SEED,
    MCMC_STATUS_NOTE,
    MCMC_TEMPERATURE_AUTO,
    MCMC_TEMPERATURE_DIVISOR,
    MAX_NUMBER,
    OUTPUT_SUFFIXES,
    PARTIAL_ORDER_SENTINEL_EXCLUSIVE,
    RANDOM_TYPE_CHOICES,
    THRESHOLD_CONSTRAINTS_MAX,
    THRESHOLD_CONSTRAINTS_MIN,
)
from maxtic_next.constraints.constraint import ConstraintSet
from maxtic_next.io.output import format_summary, write_aux_file, write_three_files
from maxtic_next.io.parsing import validate_constraint_set
from maxtic_next.random_ import RandomWrapper
from maxtic_next.ranking.edge import EdgeBuilder
from maxtic_next.ranking.greedy import order_from_graph
from maxtic_next.ranking.local_search import optimisation_locale
from maxtic_next.ranking.mixing import opt
from maxtic_next.ranking.reachability import ReachabilityMatrix
from maxtic_next.ranking.value import ValueComputer, value
from maxtic_next.robustness.sensitivity import NearOptimalCollector
from maxtic_next.tree.tree import Tree


def _binary_children(tree: Tree, node: int, who: str) -> tuple:
    """取二叉内部节点的两个子节点；非二叉时抛出带节点名与处理建议的错误。"""
    children = tree.get_children(node)
    if len(children) != 2:
        raise ValueError(
            f"{who} 只支持二叉物种树：内部节点 {tree.get_bootstrap(node)!r} 有 "
            f"{len(children)} 个子节点。请先解消多歧（polytomy）或对物种树定根/重定根"
            "，使每个内部节点恰有两个子节点。"
        )
    return children[0], children[1]


def _fmt_diag(v, spec: str = ".4f") -> str:
    """诊断量格式化（``None`` -> ``n/a``），供 ``run_metadata["params"]`` 单行摘要用。"""
    return "n/a" if v is None else format(v, spec)


def check_near_optimal_top_k(top_k) -> int:
    """校验并归一化"近优解收集容量"``top_k``。

    ``top_k`` 是稳健性/敏感性摘要能考虑多少个**去重**近优排序的上限：它过去被写死
    为 50，摘要因此永远不可能看过第 51 个序。现由 CLI ``--near-optimal-top-k`` /
    API ``api.rank(top_k=...)`` / ``Ranker.run(top_k=...)`` 三层可覆盖，并在
    ``run_metadata["params"]`` 里回显生效值。

    Args:
        top_k: 期望的保留条数（整数 ``>= 1``；整数值浮点如 ``100.0`` 亦接受）。

    Returns:
        归一化后的 ``int``。

    Raises:
        ValueError: 非整数、非有限，或 ``< 1``（0 会让摘要没有任何样本）。
    """
    if isinstance(top_k, bool) or not isinstance(top_k, (int, float)):
        raise ValueError(
            f"--near-optimal-top-k / top_k 必须是整数（局部搜索中保留的去重近优排序"
            f"条数），收到 {top_k!r}。"
        )
    if math.isnan(top_k) or math.isinf(top_k):
        raise ValueError(f"--near-optimal-top-k / top_k 必须是有限数值，收到 {top_k!r}。")
    if int(top_k) != top_k:
        raise ValueError(f"--near-optimal-top-k / top_k 必须是整数（保留条数），收到 {top_k!r}。")
    if top_k < 1:
        raise ValueError(
            f"--near-optimal-top-k / top_k 必须 >= 1：0 个近优排序无法构成稳健性/"
            f"敏感性摘要，收到 {top_k}。"
        )
    return int(top_k)


def _mcmc_diagnostics_line(diag: Dict) -> str:
    """把 MCMC 混合诊断压成一行 stdout 文本（诊断必须与统计量同屏）。"""

    def _f(key: str, spec: str = ".3f", default: str = "n/a") -> str:
        v = diag.get(key)
        return default if v is None else format(v, spec)

    coverage = diag.get("state_coverage_fraction")
    if coverage is None:
        cov_s = "状态覆盖=n/a（状态空间过大，未精确枚举）"
    else:
        cov_s = (
            f"状态覆盖={coverage * 100:.1f}%"
            f"（{diag.get('states_visited')}/{diag.get('state_space_size')}）"
        )
    return (
        "[MCMC] 收敛诊断（convergence diagnostics，按抽稀前的记录链计）："
        f"唯一样本占比={_f('unique_sample_fraction')}、"
        f"相邻重复率={_f('adjacent_duplicate_fraction')}、"
        f"最长同态连续段占比={_f('max_chain_run_fraction')}、"
        f"有效样本量 ESS={_f('effective_sample_size', '.1f')}"
        f"（{_f('ess_fraction')} 倍样本量）、{cov_s}；"
        f"能量 sd={_f('energy_sd', '.4f')}"
    )


# ----------------------------------------------------------------------
# 与原版一致的辅助函数（最大距离序、Kendall 相似度、排序树生成）
# ----------------------------------------------------------------------
def _cmp(a, b) -> int:
    return (a > b) - (a < b)


def order_from_tree(tree: Tree) -> List[str]:
    """等价原版 ``order_from_tree(tree)``：返回输入树给出的内部节点排序。"""
    nodes: List[int] = tree.get_nodes()
    root: int = tree.get_root()
    # 先收集**节点 id**（int），排完序再一次性换成 bootstrap 标签（str）。
    # 旧实现复用同一个 list 变量承载两种类型，mypy 会把 `order[i] = <str>` 判为非法赋值。
    ids: List[int] = []
    for n in nodes:
        if not tree.is_leaf(n):
            ids.append(n)

    def _by_depth(x: int, y: int) -> int:
        # 用具名函数而不是 lambda：cmp_to_key 的泛型参数需要能被推断为 int
        return _cmp(tree.distance_from(x, root), tree.distance_from(y, root))

    ids.sort(key=cmp_to_key(_by_depth))
    return [tree.get_bootstrap(n) for n in ids]


def tree_from_order(tree: Tree, order: List[str]) -> str:
    """等价原版 ``tree_from_order(order)``：按排序位置改写分支长度并返回 Newick。

    注意：本函数**原地**改写 ``tree`` 的枝长，故 ``Ranker.run()`` 只在内部副本上
    调用它。
    """
    nodes = tree.get_nodes()
    for n in nodes:
        if not tree.is_root(n):
            if tree.is_leaf(n):
                index = len(order)
            else:
                index = order.index(tree.get_bootstrap(n))
            index_parent = order.index(tree.get_bootstrap(tree.get_parent(n)))
            tree.set_length(n, index - index_parent)
    return tree.write_newick(False)


def maximum_distance(tree: Tree, root: int) -> List[List[str]]:
    """等价原版 ``maximum_distance(tree, root)``：返回两个极端拓扑序（要求二叉树）。

    实现为**显式栈后序遍历**：递归版深度等于树高，在 1200 个内部节点的
    caterpillar 树上 ``RecursionError``；迭代版在正常树上产出与递归版逐元素相同的序。
    """
    results: Dict[int, List[List[str]]] = {}
    stack: List[tuple] = [(root, False)]
    while stack:
        node, expanded = stack.pop()
        c1, c2 = _binary_children(tree, node, "maximum_distance")
        name = tree.get_bootstrap(node)
        if not expanded:
            if tree.is_leaf(c1) and tree.is_leaf(c2):
                results[node] = [[name], [name]]
            elif tree.is_leaf(c1) or tree.is_leaf(c2):
                stack.append((node, True))
                stack.append((c2 if tree.is_leaf(c1) else c1, False))
            else:
                stack.append((node, True))
                stack.append((c2, False))
                stack.append((c1, False))
            continue
        if tree.is_leaf(c1) and tree.is_leaf(c2):
            results[node] = [[name], [name]]
        elif tree.is_leaf(c1):
            oo = results[c2]
            results[node] = [[name] + oo[0], [name] + oo[1]]
        elif tree.is_leaf(c2):
            oo = results[c1]
            results[node] = [[name] + oo[0], [name] + oo[1]]
        else:
            orders1 = results[c1]
            orders2 = results[c2]
            results[node] = [[name] + orders1[0] + orders2[0], [name] + orders2[1] + orders1[1]]
    return results[root]


def kendall_distance(a: List[str], b: List[str]) -> float:
    """等价原版 ``kendall_distance(A, B)``。"""
    result = 0.0
    binv = {b[i]: i for i in range(len(b))}
    for i in range(len(a)):
        for j in range(i + 1, len(a)):
            if binv[a[i]] > binv[a[j]]:
                result += 1
    return result


def kendall_similarity(tree: Tree, a: List[str], b: List[str]) -> float:
    """等价原版 ``kendall_similarity(A, B)``。"""
    m = maximum_distance(tree, tree.get_root())
    max_dist = kendall_distance(m[0], m[1])
    if max_dist == 0.0:
        # 退化情形：单内部节点树，两极端拓扑序完全相同（只有一个内部节点），
        # 分母为零；此时任意合法排序与自身完全一致，相似度定义为 1.0（避免除零崩溃）。
        return 1.0
    return (max_dist - kendall_distance(a, b)) / max_dist


def similarity(tree: Tree, a: List[str], b: List[str]) -> float:
    """等价原版 ``similarity(A, B)``（默认 Kendall 相似度）。"""
    return kendall_similarity(tree, a, b)


def random_order(tree: Tree, root: int, rng: RandomWrapper) -> List[str]:
    """等价原版 ``random_order(tree, root)``：随机生成一个内部节点排序。

    用于 ``--random-trees`` 统计。以 ``rng.random()`` 驱动，与原版
    ``random.random()`` 消耗同一随机源（受 ``--seed`` 驱动）。在每个内部节点处随机
    交织其两棵子树的中序。

    实现为**显式栈后序遍历**（递归版深度等于树高，极端不平衡树爆栈）。
    子树处理顺序为"先左后右、再在本节点交织"，与递归版一致，故同一 ``rng`` 状态下
    消耗的随机数序列与产出排序都与递归版逐位相同。
    """
    results: Dict[int, List[str]] = {}
    stack: List[tuple] = [(root, False)]
    while stack:
        node, expanded = stack.pop()
        c1, c2 = _binary_children(tree, node, "random_order")
        if not expanded:
            if tree.is_leaf(c1) and tree.is_leaf(c2):
                results[node] = [tree.get_bootstrap(node)]
            elif tree.is_leaf(c1) or tree.is_leaf(c2):
                stack.append((node, True))
                stack.append((c2 if tree.is_leaf(c1) else c1, False))
            else:
                stack.append((node, True))
                stack.append((c2, False))  # c1 先出栈 => 与递归版同序
                stack.append((c1, False))
            continue
        name = tree.get_bootstrap(node)
        if tree.is_leaf(c1) and tree.is_leaf(c2):
            results[node] = [name]
        elif tree.is_leaf(c1):
            results[node] = [name] + results[c2]
        elif tree.is_leaf(c2):
            results[node] = [name] + results[c1]
        else:
            results[node] = [name] + _interleave(results[c1], results[c2], rng)
    return results[root]


def _interleave(order1: List[str], order2: List[str], rng: RandomWrapper) -> List[str]:
    """把 ``order2`` 随机交织进 ``order1``（原版 ``random_order`` 的双子树分支）。

    先为 ``order2`` 的每个元素抽一个不重复位置（每个位置一次 ``rng.random()``），
    再按位置依次放入 ``order2`` / ``order1`` 的元素。
    """
    pos = list(range(len(order1) + len(order2)))
    for _ in range(len(order2)):
        index = int(rng.random() * len(pos))
        del pos[index]
    order: List[str] = []
    previous = 0
    for i in pos:
        for _ in range(i - previous):
            order.append(order2[0])
            del order2[0]
        order.append(order1[0])
        del order1[0]
        previous = i + 1
    return order + order2


# ----------------------------------------------------------------------
# 结果聚合对象
# ----------------------------------------------------------------------
@dataclass
class Result:
    """``Ranker.run()`` 的聚合结果，供输出层与 Python API 共用。"""

    constraint_file: str = ""
    input_order: List[str] = field(default_factory=list)
    greedy_order: List[str] = field(default_factory=list)
    mixing_order: List[str] = field(default_factory=list)
    best_order: List[str] = field(default_factory=list)
    best_source: str = ""
    values: Dict[str, float] = field(default_factory=dict)
    uninformative: Dict[str, float] = field(default_factory=dict)
    from_leaf: float = 0.0
    trivial_conflict: float = 0.0
    total_weight: float = 0.0
    internal_node_count: int = 0
    ranked_newick: str = ""
    similarity_to_input: float = 0.0
    informative_lines: List[str] = field(default_factory=list)
    conflicting_lines: List[str] = field(default_factory=list)
    partial_lines: List[str] = field(default_factory=list)
    conflict_with_input: float = 0.0
    partial_total: float = 0.0
    informative_file: str = ""
    conflicting_file: str = ""
    partial_order_file: str = ""
    run_metadata: Dict = field(default_factory=dict)
    sensitivity_summary: Optional[Dict] = None
    html_report_file: str = ""
    dry_run_report: str = ""
    mcmc_samples: Optional[List[List[str]]] = None
    random_stats: Optional[Dict] = None
    # ---- 新增字段（均有默认值，旧构造点不受影响）----
    warnings: List[str] = field(default_factory=list)
    uninformative_percent: Optional[float] = None
    informative_count: int = 0
    removed_by_threshold_weight: float = 0.0
    removed_by_threshold_count: int = 0


# ----------------------------------------------------------------------
# 主类
# ----------------------------------------------------------------------
class Ranker:
    """MaxTiC 主类：按原版主流程编排并产出结果。"""

    def __init__(
        self,
        tree: Tree,
        cset: ConstraintSet,
        rng: RandomWrapper,
        min_transfer_distance: float = 0,
        threshold_constraints: float = 0.0,
        temperature: float = 0.001,
        random_type: int = 0,
        time_for_search: float = 0.0,
        random_trees: int = 0,
        constraint_file: str = "constraints",
        output_prefix: Optional[str] = None,
        mcmc: bool = False,
        mcmc_iters: int = DEFAULT_MCMC_ITERS,
        mcmc_temperature: float = MCMC_TEMPERATURE_AUTO,
        mcmc_burn_in: Optional[int] = None,
        mcmc_thin: Optional[int] = None,
        incremental: bool = True,
        checkpoint_path: Optional[str] = None,
        checkpoint_interval: float = 60.0,
        output_style: str = DEFAULT_OUTPUT_STYLE,
        force: bool = False,
        constraint_files_merged: Optional[List[str]] = None,
        local_search_max_iterations: int = 0,
    ) -> None:
        # ---- 参数即校验（打错数字不得静默按"未生效"跑完）----
        self._validate_params(
            min_transfer_distance=min_transfer_distance,
            threshold_constraints=threshold_constraints,
            temperature=temperature,
            random_type=random_type,
            time_for_search=time_for_search,
            random_trees=random_trees,
            mcmc_iters=mcmc_iters,
            mcmc_temperature=mcmc_temperature,
            mcmc_burn_in=mcmc_burn_in,
            mcmc_thin=mcmc_thin,
            local_search_max_iterations=local_search_max_iterations,
        )

        self.tree = tree
        self.cset = cset
        self.rng = rng
        self.min_transfer_distance = min_transfer_distance
        self.threshold_constraints = threshold_constraints
        self.temperature = temperature
        self.random_type = random_type
        self.time_for_search = time_for_search
        self.random_trees = random_trees
        self.constraint_file = constraint_file
        self.output_prefix = output_prefix if output_prefix is not None else constraint_file
        # 输出命名风格："short"（默认，简洁）/ "legacy"（与原版逐字节等价长名）
        self.output_style = output_style
        # MCMC（默认关闭，仅 --mcmc 显式开启；温度 0.0 == "auto"）
        self.mcmc = mcmc
        self.mcmc_iters = mcmc_iters
        self.mcmc_temperature = mcmc_temperature
        self.mcmc_burn_in = mcmc_burn_in
        self.mcmc_thin = mcmc_thin
        # 增量 scoring / 检查点续跑（默认关闭）
        self.incremental = incremental
        self.checkpoint_path = checkpoint_path
        self.checkpoint_interval = checkpoint_interval
        # 是否允许覆盖既有产物（，对应 CLI --force）
        self.force = force
        # 局部搜索迭代上限（0 = 仅受墙钟约束；> 0 使固定 seed 路径确定）
        self.local_search_max_iterations = local_search_max_iterations
        # 本次运行合并了哪些约束文件（产物前缀只体现第一个文件，须显式记录）
        self.constraint_files_merged: List[str] = list(constraint_files_merged or [])
        # 供可选模块（MCMC）复用的边结构（run() 中赋值）
        self.edge: Dict[str, float] = {}
        self.edge_keys: List[str] = []

    # ------------------------------------------------------------------
    # 参数校验
    # ------------------------------------------------------------------
    @staticmethod
    def _validate_params(
        min_transfer_distance: float,
        threshold_constraints: float,
        temperature: float,
        random_type: int,
        time_for_search: float,
        random_trees: int,
        mcmc_iters: int,
        mcmc_temperature: float,
        mcmc_burn_in: Optional[int],
        mcmc_thin: Optional[int],
        local_search_max_iterations: int = 0,
        near_optimal_top_k: int = DEFAULT_NEAR_OPTIMAL_TOP_K,
    ) -> None:
        def _num(name: str, v):
            if not isinstance(v, (int, float)) or isinstance(v, bool):
                raise ValueError(f"{name} 必须是数值，收到 {v!r}。")
            if math.isnan(v) or math.isinf(v):
                raise ValueError(f"{name} 必须是有限数值，收到 {v!r}。")
            return v

        _num("min_transfer_distance", min_transfer_distance)
        _num("threshold_constraints", threshold_constraints)
        _num("temperature", temperature)
        _num("time_for_search", time_for_search)
        if min_transfer_distance < 0:
            raise ValueError(
                f"--min-transfer-distance 必须 >= 0（按 phylogenetic distance 过滤的下界），"
                f"收到 {min_transfer_distance}。"
            )
        if not (THRESHOLD_CONSTRAINTS_MIN <= threshold_constraints <= THRESHOLD_CONSTRAINTS_MAX):
            raise ValueError(
                "--threshold-constraints 是「删除权重占比最低的约束」的比例，"
                f"取值域为 [{THRESHOLD_CONSTRAINTS_MIN}, {THRESHOLD_CONSTRAINTS_MAX}]，"
                f"收到 {threshold_constraints}。"
            )
        if time_for_search < 0:
            raise ValueError(
                f"--local-search 时长必须 >= 0（0 表示关闭），收到 {time_for_search}。"
            )
        if time_for_search > 0 and temperature <= 0:
            raise ValueError(
                f"--temperature 必须 > 0 才能运行局部搜索（收到 {temperature}）；"
                "如不需要随机接受准则，请把 --local-search 设为 0。"
            )
        if random_type not in RANDOM_TYPE_CHOICES:
            raise ValueError(
                f"--random-type 只能是 {list(RANDOM_TYPE_CHOICES)} 之一"
                "（0=原样，1=随机化方向，2=随机化节点），"
                f"收到 {random_type!r}。非法取值会静默退化为 0（不做随机化），"
                "从而把真实数据当对照分析，已拒绝。"
            )
        if random_trees < 0:
            raise ValueError(f"--random-trees 必须 >= 0，收到 {random_trees}。")
        if mcmc_iters is not None and mcmc_iters <= 0:
            raise ValueError(f"--mcmc-iters 必须是正整数（马尔可夫链步数），收到 {mcmc_iters}。")
        if mcmc_temperature < 0:
            raise ValueError(
                f"--mcmc-temperature 必须是正数，或 0 表示 auto（按总权重自适应），"
                f"收到 {mcmc_temperature}。"
            )
        if mcmc_burn_in is not None and mcmc_burn_in < 0:
            raise ValueError(f"--mcmc-burn-in 必须 >= 0，收到 {mcmc_burn_in}。")
        if mcmc_thin is not None and mcmc_thin < 1:
            raise ValueError(f"--mcmc-thin 必须 >= 1，收到 {mcmc_thin}。")
        if local_search_max_iterations < 0:
            raise ValueError(
                f"--local-search-max-iters 必须 >= 0（0 表示只受时长约束），"
                f"收到 {local_search_max_iterations}。"
            )
        # ：近优解收集容量必须 >= 1（同一判据也被 run(top_k=) 使用）
        check_near_optimal_top_k(near_optimal_top_k)

    # ------------------------------------------------------------------
    # 输入守卫：主流程无条件执行，不依赖 --dry-run
    # ------------------------------------------------------------------
    def _assert_binary_tree(self, tree: Tree) -> None:
        """物种树必须是**二叉树**（本项目不支持多歧）。

        ``opt`` / ``maximum_distance`` / ``random_order`` 只取 ``children[0:2]``，
        多歧节点的多余子树会被静默丢弃：残缺排序漏算约束 → 目标值人为变优 →
        在选优中不公平胜出，最终在 ``tree_from_order`` 崩溃或产出不可比数字。
        """
        bad = tree.non_binary_internal_nodes()
        if bad:
            first_label, first_n = bad[0]
            more = ""
            if len(bad) > 1:
                more = f"（另有 {len(bad) - 1} 个节点同样非二叉）"
            raise ValueError(
                f"物种树不是二叉树：内部节点 {first_label!r} 有 {first_n} 个子节点"
                f"{more}。MaxTiC-Next 不支持非二叉物种树：多歧分支会被静默丢弃，"
                "导致目标函数失真与排序崩溃。"
                "请任选其一修复：(1) 解消多歧（对多歧节点做随机/启发式二叉解析）；"
                "(2) 重新定根，使根节点仅有两个子节点（未定根 Newick 的根天然三出）；"
                "(3) 使用支持多歧的下游工具。"
            )

    def _assert_unique_leaf_names(self, tree: Tree) -> None:
        """叶子名必须唯一。

        重名叶会使 ``label_to_id`` 与 ``graph`` 键塌缩，而 ``leaves`` 列表仍含重复项，
        结果是贪婪少排一个内部节点，最终在 ``tree_from_order`` 抛无上下文的
        ``ValueError: list.index(x): x not in list``。
        """
        dup = tree.duplicate_leaf_names()
        if dup:
            detail = "; ".join(
                f"{name!r} 出现在父节点 {parents}" for name, parents in sorted(dup.items())[:5]
            )
            raise ValueError(
                f"物种树存在重复的叶子名（共 {len(dup)} 个名字重复）：{detail}。"
                "重复类群名会让节点标识塌缩，使排序丢失内部节点并崩溃。"
                "请检查是否误把基因树/带样本编号的拷贝当作物种树输入，"
                "并为每个叶子使用唯一名称。"
            )

    def _assert_permutation(self, order: List[str], internal_nodes: List[str], who: str) -> None:
        """断言 ``order`` 是"全部内部节点"的一个排列。"""
        missing = [x for x in internal_nodes if x not in order]
        extra = [x for x in order if x not in internal_nodes]
        if len(set(order)) != len(order) or missing or extra:
            raise ValueError(
                f"{who} 产出的序不是内部节点全集的合法排列："
                f"缺失 {missing[:5]}{'…' if len(missing) > 5 else ''}、"
                f"多余 {extra[:5]}{'…' if len(extra) > 5 else ''}、"
                f"重复项 {len(order) - len(set(order))} 个。"
                "这通常意味着物种树含多歧分支或标签塌缩（重名），请修正输入树后重试。"
            )

    # ------------------------------------------------------------------
    # 主流程：逐行对照原版 MaxTiC.py
    # ------------------------------------------------------------------
    def run(
        self,
        print_summary: bool = True,
        html_report: bool = True,
        seed: int = DEFAULT_SEED,
        top_k: int = DEFAULT_NEAR_OPTIMAL_TOP_K,
        mcmc: Optional[bool] = None,
        mcmc_iters: Optional[int] = None,
        mcmc_temperature: Optional[float] = None,
        mcmc_burn_in: Optional[int] = None,
        mcmc_thin: Optional[int] = None,
        stop_check: Optional[Callable[[], bool]] = None,
    ) -> Result:
        """执行完整排序流程。

        Args:
            print_summary: 是否打印与原版逐字节对齐的 stdout 摘要。
            html_report: 是否生成交互式 HTML 报告。
            seed: 重新播种随机源用的种子（幂等运行的另一半）。
            top_k: **近优解收集容量**：局部搜索期间为"稳健性/敏感性
                摘要"保留多少个**去重**近优排序。必须为 ``>= 1`` 的整数，默认
                ``config.DEFAULT_NEAR_OPTIMAL_TOP_K``（50，与历史行为一致）。50 是
                摘要的**硬上限** —— 想考察更多近优序就必须调大它，故 CLI
                ``--near-optimal-top-k`` / API ``api.rank(top_k=...)`` 均可覆盖，
                生效值记入 ``Result.run_metadata["params"]`` 与
                ``run_metadata["near_optimal_top_k"]``。
            mcmc / mcmc_iters / mcmc_temperature / mcmc_burn_in / mcmc_thin:
                本次运行覆盖构造期的 MCMC 设置（``None`` = 沿用构造值）。
            stop_check: **协作式取消回调**（``() -> bool``）：传给局部搜索
                （``optimisation_locale``）与 MCMC 链，在**迭代边界**打断计算并返回
                截至当下的结果；不抛异常、不毁掉已算出的解。默认 ``None`` 时行为与
                原来逐字节一致（stdout 基线不受影响）。取消是否真的发生可从
                ``Result.run_metadata["cancelled"]`` 读出。

        ：树上遍历（``Tree.get_leaves`` / ``_write`` / ``opt`` /
        ``maximum_distance`` / ``random_order``）全部为**显式栈迭代**实现，故深度
        不平衡（caterpillar）的物种树不再受 Python 递归上限约束；本函数因此**不再**
        临时抬高 ``sys.setrecursionlimit``（抬高上限只会把"必然的爆栈"变成更难诊断的
        C 栈溢出）。
        """
        return self._run_impl(
            print_summary=print_summary,
            html_report=html_report,
            seed=seed,
            top_k=check_near_optimal_top_k(top_k),
            mcmc=mcmc,
            mcmc_iters=mcmc_iters,
            mcmc_temperature=mcmc_temperature,
            mcmc_burn_in=mcmc_burn_in,
            mcmc_thin=mcmc_thin,
            stop_check=stop_check,
        )

    def _run_impl(
        self,
        print_summary: bool = True,
        html_report: bool = True,
        seed: int = DEFAULT_SEED,
        top_k: int = DEFAULT_NEAR_OPTIMAL_TOP_K,
        mcmc: Optional[bool] = None,
        mcmc_iters: Optional[int] = None,
        mcmc_temperature: Optional[float] = None,
        mcmc_burn_in: Optional[int] = None,
        mcmc_thin: Optional[int] = None,
        stop_check: Optional[Callable[[], bool]] = None,
    ) -> Result:
        warnings: List[str] = []
        # 直接调用 ``_run_impl`` 的路径同样要卡住容量（run() 已归一化，这里是幂等重申）
        top_k = check_near_optimal_top_k(top_k)

        # （幂等的另一半）：每次运行都从 ``seed`` **重新播种**随机源。
        # 否则第一次 run() 消耗掉的随机数会推进链状态，第二次 run() 在 mix 平局 /
        # 局部搜索 / MCMC 上走出不同路径，"同实例重复运行"就不是逐字节复现了。
        # ``RandomWrapper`` 不提供 reseed 接口，因此整对象替换为等价的新实例。
        self.rng = RandomWrapper(seed)

        # ---- 0. 输入守卫（无条件执行）----
        validate_constraint_set(self.cset, origin=self.constraint_file)
        self._assert_binary_tree(self.tree)
        self._assert_unique_leaf_names(self.tree)

        # ：run() 必须幂等 —— 一切原地改写（枝长、filter_by_distance）都发生
        # 在**副本**上，调用方的 Tree / ConstraintSet 不受影响。
        # ``diagnostics`` 也要**带上**（浅拷贝，避免把 ``distance_filter`` 段写回调用方
        # 的约束集）：它是上游适配器"丢弃了多少、按什么口径"的唯一机器可读记录。
        tree = self.tree.copy()
        cset = ConstraintSet(list(self.cset.constraints), diagnostics=dict(self.cset.diagnostics))

        nodes = tree.get_nodes()
        leaves = tree.get_leaves_names()
        root = tree.get_root()
        internal_nodes = tree.internal_node_labels()

        # 守卫：内部节点标签唯一性强制校验（主流程，非可选）
        # 重复或空标签会导致 label_to_id 字典键塌缩 → graph 建图自环 → 死循环或静默丢节点
        if len(internal_nodes) != len(set(internal_nodes)):
            dup = [x for x in internal_nodes if internal_nodes.count(x) > 1]
            raise ValueError(
                f"物种树内部节点标签不唯一（重复标签：{sorted(map(str, set(dup)))}），"
                "会导致排序死循环或静默丢失节点。请为每个内部节点设置唯一的 bootstrap 标签。"
            )
        # 空标签也需要检测（内部节点无标注时 bootstrap 为空串）
        if any(not str(lbl).strip() for lbl in internal_nodes):
            raise ValueError(
                "物种树存在内部节点无标签（bootstrap 为空），"
                "算法层无法构建唯一标识。请为所有内部节点设置标签。"
            )

        # 守卫：约束端点合法性主流程校验（非可选）
        # 未校验时，约束文件中的拼写错误会表现为无上下文的 KeyError
        valid_labels = set(str(x) for x in internal_nodes) | set(leaves)
        for c in cset.constraints:
            if c.donor not in valid_labels or c.receptor not in valid_labels:
                raise ValueError(
                    f"约束端点不在物种树中：donor={c.donor!r}, receptor={c.receptor!r}"
                    f"（约束文件/来源：{self.constraint_file}）。"
                    "请检查约束文件中的标签是否与物种树内部节点 bootstrap 或叶子名一致。"
                )

        # 1. 距离过滤（等价于原版解析时的 d 过滤）
        # ：本步骤**如实消费** ``filter_by_distance`` 返回的诊断
        # （kept / dropped / without_distance_column / self_loops_dropped / ignored），
        # 而不是自己再猜。真实行为是"无距离列的约束一律保留"，与旧预检文案
        # （"将丢弃全部约束"）恰好相反；因此 d > 0 却一条没删时必须报
        # **error 级**"该选项被忽略"（与 ``dry_run.py`` 同口径）。
        distance_report = cset.filter_by_distance(self.min_transfer_distance)
        n_before = int(distance_report.get("before", len(cset.constraints)))
        n_kept = int(distance_report.get("kept", len(cset.constraints)))
        n_dropped = int(distance_report.get("dropped", 0))
        n_without_col = int(distance_report.get("without_distance_column", 0))
        n_self_loops = int(distance_report.get("self_loops_dropped", 0))
        if self.min_transfer_distance > 0 and bool(distance_report.get("ignored")):
            if n_without_col and n_without_col == n_before:
                warnings.append(
                    f"ERROR: --min-transfer-distance={self.min_transfer_distance} "
                    f"被忽略（{n_without_col}/{n_before} 条约束不含 phylogenetic "
                    f"distance 列，无距离列的约束一律保留，即本次运行一条也没按距离"
                    f"删除）。与 --dry-run 报告同口径：这是 error 级，因为用户要求了"
                    "过滤、却没有得到过滤。需要按距离过滤请提供 4 列（空格格式）或"
                    "5 列（ALE 逗号格式）输入；不需要请去掉该选项。"
                )
            else:
                warnings.append(
                    f"WARNING: --min-transfer-distance={self.min_transfer_distance} "
                    f"未删除任何约束（{n_kept}/{n_before} 条距离均大于该阈值）。"
                )
        elif n_dropped:
            warnings.append(
                f"NOTE: --min-transfer-distance={self.min_transfer_distance} 丢弃 "
                f"{n_dropped}/{n_before} 条约束（保留 {n_kept} 条）。"
            )
        if n_self_loops:
            # ：donor == receptor 的"供体=受体谱系"事件按距离被吃掉时须计数
            warnings.append(
                f"NOTE: 上述丢弃中有 {n_self_loops} 条是**自环约束**"
                "（donor == receptor，其拓扑距离恒为 0）；它们不会进入"
                '原版口径的 "to itself" 统计。'
            )

        # 2. 聚合约束为 edge 字典
        edge = EdgeBuilder(tree, cset, self.rng, self.random_type).build()
        total_transfers = sum(edge.values())

        # 3. 初始化 graph / degre_entrant（内部节点 + 叶子）
        graph: Dict[str, List[str]] = {}
        degre_entrant: Dict[str, List[str]] = {}
        for n in nodes:
            if not tree.is_leaf(n):
                node = tree.get_bootstrap(n)
                graph[node] = []
                degre_entrant[node] = []
            else:
                node = tree.get_name(n)
                graph[node] = []

        # 4. 分支初始化：把树种系边加入 graph，并将约束中的树边标记为 MAX_NUMBER
        to_desc = 0.0
        uninformative = 0.0
        sentinel_hits: List[tuple] = []  # (边键, 原权重)：与谱系同向的约束
        for n in nodes:
            if not tree.is_leaf(n):
                if not tree.is_root(n):
                    parent = tree.get_bootstrap(tree.get_parent(n))
                    child = tree.get_bootstrap(n)
                    graph[parent].append(child)
                    key = parent + "," + child
                    if key in edge:
                        to_desc += edge[key]
                        uninformative += edge[key]
                        sentinel_hits.append((key, edge[key]))
                    edge[key] = MAX_NUMBER
            else:
                parent = tree.get_bootstrap(tree.get_parent(n))
                child = tree.get_name(n)
                graph[parent].append(child)
        if sentinel_hits:
            # 原版 MaxTiC.py:444 有逐条命中打印（在原版执行路径中不可达）；此处聚合为
            # 一条 NOTE，说明这些约束已被置为不可违反的拓扑硬约束。
            # 措辞刻意用"上方…已计入"：本 NOTE 随 warnings 打在摘要**末尾**，
            # uninformative 统计行在它上方；该权重在 total_weight 分母里只计一次。
            sent_weight = sum(w for _, w in sentinel_hits)
            warnings.append(
                f"NOTE: {len(sentinel_hits)} 条约束与物种树谱系边同向"
                f"（权重合计 {sent_weight}），已被置为不可违反的拓扑硬约束，"
                "并已计入上方 uninformative 统计（在 total_weight 分母中只计一次，"
                "不重复计权）。"
            )

        # 5. 权重阈值过滤（THRESHOLD_CONSTRAINTS）
        # ：树种系边已被置为 MAX_NUMBER，**必须**排除在可删除集合之外，
        # 否则 --ts 足够大时会连拓扑硬约束一起删掉，产出"不再受物种树约束"的排名。
        # ：统计口径在过滤之后计算（分母不含被删除的权重）。
        thresholdable = [k for k in edge.keys() if edge[k] < MAX_NUMBER]
        thresholdable.sort(key=lambda k: edge[k])
        deletable_total = sum(edge[k] for k in thresholdable)
        target = deletable_total * self.threshold_constraints
        sub_total = 0.0
        removed_by_threshold_weight = 0.0
        removed_by_threshold_count = 0
        for key in thresholdable:
            if sub_total >= target:
                break
            sub_total += edge[key]
            removed_by_threshold_weight += edge[key]
            removed_by_threshold_count += 1
            del edge[key]
        total_transfers = total_transfers - removed_by_threshold_weight

        # 6. 构建 degre_entrant（仅保留 < MAX_NUMBER 的边）
        for k in edge.keys():
            if edge[k] < MAX_NUMBER:
                first = k.split(",")[0]
                second = k.split(",")[1]
                if second in degre_entrant:
                    degre_entrant[second].append(first)

        # 6.5 初始化可达性矩阵
        # 此处 graph 仅含树种系边，后续 step 11 贪婪启发式中动态增量更新
        all_labels = list(graph.keys())
        tree_edges_reach = []
        for parent_label, children in graph.items():
            for child_label in children:
                tree_edges_reach.append((parent_label, child_label))
        reach = ReachabilityMatrix(all_labels, tree_edges_reach)

        # 7. 信息性分类
        trivial_conflict = 0.0
        to_itself = 0.0
        to_leaf = 0.0
        to_anc = 0.0
        from_leaf = 0.0
        edge_keys = sorted(edge.keys(), key=lambda k: edge[k], reverse=True)
        for e in edge_keys:
            first = e.split(",")[0]
            second = e.split(",")[1]
            if first == second and edge[e] < MAX_NUMBER:
                uninformative += edge[e]
                to_itself += edge[e]
            elif reach.can_reach(first, second) and edge[e] < MAX_NUMBER:
                uninformative += edge[e]
                to_desc += edge[e]
            elif second in leaves and edge[e] < MAX_NUMBER:
                uninformative += edge[e]
                to_leaf += edge[e]
            elif first in leaves and edge[e] < MAX_NUMBER:
                uninformative += edge[e]
                from_leaf += edge[e]
            elif reach.can_reach(second, first) and edge[e] < MAX_NUMBER:
                uninformative += edge[e]
                trivial_conflict += edge[e]
                to_anc += edge[e]

        # 8. 互反约束减法（累加 trivial_conflict）
        edge_keys = sorted(edge.keys(), key=lambda k: edge[k], reverse=True)
        for e in edge_keys:
            words = e.split(",")
            opposite = words[1] + "," + words[0]
            if (
                edge[e] < MAX_NUMBER
                and opposite in edge
                and edge[opposite] < MAX_NUMBER
                and edge[e] >= edge[opposite]
            ):
                trivial_conflict += edge[opposite]

        # 移除权重为 0 的边
        for k in list(edge.keys()):
            if edge[k] == 0:
                del edge[k]

        # 8.5 输出信息性约束（必须在删除 0 权重边**之后**，与原版顺序一致）
        edge_keys = sorted(edge.keys(), key=lambda k: edge[k], reverse=True)
        informative_lines: List[str] = []
        for k in edge_keys:
            first = k.split(",")[0]
            second = k.split(",")[1]
            if (
                edge[k] < MAX_NUMBER
                and first != second
                and not reach.can_reach(first, second)
                and second not in leaves
                and first not in leaves
                and not reach.can_reach(second, first)
            ):
                informative_lines.append(f"{k} {edge[k]}")

        # 9. 输入树排序 + value
        order_input = order_from_tree(tree)
        value_input = value(order_input, edge, edge_keys)

        # 10. 贪婪启发式
        rejected = 0.0
        for e in edge_keys:
            words = e.split(",")
            if words[0] != words[1] and reach.can_reach(words[1], words[0]):
                rejected += edge[e]
            else:
                graph[words[0]].append(words[1])
                reach.add_edge(words[0], words[1])
        order_greedy = order_from_graph(graph, leaves)
        value_greedy = value(order_greedy, edge, edge_keys)

        # 11. 混合启发式
        order_heuristic = opt(tree, root, degre_entrant, edge, self.rng)
        order_heuristic.reverse()
        value_heuristic = value(order_heuristic, edge, edge_keys)

        # 12. 选优（NaN 安全比较：，这里再显式排序）
        if value_greedy <= value_heuristic:
            order = list(order_greedy)
            pre_search_source = "greedy heuristic"
        else:
            order = list(order_heuristic)
            pre_search_source = "mixing heuristic"

        # 12.5 三个候选序都必须是"内部节点全集的一个排列"
        self._assert_permutation(order_input, internal_nodes, "order_from_tree(输入树)")
        self._assert_permutation(order_greedy, internal_nodes, "order_from_graph(贪婪)")
        self._assert_permutation(order_heuristic, internal_nodes, "opt/mix(混合)")

        # 13. 生成排序树 Newick（在树上，枝长被覆写；已在副本上，安全）
        ranked_newick = tree_from_order(tree, order)

        # 14. 与输入树的相似度
        sim = similarity(tree, order_input, order)

        # 14.5 口径闭合告警（不得"空数据得满分"）
        if not cset.constraints:
            warnings.append(
                "WARNING: 输入约束集为空（0 条时间约束）。"
                "本次输出的所有 value 与百分比恒为 0.0，这只表示"
                "「没有可违反的约束」，不是「完美时间一致性」的证据。"
            )
        if self.threshold_constraints > 0 and not informative_lines:
            warnings.append(
                f"WARNING: --threshold-constraints={self.threshold_constraints} 删除了 "
                f"{removed_by_threshold_count} 条约束（权重合计 {removed_by_threshold_weight}）"
                "后，剩余信息性约束 0 条。此时所有 value 必为 0.0、相似度为平凡值，"
                "这不是「完美时间一致性」的证据，而是「无数据可检验」。请调低该阈值。"
            )
        elif total_transfers == 0.0 and cset.constraints:
            warnings.append(
                "WARNING: 过滤后总约束权重为 0.0（无任何信息性约束），"
                "本次输出的所有 value 与百分比全部退化为 0.0，不构成一致性证据。"
            )

        # 15. 局部搜索
        sensitivity_summary = None
        local_search_value: Optional[float] = None
        ls_stats: Dict = {}
        if self.time_for_search > 0:
            if print_summary:
                print(
                    f"attempting a local search from the best found order, please wait "
                    f"{self.time_for_search} seconds"
                )
            value_computer = ValueComputer(edge, edge_keys)
            collector = NearOptimalCollector(top_k=top_k)
            if not self.local_search_max_iterations:
                warnings.append(
                    "WARNING: --local-search > 0 以墙钟时长为停止条件，实际迭代次数依机器"
                    "负载而变，因此固定 --seed **不保证跨机器复现**（默认 --ls 0 路径完全"
                    "可复现）。需要确定性请设置 --local-search-max-iters。"
                )
            order = optimisation_locale(
                order,
                edge,
                edge_keys,
                self.rng,
                self.temperature,
                self.time_for_search,
                order_input=order_input,
                similarity=lambda a, b: similarity(tree, a, b),
                collector=collector,
                value_computer=value_computer,
                incremental=self.incremental,
                checkpoint_path=self.checkpoint_path,
                checkpoint_interval=self.checkpoint_interval,
                max_iterations=self.local_search_max_iterations,
                stats_out=ls_stats,
                stop_check=stop_check,
            )
            if ls_stats.get("cancelled"):
                # 取消不是错误：如实说明交付的是"截至取消点"的解
                warnings.append(
                    "NOTE: 本次运行按请求**取消**：局部搜索在第 "
                    f"{ls_stats.get('iterations')} 次迭代边界停止（计划的搜索时长为 "
                    f"{self.time_for_search} 秒，未跑满）。以下交付的是**截至取消点**"
                    "的历史最优序。"
                )
            # 墙钟停止条件让"同种子不同机"结果不同（README 已如实声明）；
            # 这里把**本次实际迭代数**回灌成可执行的复现配方，
            # 且只写 stderr —— stdout 与原版逐字节对齐是硬契约，不能加行。
            iters = ls_stats.get("iterations")
            elapsed = ls_stats.get("elapsed_sec")
            if (
                self.time_for_search > 0
                and not self.local_search_max_iterations
                and iters is not None
            ):
                hint = (
                    f"NOTE: 本次局部搜索实际执行 {iters} 次迭代"
                    + (f"（{float(elapsed):.2f} 秒）" if elapsed is not None else "")
                    + f"，停止条件是墙钟 {self.time_for_search} 秒。"
                    "跨机器复现请改用 --local-search-max-iters "
                    + str(iters)
                    + "（迭代数一致时结果逐字节一致）。"
                )
                warnings.append(hint)
                print(hint, file=sys.stderr)
            sensitivity_summary = collector.summary()
            # ：回写改进后的 value / best_source，使 stdout 不与交付树矛盾
            local_search_value = value(order, edge, edge_keys)
            ranked_newick = tree_from_order(tree, order)
            sim = similarity(tree, order_input, order)
            best_source = (
                f"{pre_search_source} + local search"
                if local_search_value
                < (value_greedy if pre_search_source == "greedy heuristic" else value_heuristic)
                else pre_search_source
            )
        else:
            best_source = pre_search_source

        self._assert_permutation(order, internal_nodes, f"最终交付序（{best_source}）")

        # 16.5 持久化边结构，供后续可选模块（MCMC）复用
        self.edge = edge
        self.edge_keys = edge_keys

        # ---- 16.6 --random-trees 置换检验 ----
        # p 值用 (k+1)/(n+1) 校正（Phipson & Smyth 2010），且检验对象为
        # **最终交付序**（局部搜索之后），而非搜索前的启发式序。
        random_stats = None
        if self.random_trees > 0:
            # 原版在 --random-trees 时打印进度（MaxTiC.py:652-656）；这里同样给出反馈，
            # 但只写 stderr，以保持 stdout 仍是可机读的摘要。
            sys.stderr.write(
                f"evaluating {self.random_trees} random orders"
                " (permutation test, --random-trees) ...\n"
            )
            sys.stderr.flush()
            time0 = time.time()
            dist_path = self.output_prefix + OUTPUT_SUFFIXES[self.output_style][3]
            vmin: float = MAX_NUMBER
            vmax = 0.0
            k_value = 0
            sim_min: float = MAX_NUMBER
            sim_max = 0.0
            k_sim = 0
            n_rand = self.random_trees
            best_value = value(order, edge, edge_keys)
            best_sim = similarity(tree, order_input, order)
            dist_lines: List[str] = []
            for _ in range(n_rand):
                rorder = random_order(tree, root, self.rng)
                val = value(rorder, edge, edge_keys)
                sval = similarity(tree, rorder, order_input)
                dist_lines.append(f"{val} {sval}")
                vmin = min(val, vmin)
                vmax = max(val, vmax)
                if val <= best_value:
                    k_value += 1
                sim_min = min(sval, sim_min)
                sim_max = max(sval, sim_max)
                if sval >= best_sim:
                    k_sim += 1
            write_aux_file(dist_path, dist_lines, force=self.force)
            random_stats = {
                "n": n_rand,
                "value_min": vmin,
                "value_max": vmax,
                "value_k": k_value,
                "value_pvalue": (k_value + 1) / (n_rand + 1),
                "sim_min": sim_min,
                "sim_max": sim_max,
                "sim_k": k_sim,
                "sim_pvalue": (k_sim + 1) / (n_rand + 1),
                "distribution_file": dist_path,
                "elapsed_sec": time.time() - time0,
                "tested_order": "delivered best order (after local search)",
                "correction": "(k+1)/(n+1)",
            }

        # ---- 线性扩展上的 Metropolis–Hastings（默认关闭）----
        # ：这不是"严格后验采样器"。温度自适应 + burn-in + thin 后仍只是
        # 一条 MH 链；收敛诊断未经验证，故状态说明固定为 MCMC_STATUS_NOTE。
        do_mcmc = self.mcmc if mcmc is None else mcmc
        eff_iters = self.mcmc_iters if mcmc_iters is None else mcmc_iters
        eff_temp = self.mcmc_temperature if mcmc_temperature is None else mcmc_temperature
        eff_burn_in = self.mcmc_burn_in if mcmc_burn_in is None else mcmc_burn_in
        eff_thin = self.mcmc_thin if mcmc_thin is None else mcmc_thin
        mcmc_diagnostics_dict: Optional[Dict] = None
        chain_diag: Optional[Dict] = None
        if do_mcmc:
            if not eff_iters or eff_iters <= 0:
                raise ValueError(f"--mcmc-iters 必须是正整数（马尔可夫链步数），收到 {eff_iters}。")
            temp_is_auto = eff_temp is None or eff_temp <= 0
            eff_temp_auto = temp_is_auto
            if temp_is_auto:
                eff_temp = max(total_transfers, 1.0) / MCMC_TEMPERATURE_DIVISOR
            if eff_burn_in is None:
                eff_burn_in = int(eff_iters * DEFAULT_MCMC_BURN_IN_FRACTION)
            if eff_thin is None:
                eff_thin = max(1, eff_iters // DEFAULT_MCMC_THIN_DIVISOR)
            from maxtic_next.robustness.mcmc import (
                MCMCSampler,
                mcmc_diagnostics as compute_diagnostics,
            )

            sampler = MCMCSampler(
                tree=tree,
                edge=self.edge,
                edge_keys=self.edge_keys,
                rng=self.rng,
                temperature=eff_temp,
                initial_order=order,
            )
            # burn-in：先推进链但丢弃样本
            mcmc_stopped_early = False
            if eff_burn_in > 0:
                sampler.sample(eff_burn_in, burn_in=0, stop_check=stop_check)
            mcmc_samples = sampler.sample(eff_iters, thin=max(1, eff_thin), stop_check=stop_check)
            mcmc_stopped_early = bool(getattr(sampler, "last_stopped_early", False))
            if not mcmc_samples:
                raise RuntimeError(
                    "MCMC 未产出任何样本（iters="
                    f"{eff_iters}, thin={eff_thin}）"
                    + (
                        "：链在产出第一个样本之前就被**取消**（stop_check），"
                        "本次运行没有可交付的后验近似样本。"
                        if mcmc_stopped_early
                        else "。"
                    )
                )
            energies = [sampler.value_computer.value(s) for s in mcmc_samples]
            mean_e = sum(energies) / len(energies)
            min_e = min(energies)
            n_unique = len({tuple(s) for s in mcmc_samples})
            # 混合诊断（诊断必须与统计量一并展示，两处各有其用途）
            #  * chain_diag：**抽稀前**的整条记录链 —— 相邻重复率 / ESS / 状态覆盖
            #    只有在未抽稀的链上才有意义（抽稀本身就是为削弱自相关）；
            #  * mcmc_diagnostics_dict：交付给用户的 ``mcmc_samples``（已抽稀）的
            #    诊断，与 stdout"保留 N 个样本，其中唯一排序 M 个"同一口径。
            chain_diag = sampler.diagnostics()
            mcmc_diagnostics_dict = compute_diagnostics(mcmc_samples, energies)
            for key in (
                "temperature",
                "burn_in",
                "thin",
                "status_note",
                "state_space_size",
                "state_coverage_fraction",
                "states_visited",
            ):
                if key in chain_diag:
                    mcmc_diagnostics_dict[key] = chain_diag[key]
            mcmc_diagnostics_dict["n_recorded_samples"] = len(mcmc_samples)
            mcmc_diagnostics_dict["n_unique_recorded"] = n_unique
            mcmc_diagnostics_dict["chain_n_samples"] = chain_diag.get("n_samples")
            mcmc_diagnostics_dict["chain_n_unique_samples"] = chain_diag.get("n_unique_samples")
            print(f"[MCMC] {MCMC_STATUS_NOTE}")
            print(
                "[MCMC] 口径声明：这是一条**初步实现**的 Metropolis–Hastings 链，"
                "收敛诊断未经验证；链上样本彼此自相关（非独立观测），"
                "未通过下列诊断前不得用作不确定性度量"
            )
            print(
                f"[MCMC] 温度 T={eff_temp}"
                + (f"（auto：max(总权重,1)/{MCMC_TEMPERATURE_DIVISOR:g}）" if eff_temp_auto else "")
                + f"；iters={eff_iters}、burn-in={eff_burn_in}、thin={eff_thin}"
            )
            print(
                f"[MCMC] 保留 {len(mcmc_samples)} 个样本，其中唯一排序 {n_unique} 个；"
                f"能量 mean={mean_e:.4f}, min={min_e:.4f}"
                "（唯一样本数远低于样本数说明链未混合，均值不可用作不确定性度量）"
            )
            print(_mcmc_diagnostics_line(chain_diag))
        else:
            mcmc_samples = None
            mcmc_stopped_early = False
            # MCMC 关闭：不虚构"有效温度"（0.0 即 auto，从未应用于任何链）
            eff_temp_auto = eff_temp is None or eff_temp <= 0
            if eff_burn_in is None:
                eff_burn_in = int((eff_iters or 0) * DEFAULT_MCMC_BURN_IN_FRACTION)
            if eff_thin is None:
                eff_thin = max(1, (eff_iters or 0) // DEFAULT_MCMC_THIN_DIVISOR)

        # 17. 与最优序冲突的约束文件
        index = {}
        for i in range(len(order)):
            index[order[i]] = i
        conflicting_lines: List[str] = []
        for e in edge_keys:
            sommets = e.split(",")
            if (
                sommets[0] in index
                and sommets[1] in index
                and index[sommets[0]] > index[sommets[1]]
            ):
                conflicting_lines.append(f"{e} {edge[e]}")

        # 18. 偏序文件（哨兵判定用 MAX_NUMBER，不再是魔法数 100000）
        conflict_with_input = 0.0
        partial_total = 0.0
        partial_lines: List[str] = []
        for e in edge_keys:
            sommets = e.split(",")
            a, b = sommets[0], sommets[1]
            if (
                a in index
                and b in index
                and index[a] < index[b]
                and edge[e] < PARTIAL_ORDER_SENTINEL_EXCLUSIVE
            ):
                partial_total += edge[e]
                if order_input.index(a) < order_input.index(b):
                    partial_lines.append(f"{a} {b} {edge[e]} black")
                else:
                    partial_lines.append(f"{a} {b} {edge[e]} green")
                    conflict_with_input += edge[e]

        # 19. 写入三文件
        informative_file, conflicting_file, partial_order_file = write_three_files(
            self.output_prefix,
            informative_lines,
            conflicting_lines,
            partial_lines,
            output_style=self.output_style,
            force=self.force,
        )

        # 20. 组装 Result
        values: Dict[str, float] = {
            "input": value_input,
            "greedy": value_greedy,
            "mixing": value_heuristic,
        }
        if local_search_value is not None:
            values["local_search"] = local_search_value
        values["best"] = value(order, edge, edge_keys)
        params = (
            f"--ls {self.time_for_search} --t {self.temperature} "
            f"--r {self.random_type} --d {self.min_transfer_distance} "
            f"--ts {self.threshold_constraints} "
            f"--near-optimal-top-k {top_k} "
            f"--mcmc {int(bool(do_mcmc))} --mcmc-iters {eff_iters} "
            f"--mcmc-temperature {eff_temp}"
            f"{' (auto)' if eff_temp_auto else ''} "
            f"--mcmc-burn-in {eff_burn_in} --mcmc-thin {eff_thin}"
        )
        if mcmc_diagnostics_dict is not None:
            # ：采样器诊断必须与 --mcmc-* 参数一并落到 params 里（可追溯）。
            # 混合类量取自**抽稀前的记录链**（抽稀本身就是为了削弱自相关）。
            # `chain` 只是**读数用的局部别名**：`chain_diag` 本身仍按原样（可能为 None）
            # 落进 run_metadata，避免把 None 归一成 {} 而改变产物元数据摘要。
            chain = chain_diag or {}
            params += (
                " --mcmc-diagnostics"
                f" recorded={mcmc_diagnostics_dict.get('n_samples')}"
                f" unique_recorded={mcmc_diagnostics_dict.get('n_unique_samples')}"
                f" chain_unique={chain.get('n_unique_samples')}"
                f"/{chain.get('n_samples')}"
                f" adjacent_dup={_fmt_diag(chain.get('adjacent_duplicate_fraction'))}"
                f" ess={_fmt_diag(chain.get('effective_sample_size'), '.1f')}"
                f" state_cov={_fmt_diag(chain.get('state_coverage_fraction'))}"
                f" status={MCMC_STATUS_NOTE}"
            )
        run_metadata = {
            "seed": seed,
            "params": params,
            "temperature": self.temperature,
            "random_type": self.random_type,
            "min_transfer_distance": self.min_transfer_distance,
            "threshold_constraints": self.threshold_constraints,
            "time_for_search": self.time_for_search,
            "mcmc": bool(do_mcmc),
            "mcmc_iters": eff_iters,
            "mcmc_temperature": eff_temp if do_mcmc else None,
            "mcmc_temperature_is_auto": bool(eff_temp_auto),
            "mcmc_burn_in": eff_burn_in,
            "mcmc_thin": eff_thin,
            "mcmc_diagnostics": mcmc_diagnostics_dict,
            "mcmc_chain_diagnostics": chain_diag,
            "mcmc_status": MCMC_STATUS_NOTE if do_mcmc else "disabled",
            "mcmc_stopped_early": bool(mcmc_stopped_early),
            # ：近优解收集容量的**生效值**（stdout 之外的机器可读入口）
            "near_optimal_top_k": int(top_k),
            "local_search_max_iterations": self.local_search_max_iterations,
            "local_search_iterations": ls_stats.get("iterations"),
            "local_search_elapsed_sec": ls_stats.get("elapsed_sec"),
            # 机器可读地告诉下游"这次运行是否确定性的"
            "local_search_deterministic": bool(self.local_search_max_iterations)
            or self.time_for_search <= 0,
            # ：协作式取消是否真的在下层执行（stop_check 请求 vs 实际打断）
            "cancelled": bool(ls_stats.get("cancelled")) or bool(mcmc_stopped_early),
            "stop_check_requested": stop_check is not None,
            "removed_by_threshold_weight": removed_by_threshold_weight,
            "removed_by_threshold_count": removed_by_threshold_count,
            "informative_count": len(informative_lines),
            "constraint_files_merged": self.constraint_files_merged,
            # ：上游适配器的机器可读诊断 + 本次距离过滤实况
            "adapter_diagnostics": dict(cset.diagnostics),
            "distance_filter": dict(distance_report),
        }
        uninformative_percent = (
            uninformative * 100 / total_transfers if total_transfers != 0 else 0.0
        )
        result = Result(
            constraint_file=self.constraint_file,
            input_order=order_input,
            greedy_order=order_greedy,
            mixing_order=order_heuristic,
            best_order=order,
            best_source=best_source,
            values=values,
            uninformative={
                "total": uninformative,
                "to_desc": to_desc,
                "to_leaf": to_leaf,
                "to_anc": to_anc,
                "to_itself": to_itself,
            },
            from_leaf=from_leaf,
            trivial_conflict=trivial_conflict,
            total_weight=total_transfers,
            internal_node_count=len(internal_nodes),
            ranked_newick=ranked_newick,
            similarity_to_input=sim,
            informative_lines=informative_lines,
            conflicting_lines=conflicting_lines,
            partial_lines=partial_lines,
            conflict_with_input=conflict_with_input,
            partial_total=partial_total,
            informative_file=informative_file,
            conflicting_file=conflicting_file,
            partial_order_file=partial_order_file,
            run_metadata=run_metadata,
            sensitivity_summary=sensitivity_summary,
            html_report_file="",
            mcmc_samples=mcmc_samples,
            random_stats=random_stats,
            warnings=warnings,
            uninformative_percent=uninformative_percent,
            informative_count=len(informative_lines),
            removed_by_threshold_weight=removed_by_threshold_weight,
            removed_by_threshold_count=removed_by_threshold_count,
        )

        # 20.5 交互式 HTML 报告
        if html_report:
            try:
                from maxtic_next.report.html import write_html_report

                result.html_report_file = write_html_report(
                    self.output_prefix,
                    result,
                    sensitivity_summary=sensitivity_summary,
                    force=self.force,
                )
            except Exception as exc:  # noqa: BLE001 - 报告为可选特性，失败不阻断
                # 降级不阻断主流程，但不能只落在 stdout 上 ——
                # 批处理与图形前端都需要机器可读的入口。
                msg = (
                    f"HTML 报告生成失败，已跳过：{type(exc).__name__}: {exc}"
                    "（三件套与退出码不受影响；这仍是刻意设计的降级路径）"
                )
                warnings.append("WARNING: " + msg)
                result.run_metadata["html_report_error"] = f"{type(exc).__name__}: {exc}"
                print(f"[警告] {msg}")

        if print_summary:
            print(format_summary(result))
        return result
