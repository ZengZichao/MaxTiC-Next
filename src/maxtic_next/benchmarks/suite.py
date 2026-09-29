"""性能基准套件。

目标：在算法等价前提下，量化 MaxTiC-Next 在**内部节点数 n** 与**约束数 |E|** 上的可扩展性。

计时与计数口径：

* **CPU 时间**用 ``time.process_time()``（进程消耗的 CPU 时间），**墙钟时间**用
  ``time.perf_counter()`` 单独报告 —— 墙钟与 CPU 时间必须名实相符；
* **只计时排序内核**：输出文件写入不在计时窗口内（用同一份内容**另行**写一次来测 I/O
  开销，得到 ``io_time``，并把 ``kernel_cpu_time = cpu_time - io_time`` 一并报告）；
  ``tracemalloc`` 也不在计时窗口内——峰值内存由**独立的内存测量趟**取得，避免把
  追踪开销算进耗时；
* **预热 + 重复**：默认 1 次预热（不计入）+ ``n_repeats``（默认 3）次计时重复，
  报告**中位数与四分位距 (IQR)**，而不是单次测量；
* **诚实计数**：横轴的 ``|E|`` 报告"请求数 / 实际解析出的约束行数 / 聚合后的有效边键数 /
  被树种系哨兵占用的边数"，不再把注释行、FRQ 头行算作约束，也不再"有放回抽样"后
  把聚合掉的重复对当成新约束；
* **种子有效**：合成树的形状与标签确实由 ``seed`` 决定；
* **正确性断言**：除计时外还独立重算交付序的 ``best_value``（与 ``values["best"]``、
  与冲突约束权重和、与"不劣于任何启发式"三点对比），记为 ``correctness_ok``。

数据集为合成的**二叉树**（n 个内部节点、n+1 个叶子）与随机约束（|E| 条），不依赖真实
ALE 或任何外部数据，保证可复现（受 ``seed`` 驱动）。

此外，本模块还提供**真实数据集基准**：

* :func:`run_real_benchmark`：对单组真实 tree + constraints 运行完整排序；
* :func:`run_real_suite`：扫描目录下所有 tree+constraints 对，批量运行。

用法：
    python -m maxtic_next.benchmarks.suite            # 运行默认合成扫描
    python -m maxtic_next.benchmarks.suite --max-n 200 --max-e 2000
    python -m maxtic_next.benchmarks.suite --repeats 5 --warmup 2
    python -m maxtic_next.benchmarks.suite --real                     # 用 examples/ 真实数据
    python -m maxtic_next.benchmarks.suite --real --data-dir /path/to/datasets
"""

import argparse
import math
import os
import sys
import time
import tracemalloc
from typing import Any, Dict, List, Optional, Sequence, Tuple

from maxtic_next.config import OUTPUT_SUFFIXES
from maxtic_next.constraints.constraint import Constraint, ConstraintSet
from maxtic_next.random_ import RandomWrapper
from maxtic_next.ranking.ranker import Ranker
from maxtic_next.tree.tree import Tree


def _median_iqr(samples: Sequence[float]) -> Tuple[Optional[float], Optional[float], Optional[Any]]:
    """返回 ``(中位数, IQR, (Q1, Q3))``（纯标准库，样本为空时全 ``None``）。"""
    xs = sorted(float(x) for x in samples)
    if not xs:
        return None, None, None

    def _quantile(p: float) -> float:
        if len(xs) == 1:
            return xs[0]
        pos = p * (len(xs) - 1)
        lo = int(math.floor(pos))
        hi = min(lo + 1, len(xs) - 1)
        frac = pos - lo
        return xs[lo] * (1.0 - frac) + xs[hi] * frac

    q1, med, q3 = _quantile(0.25), _quantile(0.5), _quantile(0.75)
    return med, q3 - q1, (q1, q3)


def _stats(samples: Sequence[float], prefix: str) -> Dict[str, Any]:
    """把一组重复测量整理成 ``{prefix}_median`` / ``_iqr`` / ``_min`` / ``_max`` / 原样本。"""
    med, iqr, _ = _median_iqr(samples)
    return {
        f"{prefix}_median": med,
        f"{prefix}_iqr": iqr,
        f"{prefix}_min": (min(samples) if samples else None),
        f"{prefix}_max": (max(samples) if samples else None),
        f"{prefix}_samples": [float(s) for s in samples],
    }


def make_synthetic_tree(n_internal: int, seed: int = 0) -> Tuple[str, List[str]]:
    """生成一棵含 ``n_internal`` 个内部节点的二叉树 Newick 与内部节点标签列表。

    ：``seed`` **确实生效**——它同时决定 (a) 每个内部节点左右子树的规模分配
    （在平衡点附近抖动，抖动幅度受控以保证树不至于退化成链条）与 (b) 内部节点标签的
    编号顺序。不同的 ``seed`` 给出不同拓扑，相同 ``seed`` 逐字符可复现。

    Args:
        n_internal: 内部节点数量（平衡二叉树，叶子数 = n_internal + 1）。
        seed: 随机种子（驱动树形与标签命名，保证可复现）。

    Returns:
        (newick, internal_labels)，其中 internal_labels 为该树全部内部节点 bootstrap 标签。
    """
    rng = RandomWrapper(seed * 7 + 1)
    leaves_pool = [f"L{k}" for k in range(2 * (n_internal + 1) + 1)]
    leaf_idx = [0]

    def _newick_subtree(nodes_left: int) -> str:
        # 生成一棵恰好含 nodes_left 个内部节点的子树
        if nodes_left <= 0:
            name = leaves_pool[leaf_idx[0]]
            leaf_idx[0] += 1
            return name
        # 把剩余内部节点分给左右子树：以平衡点为中心、受 seed 驱动的抖动
        half = (nodes_left - 1) // 2
        jitter = min(half, max(1, nodes_left // 8))
        left = half + rng.randint(2 * jitter + 1) - jitter
        left = max(0, min(nodes_left - 1, left))
        right = nodes_left - 1 - left
        return f"({_newick_subtree(left)},{_newick_subtree(right)})"

    # 先用占位构建结构，再回填内部标签
    structure = _newick_subtree(n_internal)
    # seed 驱动的标签置换：同一 seed 逐字符可复现，不同 seed 改变"哪个号在哪个位置"
    label_pool = [str(k) for k in range(1, n_internal + 1)]
    rng.shuffle(label_pool)
    counter = [0]
    internal_labels: List[str] = []

    def _relabel(s: str) -> str:
        # 把结构中的空内部节点（"()" 形式）替换为带标签的 "(...)<label>"
        # 采用递归解析：遇到 '(' 进入子树，遇到 ')' 之后回填标签
        i = 0
        n = len(s)

        def _parse() -> str:
            # 解析一个子树，返回带标签的 Newick 片段
            nonlocal i
            i += 1  # 跳过 '('
            children: List[str] = []
            while i < n and s[i] != ")":
                if s[i] == "(":
                    children.append(_parse())
                elif s[i] == ",":
                    i += 1
                else:
                    # 叶子名
                    start = i
                    while i < n and s[i] not in ",)":
                        i += 1
                    children.append(s[start:i])
            i += 1  # 跳过 ')'
            k = counter[0]
            counter[0] += 1
            label = label_pool[k] if k < len(label_pool) else str(k + 1)
            internal_labels.append(label)
            return "(" + ",".join(children) + ")" + label

        return _parse()

    newick = _relabel(structure) + ";"
    return newick, internal_labels


def comparable_pairs(labels: Sequence[str], tree: Optional[Tree] = None) -> List[Tuple[str, str]]:
    """列出**不受树种系拓扑强制**的有序节点对（可用作信息性约束的方向对）。

    ``tree is None`` 时返回全部有序对；给定树时剔除祖先-后代对（这些对的方向在任何
    合法排序中都相同，写进约束也会被置为 ``MAX_NUMBER`` 哨兵，属于"假 |E|"）。
    """
    nodes = list(labels)
    anc: Optional[set] = None
    if tree is not None:
        label_to_id = {}
        for nid in tree.get_nodes():
            if not tree.is_leaf(nid):
                label_to_id[tree.get_bootstrap(nid)] = nid
        anc = set()
        for a in nodes:
            for b in nodes:
                if (
                    a != b
                    and a in label_to_id
                    and b in label_to_id
                    and tree.is_ancestor(label_to_id[a], label_to_id[b])
                ):
                    anc.add((a, b))
    pairs: List[Tuple[str, str]] = []
    for a in nodes:
        for b in nodes:
            if a == b:
                continue
            if anc is not None and ((a, b) in anc or (b, a) in anc):
                continue
            pairs.append((a, b))
    return pairs


def make_synthetic_constraints(
    internal_labels: List[str],
    n_edges: int,
    seed: int = 0,
    tree: Optional[Tree] = None,
    stats_out: Optional[Dict] = None,
) -> ConstraintSet:
    """生成至多 ``n_edges`` 条随机约束（donor/receptor 取自内部节点标签，权重 1.0）。

    ：改为**无放回**抽样（不再"有放回取对"后靠聚合去重，导致 |E|=1000 实际
    只有 501 条边键），并可在给定 ``tree`` 时只取**不被拓扑强制**的方向对（否则这些
    约束会被置成 ``MAX_NUMBER`` 树边，是"假规模"）。可用对不足时按实际数量生成，
    真实口径写入 ``stats_out``。

    Args:
        internal_labels: 内部节点标签列表。
        n_edges: 期望的约束条数。
        seed: 随机种子。
        tree: 可选物种树；提供时剔除祖先-后代平凡对。
        stats_out: 可选字典，回填 ``n_edges_requested`` / ``n_edges_generated`` /
            ``n_pairs_available`` / ``n_edge_keys_effective``。
    """
    rng = RandomWrapper(seed * 13 + 5)
    cset = ConstraintSet()
    pairs = comparable_pairs(internal_labels, tree)
    # 打乱后再顺序取，等价"无放回均匀抽样"且可复现
    order = list(range(len(pairs)))
    rng.shuffle(order)
    picked = order[: min(n_edges, len(pairs))]
    for idx in picked:
        a, b = pairs[idx]
        cset.add(
            Constraint(
                donor=a,
                receptor=b,
                weight=1.0,
                metadata={"family": "synthetic", "support": 1.0, "distance": None},
            )
        )
    if stats_out is not None:
        stats_out["n_pairs_available"] = len(pairs)
        stats_out["n_edges_requested"] = n_edges
        stats_out["n_edges_generated"] = len(cset.constraints)
        stats_out["n_edge_keys_effective"] = len({(c.donor, c.receptor) for c in cset.constraints})
    return cset


def verify_delivered_value(result) -> Dict:
    """**独立重算**交付序的目标值并做三点一致性断言（不能只断言计时）。

    MaxTiC 的目标函数是"该序被违反的约束权重和"。本函数不信任 ``Result.values``，
    而是从 ``informative_lines`` **与** ``conflicting_lines``（同为
    ``"donor,receptor weight"`` 行）重建边集，用 ``ranking.value.value`` 重算
    ``best_order`` 的值，再检查三件事：

    1. 与 ``values["best"]``（若上层已给出）一致；
    2. 与 ``conflicting_lines`` 的权重和一致（被违反集合就是冲突文件的内容）；
    3. **不劣于**任何启发式给出的候选值（``input`` / ``greedy`` / ``mixing`` /
       ``local_search``）——择优与局部搜索没有把结果改坏。

    为什么也要读 ``conflicting_lines``：``informative_lines`` 只含"信息性"约束，
    而 ``value()`` / 冲突文件还涵盖**平凡冲突**（donor 是受体后代的"反向"约束，任何
    合法序都会违反它）。只用 informative 重建边集会在含此类约束的真实数据上得到
    "重算值 < 冲突权重和"，把**正确**的结果误判为 BAD。两边取并集仍是"从交付序独立
    重算"，判定力不降（伪造的冲突行与序不符时依然对不上）。

    Returns:
        dict，含各检查项与总判定 ``ok``；``ok is False`` 表示"计时看起来正常、
        但算出来的东西不对"，基准必须让它失败。
    """
    from maxtic_next.ranking.value import value as _value

    edge: Dict[str, float] = {}
    keys: List[str] = []
    n_informative = 0
    for source_attr, is_informative in (("informative_lines", True), ("conflicting_lines", False)):
        for line in getattr(result, source_attr, []) or []:
            parts = str(line).split()
            if len(parts) < 2:
                continue
            key = parts[0]
            try:
                weight = float(parts[1])
            except ValueError:
                continue
            if key not in edge:
                edge[key] = weight
                keys.append(key)
                if is_informative:
                    n_informative += 1
            elif is_informative:
                n_informative += 1
    delivered = _value(list(result.best_order), edge, keys)
    conflict_sum = 0.0
    for line in getattr(result, "conflicting_lines", []) or []:
        parts = str(line).split()
        if len(parts) >= 2:
            try:
                conflict_sum += float(parts[1])
            except ValueError:
                pass
    values = getattr(result, "values", {}) or {}
    candidates = [
        v
        for k, v in values.items()
        if k in ("input", "greedy", "mixing", "local_search")
        and isinstance(v, (int, float))
        and math.isfinite(v)
    ]
    reported = values.get("best")
    tol = 1e-6 * max(1.0, abs(delivered))
    checks = {
        "recomputed_value": delivered,
        "reported_best_value": reported,
        "conflicting_weight_sum": conflict_sum,
        "informative_count": n_informative,
        "edges_considered": len(keys),
        "finite_value": bool(math.isfinite(delivered)),
        "matches_reported_best": (
            not isinstance(reported, (int, float)) or abs(delivered - reported) <= tol
        ),
        "matches_conflicting_sum": abs(delivered - conflict_sum) <= tol,
        "not_worse_than_heuristics": (not candidates or delivered <= min(candidates) + tol),
    }
    checks["ok"] = bool(
        checks["finite_value"]
        and checks["matches_reported_best"]
        and checks["matches_conflicting_sum"]
        and checks["not_worse_than_heuristics"]
    )
    return checks


def _make_ranker(
    tree: Tree, cset: ConstraintSet, seed: int, output_prefix: str, **overrides
) -> Ranker:
    """基准用的固定配置 ``Ranker``（关闭 LS / 随机树 / MCMC，只测排序内核）。"""
    kwargs: Dict[str, Any] = dict(
        tree=tree,
        cset=cset,
        rng=RandomWrapper(seed),
        min_transfer_distance=0,
        threshold_constraints=0.0,
        temperature=0.001,
        random_type=0,
        time_for_search=0.0,
        random_trees=0,
        constraint_file="benchmark",
        output_prefix=output_prefix,
    )
    kwargs.update(overrides)
    return Ranker(**kwargs)


def _measure_io_seconds(result, dest_prefix: str) -> Tuple[float, float]:
    """把**同一份内容**再写一次，测出输出 I/O 的 ``(墙钟, CPU)`` 秒数。

    基准关心的是排序内核的可扩展性，但 ``Ranker.run()`` 内联写三文件；此处在计时窗口
    **之外**重放同样的写入，得到可从总耗时中扣除的 I/O 成本估计。
    """
    from maxtic_next.io.output import write_three_files

    os.makedirs(os.path.dirname(dest_prefix) or ".", exist_ok=True)
    cpu0, wall0 = time.process_time(), time.perf_counter()
    write_three_files(
        dest_prefix,
        list(result.informative_lines),
        list(result.conflicting_lines),
        list(result.partial_lines),
        output_style="short",
        force=True,
    )
    return (time.perf_counter() - wall0, time.process_time() - cpu0)


def _measure_peak_mb(make_inputs, prefix_root: str, seed: int) -> Optional[float]:
    """在**独立的一趟**里用 ``tracemalloc`` 测峰值内存（不污染计时窗口）。"""
    tree, cset, size_stats = make_inputs()
    sub = os.path.join(prefix_root, "mem")
    os.makedirs(sub, exist_ok=True)
    tracemalloc.start()
    try:
        ranker = _make_ranker(tree, cset, seed, os.path.join(sub, "benchmark"))
        ranker.run(print_summary=False, html_report=False)
        _cur, peak = tracemalloc.get_traced_memory()
    finally:
        tracemalloc.stop()
    return peak / (1024 * 1024)


def run_benchmark(
    n_internal: int,
    n_edges: int,
    seed: int = 0,
    time_budget: float = 30.0,
    output_prefix: Optional[str] = None,
    n_repeats: int = 3,
    warmup: int = 1,
    **ranker_overrides,
) -> Dict:
    """在给定 (n, |E|) 上运行基准，返回 CPU / 墙钟（中位数 + IQR）、内存与正确性判定。

    Args:
        n_internal: 内部节点数。
        n_edges: **请求**的约束数 |E|；实际生成/有效边键数见返回的
            ``n_edges_generated`` / ``n_edge_keys_effective``。
        seed: 随机种子（同时驱动合成树形状与约束抽样）。
        time_budget: 单次墙钟预算（秒）；中位数超过则标记 ``exceeded``。
        output_prefix: 输出目录前缀；``None`` 时用系统临时目录并在结束后清理。
        n_repeats: 计时重复次数（>= 1），报告**中位数与 IQR**。
        warmup: 预热次数（其结果**不**进入统计，仅消除首次导入/缓存效应）。
        **ranker_overrides: 透传给 ``Ranker`` 的额外参数（如 ``time_for_search``）。

    Returns:
        含以下键的 dict（保留旧键名以兼容既有调用方/测试）：
        ``n_internal`` / ``n_edges``（请求数）/ ``cpu_time``（= CPU 中位数，
        ``time.process_time``）/ ``wall_time`` / ``cpu_time_iqr`` / ``wall_time_iqr`` /
        ``io_time`` / ``kernel_cpu_time`` / ``n_repeats`` / ``warmup`` /
        ``peak_memory_mb`` / ``best_value`` / ``delivered_value`` /
        ``correctness_ok`` / ``correctness`` / ``exceeded`` / ``error`` /
        ``n_edges_generated`` / ``n_edge_keys_effective`` / ``n_pairs_available``。
    """
    import shutil
    import tempfile

    result: Dict = {
        "n_internal": n_internal,
        "n_edges": n_edges,
        "cpu_time": None,
        "wall_time": None,
        "cpu_time_median": None,
        "cpu_time_iqr": None,
        "wall_time_median": None,
        "wall_time_iqr": None,
        "io_time_median": None,
        "kernel_cpu_time": None,
        "n_repeats": max(1, int(n_repeats)),
        "warmup": max(0, int(warmup)),
        "peak_memory_mb": None,
        "best_value": None,
        "delivered_value": None,
        "correctness_ok": None,
        "correctness": None,
        "exceeded": False,
        "error": None,
    }

    _tmpdir_to_clean = None
    if output_prefix is None:
        _tmpdir_to_clean = tempfile.mkdtemp(prefix="maxtic_bench_")
        prefix_root = os.path.join(_tmpdir_to_clean, "out")
    else:
        prefix_root = os.path.join(output_prefix + "_bench_repeats")
    os.makedirs(prefix_root, exist_ok=True)

    def make_inputs():
        newick, labels = make_synthetic_tree(n_internal, seed=seed)
        tree = Tree()
        tree.read_newick(newick)
        stats: Dict = {}
        cset = make_synthetic_constraints(labels, n_edges, seed=seed, tree=tree, stats_out=stats)
        return tree, cset, stats

    cpu_samples: List[float] = []
    wall_samples: List[float] = []
    io_wall_samples: List[float] = []
    io_cpu_samples: List[float] = []
    last_res = None
    try:
        tree, cset, size_stats = make_inputs()
        result["n_edges_requested"] = n_edges
        result["n_edge_keys_effective"] = size_stats.get("n_edge_keys_effective")
        result["n_edges_generated"] = size_stats.get("n_edges_generated")
        result["n_pairs_available"] = size_stats.get("n_pairs_available")

        # ---- 预热（不计入统计）----
        for w in range(max(0, int(warmup))):
            sub = os.path.join(prefix_root, f"warmup{w}")
            os.makedirs(sub, exist_ok=True)
            ranker = _make_ranker(
                tree, cset, seed, os.path.join(sub, "benchmark"), **ranker_overrides
            )
            ranker.run(print_summary=False, html_report=False)

        # ---- 计时重复 ----
        for r in range(max(1, int(n_repeats))):
            sub = os.path.join(prefix_root, f"rep{r}")
            os.makedirs(sub, exist_ok=True)
            ranker = _make_ranker(
                tree, cset, seed, os.path.join(sub, "benchmark"), **ranker_overrides
            )
            cpu0 = time.process_time()
            wall0 = time.perf_counter()
            res = ranker.run(print_summary=False, html_report=False)
            wall = time.perf_counter() - wall0
            cpu = time.process_time() - cpu0
            cpu_samples.append(cpu)
            wall_samples.append(wall)
            last_res = res
            io_wall, io_cpu = _measure_io_seconds(
                res, os.path.join(prefix_root, f"io_rep{r}", "benchmark")
            )
            io_wall_samples.append(io_wall)
            io_cpu_samples.append(io_cpu)

        result.update(_stats(cpu_samples, "cpu_time"))
        result.update(_stats(wall_samples, "wall_time"))
        result.update(_stats(io_wall_samples, "io_time"))
        med_cpu, med_iqr, _ = _median_iqr(cpu_samples)
        med_wall, wall_iqr, _ = _median_iqr(wall_samples)
        med_io = _median_iqr(io_cpu_samples)[0]
        result["cpu_time"] = med_cpu  # 真正的 CPU 时间（process_time 中位数）
        result["wall_time"] = med_wall
        result["cpu_time_iqr"] = med_iqr
        result["wall_time_iqr"] = wall_iqr
        result["io_time"] = med_io
        result["kernel_cpu_time"] = (
            None if med_cpu is None or med_io is None else max(0.0, med_cpu - med_io)
        )
        result["exceeded"] = bool(med_wall is not None and med_wall > time_budget)

        # ---- 峰值内存：独立一趟，tracemalloc 不进计时窗口 ----
        result["peak_memory_mb"] = _measure_peak_mb(make_inputs, prefix_root, seed)

        # ---- 正确性断言（不只是"跑完了没报错"）----
        if last_res is not None:
            checks = verify_delivered_value(last_res)
            result["correctness"] = checks
            result["correctness_ok"] = checks["ok"]
            result["delivered_value"] = checks["recomputed_value"]
            result["best_value"] = min(
                last_res.values.get("greedy", float("inf")),
                last_res.values.get("mixing", float("inf")),
            )
    except Exception as exc:  # noqa: BLE001 - 基准需捕获一切异常
        result["error"] = f"{type(exc).__name__}: {exc}"
    finally:
        if _tmpdir_to_clean is not None:
            shutil.rmtree(_tmpdir_to_clean, ignore_errors=True)
    return result


def run_suite(
    n_list: List[int],
    e_list: List[int],
    seed: int = 0,
    time_budget: float = 30.0,
    n_repeats: int = 3,
    warmup: int = 1,
) -> List[Dict]:
    """在 (n, |E|) 网格上运行基准，返回结果列表。"""
    rows: List[Dict] = []
    for n in n_list:
        for e in e_list:
            rows.append(
                run_benchmark(
                    n, e, seed=seed, time_budget=time_budget, n_repeats=n_repeats, warmup=warmup
                )
            )
    return rows


# ----------------------------------------------------------------------
# 真实数据集基准
# ----------------------------------------------------------------------
# 真实数据集默认目录：项目根目录下的 examples/
_PROJECT_ROOT = os.path.dirname(
    os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
)
_DEFAULT_REAL_DATA_DIR = os.path.join(_PROJECT_ROOT, "examples")

# 树文件扩展名与约束文件扩展名
_TREE_EXTS = (".tree", ".nwk", ".newick", ".tre")
_CONST_EXTS = (".tsv", ".txt", ".constraints", ".csv")

# MaxTiC 自身产出的输出文件后缀（legacy + short 两种命名风格，均从 config 派生），
# 数据集扫描须排除之，避免把输出当作输入约束重新摄入（尤其 short 风格以 .tsv 结尾）。
_OUTPUT_SUFFIXES_ALL = tuple(s.lower() for group in OUTPUT_SUFFIXES.values() for s in group) + (
    ".html",
)


def _is_maxtic_output(name: str) -> bool:
    """判断文件名是否为 MaxTiC 输出（三文件 / distribution / HTML 报告）。"""
    low = name.lower()
    return any(low.endswith(suf) for suf in _OUTPUT_SUFFIXES_ALL)


def _constraint_stats(path: str) -> Dict[str, int]:
    """约束文件的**诚实规模**。

    旧实现数"非空行数"，于是 ``#`` 注释行与 ALE 的 ``FRQ`` 头行也被当成约束（实测
    报 4 条、真实 2 条），"规模–耗时"曲线的横轴因此失真。现用真正的解析器计数，
    并同时报告去重后的不同 ``(donor, receptor)`` 对数（聚合进 ``edge`` 的实际边键数）。
    """
    out = {
        "n_lines_nonempty": 0,
        "n_constraints_parsed": 0,
        "n_distinct_pairs": 0,
        "n_comment_or_header_lines": 0,
    }
    try:
        with open(path, "r", encoding="utf-8") as fh:
            lines = fh.readlines()
    except (OSError, UnicodeDecodeError):
        return out
    nonempty = [ln for ln in lines if ln.strip()]
    out["n_lines_nonempty"] = len(nonempty)
    out["n_comment_or_header_lines"] = len(nonempty) - sum(
        1 for ln in nonempty if not ln.lstrip().startswith("#") and "FRQ" not in ln
    )
    try:
        from maxtic_next.constraints.parsers import parse_constraints_file

        cset = parse_constraints_file(path)
        out["n_constraints_parsed"] = len(cset.constraints)
        out["n_distinct_pairs"] = len({(c.donor, c.receptor) for c in cset.constraints})
    except Exception:  # noqa: BLE001 - 解析失败退回保守的行计数口径
        out["n_constraints_parsed"] = sum(
            1 for ln in nonempty if not ln.lstrip().startswith("#") and "FRQ" not in ln
        )
        out["n_distinct_pairs"] = out["n_constraints_parsed"]
    return out


def _count_constraints(path: str) -> int:
    """约束文件中的**有效约束条数**（解析后，不含 ``#`` 注释与 FRQ 头行）。"""
    return _constraint_stats(path)["n_constraints_parsed"]


def _count_internal_nodes(tree_path: str) -> int:
    """统计 Newick 树文件中的内部节点数。"""
    try:
        from maxtic_next.tree.tree import Tree

        tree = Tree()
        with open(tree_path, "r", encoding="utf-8") as fh:
            tree.read_newick(fh.readline())
        return len(tree.internal_node_labels())
    except Exception:  # noqa: BLE001
        return 0


def find_dataset_pairs(data_dir: str, quiet: bool = False) -> List[Tuple[str, str]]:
    """扫描目录，找到所有 (tree_file, constraints_file) 对。

    配对策略（按优先级）：
    1. 名称匹配（去掉扩展名后相同）：``foo.tree`` + ``foo.tsv``
    2. 剩余文件若树数 == 约束数 -> 按名序一对一配对（策略 1 没吃掉的）
    3. 剩余树与剩余约束数量不等 -> 取笛卡尔积
    4. **只剩树**（全部约束都已被同名配对吃掉）-> 每棵剩余树与**全部**约束文件配对
       （名称不匹配的兜底），而不是把树丢掉

    ：旧实现在"已经出现任一同名配对"后就**静默丢弃**其余未匹配文件
    （``elif not pairs`` 分支），导致批量基准的数据组数被悄悄少算。现在
    **任何**树都会进入配对，且每一次"非同名兜底"都在 stderr 明说配对关系是
    怎么来的（含文件名），批量基准的行数因此可核对；确实想只跑同名配对时，
    请把目录整理成一对一命名或用 ``--data-dir`` 指向子目录。

    Args:
        data_dir: 数据集目录。
        quiet: 不打 stderr 的配对说明（诊断信息仍照常返回）。

    Returns:
        (tree_path, constraints_path) 元组列表（按路径排序，确定性）。
    """
    if not os.path.isdir(data_dir):
        return []

    tree_files = []
    cons_files = []
    for name in sorted(os.listdir(data_dir)):
        full = os.path.join(data_dir, name)
        if not os.path.isfile(full):
            continue
        low = name.lower()
        if _is_maxtic_output(low):
            continue  # 跳过 MaxTiC 自身输出，避免被当作输入约束
        if low.endswith(_TREE_EXTS):
            tree_files.append(full)
        elif low.endswith(_CONST_EXTS):
            cons_files.append(full)

    if not tree_files or not cons_files:
        return []

    # 策略 1：按去扩展名后的名称匹配
    def _stem(path):
        return os.path.splitext(os.path.basename(path))[0]

    cons_by_stem: Dict[str, List[str]] = {}
    for cf in cons_files:
        cons_by_stem.setdefault(_stem(cf), []).append(cf)

    pairs: List[Tuple[str, str]] = []
    used_trees = set()
    used_cons = set()
    for tf in tree_files:
        for cf in cons_by_stem.get(_stem(tf), []):
            pairs.append((tf, cf))
            used_trees.add(tf)
            used_cons.add(cf)

    # 策略 2 / 3 / 4：剩余文件一律参与配对，绝不静默丢弃
    remaining_trees = [t for t in tree_files if t not in used_trees]
    remaining_cons = [c for c in cons_files if c not in used_cons]
    fallback: List[Tuple[str, str]] = []
    mode = ""
    if remaining_trees and remaining_cons:
        if len(remaining_trees) == len(remaining_cons):
            mode = "按名序一对一"
            fallback.extend(zip(sorted(remaining_trees), sorted(remaining_cons)))
        else:
            mode = "笛卡尔积（剩余树 × 剩余约束）"
            for tf in remaining_trees:
                for cf in remaining_cons:
                    fallback.append((tf, cf))
    elif remaining_trees:
        # 只剩树（所有约束都已被同名配对吃掉）：与**全部**约束文件配对，
        # 否则这些树会被静默丢弃 —— 正是
        mode = "兜底（与目录内全部约束文件配对）"
        for tf in remaining_trees:
            for cf in cons_files:
                fallback.append((tf, cf))
    if fallback and not quiet:
        print(
            f"[benchmark] {data_dir}：{len(remaining_trees)} 棵树"
            f"（{', '.join(os.path.basename(t) for t in remaining_trees)}）"
            f"没有同名约束文件，已按{mode}新增 {len(fallback)} 组配对。"
            "这种配对**不保证**树与约束来自同一数据集，批量基准的结论请先核对"
            "配对关系（或把目录整理成同名一对一）。",
            file=sys.stderr,
            flush=True,
        )
    pairs.extend(fallback)
    # 只剩约束（无树可配）时同样要说出来，而不是静默少算数据组
    if remaining_cons and not remaining_trees and not quiet:
        print(
            f"[benchmark] {data_dir}：{len(remaining_cons)} 个约束文件没有任何"
            f"树可配（已跳过）："
            f"{', '.join(os.path.basename(c) for c in remaining_cons)}",
            file=sys.stderr,
            flush=True,
        )
    return sorted(set(pairs))


def run_real_benchmark(
    tree_path: str,
    constraints_path: str,
    seed: int = 42,
    time_budget: float = 120.0,
    output_prefix: Optional[str] = None,
    n_repeats: int = 3,
    warmup: int = 1,
) -> Dict:
    """使用真实数据运行完整排序，测量 CPU / 墙钟时间（中位数 + IQR）与峰值内存。

    与 :func:`run_benchmark`（合成数据）不同，本函数接受真实的物种树
    (Newick) 与约束文件路径，通过 ``rank`` API 运行完整排序流程。计时口径与合成
    基准一致：CPU 用 ``time.process_time``、墙钟单独报告、输出 I/O 在
    计时窗口外另行测量并扣除、``tracemalloc`` 只用于独立的一趟内存测量，并对交付序
    做独立重算的 ``best_value`` 正确性断言。

    Args:
        tree_path: 物种树文件路径（Newick）。
        constraints_path: 约束文件路径（TSV / 文本双格式）。
        seed: 随机种子。
        time_budget: 单次墙钟预算（秒）；中位数超过则标记 ``exceeded``。
        output_prefix: 输出文件前缀；``None`` 时使用临时目录。
        n_repeats: 计时重复次数（>= 1）。
        warmup: 预热次数（不计入统计）。

    Returns:
        含 ``tree_path`` / ``constraints_path`` / ``n_internal`` /
        ``n_constraints`` / ``n_distinct_pairs`` / ``cpu_time`` / ``wall_time`` /
        ``cpu_time_iqr`` / ``io_time`` / ``kernel_cpu_time`` / ``peak_memory_mb`` /
        ``best_value`` / ``delivered_value`` / ``correctness_ok`` / ``best_source`` /
        ``similarity_to_input`` / ``exceeded`` / ``error`` 的 dict。
    """
    from maxtic_next.api import rank

    result: Dict = {
        "tree_path": tree_path,
        "constraints_path": constraints_path,
        "n_internal": _count_internal_nodes(tree_path),
        "n_constraints": _count_constraints(constraints_path),
        "cpu_time": None,
        "wall_time": None,
        "cpu_time_median": None,
        "cpu_time_iqr": None,
        "wall_time_median": None,
        "wall_time_iqr": None,
        "io_time": None,
        "io_time_median": None,
        "kernel_cpu_time": None,
        "n_repeats": max(1, int(n_repeats)),
        "warmup": max(0, int(warmup)),
        "peak_memory_mb": None,
        "best_value": None,
        "delivered_value": None,
        "correctness_ok": None,
        "correctness": None,
        "best_source": None,
        "similarity_to_input": None,
        "exceeded": False,
        "error": None,
    }
    # 诚实规模口径：注释/FRQ 头行不算约束，同时给出去重后的边键数
    stats = _constraint_stats(constraints_path)
    result.update({f"cons_{k}": v for k, v in stats.items()})
    result["n_distinct_pairs"] = stats["n_distinct_pairs"]

    import tempfile
    import shutil

    _tmpdir_to_clean = None
    if output_prefix is None:
        _tmpdir_to_clean = tempfile.mkdtemp(prefix="maxtic_real_bench_")
        prefix_root = os.path.join(_tmpdir_to_clean, "real")
    else:
        prefix_root = output_prefix + "_bench_repeats"
    os.makedirs(prefix_root, exist_ok=True)

    cpu_samples: List[float] = []
    wall_samples: List[float] = []
    io_wall: List[float] = []
    io_cpu: List[float] = []
    last_res = None
    try:

        def _one(prefix: str):
            return rank(
                tree_path,
                constraints_path,
                seed=seed,
                output_prefix=prefix,
                print_summary=False,
                html_report=False,
            )

        for w in range(max(0, int(warmup))):
            sub = os.path.join(prefix_root, f"warmup{w}")
            os.makedirs(sub, exist_ok=True)
            _one(sub)

        for r in range(max(1, int(n_repeats))):
            sub = os.path.join(prefix_root, f"rep{r}")
            os.makedirs(sub, exist_ok=True)
            cpu0, wall0 = time.process_time(), time.perf_counter()
            res = _one(sub)
            wall = time.perf_counter() - wall0
            cpu = time.process_time() - cpu0
            cpu_samples.append(cpu)
            wall_samples.append(wall)
            last_res = res
            iw, ic = _measure_io_seconds(res, os.path.join(prefix_root, f"io_rep{r}", "real"))
            io_wall.append(iw)
            io_cpu.append(ic)

        result.update(_stats(cpu_samples, "cpu_time"))
        result.update(_stats(wall_samples, "wall_time"))
        result.update(_stats(io_wall, "io_time"))
        med_cpu, med_iqr, _ = _median_iqr(cpu_samples)
        med_wall, wall_iqr, _ = _median_iqr(wall_samples)
        med_io = _median_iqr(io_cpu)[0]
        result["cpu_time"] = med_cpu
        result["wall_time"] = med_wall
        result["cpu_time_iqr"] = med_iqr
        result["wall_time_iqr"] = wall_iqr
        result["io_time"] = med_io
        result["kernel_cpu_time"] = (
            None if med_cpu is None or med_io is None else max(0.0, med_cpu - med_io)
        )
        result["exceeded"] = bool(med_wall is not None and med_wall > time_budget)

        # 峰值内存：独立一趟，tracemalloc 不进计时窗口
        mem_sub = os.path.join(prefix_root, "mem")
        os.makedirs(mem_sub, exist_ok=True)
        tracemalloc.start()
        try:
            _one(mem_sub)
            _cur, peak = tracemalloc.get_traced_memory()
        finally:
            tracemalloc.stop()
        result["peak_memory_mb"] = peak / (1024 * 1024)

        if last_res is not None:
            checks = verify_delivered_value(last_res)
            result["correctness"] = checks
            result["correctness_ok"] = checks["ok"]
            result["delivered_value"] = checks["recomputed_value"]
            greedy = last_res.values.get("greedy", float("inf"))
            mixing = last_res.values.get("mixing", float("inf"))
            result["best_value"] = min(greedy, mixing)
            result["best_source"] = last_res.best_source
            result["similarity_to_input"] = last_res.similarity_to_input
    except Exception as exc:  # noqa: BLE001 - 基准需捕获一切异常
        result["error"] = f"{type(exc).__name__}: {exc}"
    finally:
        if _tmpdir_to_clean is not None:
            shutil.rmtree(_tmpdir_to_clean, ignore_errors=True)
    return result


def run_real_suite(
    data_dir: str, seed: int = 42, time_budget: float = 120.0, n_repeats: int = 3, warmup: int = 1
) -> List[Dict]:
    """扫描目录下所有 tree+constraints 对，批量运行真实数据基准。

    Args:
        data_dir: 数据集目录路径。目录内应包含 Newick 树文件
            (``*.tree`` / ``*.nwk`` / ``*.newick``) 与约束文件
            (``*.tsv`` / ``*.txt``)。
        seed: 随机种子。
        time_budget: 单次墙钟预算（秒）。
        n_repeats: 每组数据的计时重复次数（报告中位数 + IQR）。
        warmup: 每组数据的预热次数（不计入统计）。

    Returns:
        每组数据的基准结果 dict 列表（同 :func:`run_real_benchmark` 返回格式）。
    """
    pairs = find_dataset_pairs(data_dir)
    if not pairs:
        print(f"警告：目录 {data_dir} 中未找到 tree+constraints 数据对。")
        print(f"  树文件扩展名: {_TREE_EXTS}")
        print(f"  约束文件扩展名: {_CONST_EXTS}")
        return []

    results: List[Dict] = []
    for tree_path, cons_path in pairs:
        name = os.path.basename(tree_path)
        print(f"  运行: {name} + {os.path.basename(cons_path)} ...", end=" ", flush=True)
        r = run_real_benchmark(
            tree_path,
            cons_path,
            seed=seed,
            time_budget=time_budget,
            n_repeats=n_repeats,
            warmup=warmup,
        )
        if r["error"]:
            print(f"ERROR: {r['error']}")
        elif r["exceeded"]:
            print(f"EXCEEDED (wall {r['wall_time']:.1f}s)")
        else:
            print(
                f"ok (cpu {r['cpu_time']:.3f}s ±IQR {r['cpu_time_iqr']:.3f}, "
                f"wall {r['wall_time']:.3f}s, "
                f"{r['peak_memory_mb']:.2f}MB, "
                f"value={r['delivered_value']}, "
                f"correct={'yes' if r['correctness_ok'] else 'NO'})"
            )
        results.append(r)
    return results


def _fmt(x: Any, spec: str = ".3f") -> str:
    """按格式打印数值，``None`` / 异常值统一显示为 ``-``。"""
    if x is None:
        return "-"
    try:
        return format(x, spec)
    except (TypeError, ValueError):
        return str(x)


def _correctness_cell(r: Dict) -> str:
    """正确性列：``ok`` / ``BAD``（基准不能只断言计时）。"""
    if r.get("error"):
        return "-"
    if r.get("correctness_ok") is None:
        return "-"
    return "ok" if r["correctness_ok"] else "BAD"


def _format_real_table(rows: List[Dict]) -> str:
    """格式化真实数据基准结果表（CPU = ``process_time`` 中位数；墙钟与 IQR 分列）。"""
    header = (
        f"{'dataset':>30} {'n_int':>6} {'|E|':>8} {'|E|distinct':>12} "
        f"{'CPU(s)':>10} {'IQR':>9} {'wall(s)':>9} {'kernel(s)':>10} "
        f"{'PeakMB':>10} {'best_val':>12} {'correct':>8} {'sim':>8} "
        f"{'status':>12}"
    )
    lines = [header, "-" * len(header)]
    for r in rows:
        name = os.path.basename(r.get("tree_path", "?"))[:30]
        if r["error"]:
            status = "ERROR"
            cpu = iqr = wall = kernel = mem = bval = sim = "-"
        else:
            status = "EXCEEDED" if r["exceeded"] else "ok"
            cpu = _fmt(r.get("cpu_time"))
            iqr = _fmt(r.get("cpu_time_iqr"))
            wall = _fmt(r.get("wall_time"))
            kernel = _fmt(r.get("kernel_cpu_time"))
            mem = _fmt(r.get("peak_memory_mb"), ".2f")
            bval = _fmt(r.get("delivered_value", r.get("best_value")), ".2f")
            sim = _fmt(r.get("similarity_to_input"), ".4f")
        lines.append(
            f"{name:>30} {r['n_internal']:>6} {r['n_constraints']:>8} "
            f"{r.get('n_distinct_pairs', r['n_constraints']):>12} "
            f"{cpu:>10} {iqr:>9} {wall:>9} {kernel:>10} {mem:>10} {bval:>12} "
            f"{_correctness_cell(r):>8} {sim:>8} {status:>12}"
        )
    return "\n".join(lines)


def _format_table(rows: List[Dict]) -> str:
    """格式化合成扫描结果表：n / 请求|E| / 实际|E| / CPU 中位数 / IQR / 墙钟 / 内核 / 内存。"""
    header = (
        f"{'n':>6} {'|E|req':>8} {'|E|real':>8} {'CPU(s)':>10} {'IQR':>9} "
        f"{'wall(s)':>9} {'kernel(s)':>10} {'io(s)':>8} {'PeakMB':>10} "
        f"{'best_val':>10} {'correct':>8} {'reps':>5} {'status':>12}"
    )
    lines = [header, "-" * len(header)]
    for r in rows:
        if r["error"]:
            status = "ERROR"
            cpu = iqr = wall = kernel = io_t = mem = bval = "-"
        else:
            status = "EXCEEDED" if r["exceeded"] else "ok"
            cpu = _fmt(r.get("cpu_time"))
            iqr = _fmt(r.get("cpu_time_iqr"))
            wall = _fmt(r.get("wall_time"))
            kernel = _fmt(r.get("kernel_cpu_time"))
            io_t = _fmt(r.get("io_time"), ".4f")
            mem = _fmt(r.get("peak_memory_mb"), ".2f")
            bval = _fmt(r.get("delivered_value", r.get("best_value")))
        n_real = r.get("n_edge_keys_effective", r.get("n_edges_generated", r["n_edges"]))
        lines.append(
            f"{r['n_internal']:>6} {r['n_edges']:>8} {n_real:>8} {cpu:>10} "
            f"{iqr:>9} {wall:>9} {kernel:>10} {io_t:>8} {mem:>10} {bval:>10} "
            f"{_correctness_cell(r):>8} {r.get('n_repeats', 1):>5} {status:>12}"
        )
    return "\n".join(lines)


def main(argv: Optional[List[str]] = None) -> int:
    """基准 CLI 入口。"""
    parser = argparse.ArgumentParser(description="MaxTiC-Next 性能基准（纯标准库）")
    parser.add_argument("--max-n", type=int, default=200, help="最大内部节点数")
    parser.add_argument("--max-e", type=int, default=2000, help="最大约束数")
    parser.add_argument("--seed", type=int, default=0, help="随机种子（驱动合成树形状与约束抽样）")
    parser.add_argument("--time-budget", type=float, default=30.0, help="单次墙钟时间预算（秒）")
    parser.add_argument("--steps", type=int, default=4, help="n 与 |E| 的扫描步数")
    parser.add_argument(
        "--repeats", type=int, default=3, help="计时重复次数（>= 1，报告中位数与 IQR）"
    )
    parser.add_argument(
        "--warmup", type=int, default=1, help="预热次数（不计入统计，消除首次导入/缓存效应）"
    )
    parser.add_argument(
        "--real", action="store_true", help="使用真实数据集进行基准测试（而非合成数据）"
    )
    parser.add_argument(
        "--data-dir",
        type=str,
        default=None,
        help="真实数据集目录路径（配合 --real 使用）；默认使用项目根目录下的 examples/",
    )
    args = parser.parse_args(argv)
    if args.repeats < 1:
        parser.error(f"--repeats 必须 >= 1，收到 {args.repeats}")
    if args.warmup < 0:
        parser.error(f"--warmup 必须 >= 0，收到 {args.warmup}")

    # ---- 真实数据集基准模式 ----
    if args.real:
        data_dir = args.data_dir or _DEFAULT_REAL_DATA_DIR
        print(
            f"真实数据集基准：data_dir={data_dir}, seed={args.seed}, "
            f"墙钟预算={args.time_budget}s, 重复={args.repeats}, "
            f"预热={args.warmup}"
        )
        rows = run_real_suite(
            data_dir,
            seed=args.seed,
            time_budget=args.time_budget,
            n_repeats=args.repeats,
            warmup=args.warmup,
        )
        if rows:
            print()
            print(_format_real_table(rows))
            feasible = [r for r in rows if not r["exceeded"] and not r["error"]]
            bad = [r for r in rows if r["correctness_ok"] is False]
            if bad:
                print(
                    f"\n错误：{len(bad)} 组数据的交付值独立重算不一致（correct=BAD），基准不可信。"
                )
            if feasible:
                print(f"\n真实数据基准完成：{len(feasible)}/{len(rows)} 组数据在预算内完成。")
            else:
                print("\n警告：所有真实数据组均超出时间预算或出错。")
            return 1 if bad else 0
        return 0

    # ---- 合成数据基准模式（默认） ----
    n_list = [max(1, args.max_n // (2**k)) for k in range(args.steps)][::-1]
    n_list = sorted(set(n_list))
    e_list = [max(1, args.max_e // (2**k)) for k in range(args.steps)][::-1]
    e_list = sorted(set(e_list))

    print(
        f"基准扫描：n in {n_list}, |E| in {e_list}, 墙钟预算 {args.time_budget}s, "
        f"重复 {args.repeats} 次（预热 {args.warmup} 次，取中位数 ± IQR）"
    )
    rows = run_suite(
        n_list,
        e_list,
        seed=args.seed,
        time_budget=args.time_budget,
        n_repeats=args.repeats,
        warmup=args.warmup,
    )
    print(_format_table(rows))

    # 标出可处理上限：在预算内完成、且交付值通过正确性断言的最大 (n, |E|) 组合
    feasible = [
        r
        for r in rows
        if not r["exceeded"] and not r["error"] and r.get("correctness_ok") is not False
    ]
    if feasible:
        best = max(feasible, key=lambda r: (r["n_internal"], r["n_edges"]))
        print(
            f"可处理上限（预算内且交付值正确）：n={best['n_internal']}, "
            f"|E|={best['n_edges']}（有效边键 {best.get('n_edge_keys_effective')}）"
        )
    else:
        print("警告：所有组合均超出时间预算、出错或未通过交付值正确性断言。")
    bad = [r for r in rows if r["correctness_ok"] is False]
    if bad:
        print(f"错误：{len(bad)} 个组合的交付值独立重算不一致（correct=BAD）。")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
