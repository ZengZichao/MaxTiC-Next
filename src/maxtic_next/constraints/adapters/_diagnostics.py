"""上游适配器共享的**诊断结构**与**诊断异常**。

为什么需要本模块：五个上游格式适配器很容易把"静默丢弃"当作常规控制流——
标签没命中、支持度不够、家族规模探测失败、XML 解析崩掉，全都表现为"该工具没有
检出转移"，用户无从分辨。本模块提供三样东西：

1. ``new_diagnostics`` / ``merge_diagnostics`` / ``summary_line``：一个**机器可读**的
   计数字典，随 ``ConstraintSet.diagnostics`` 返回，并在 ``convert()`` 结束时打印
   一行 stderr 摘要；
2. ``donor_endpoint_convention`` 与 ``weight_semantics```两个口径标签：
   前者记录"供体端点取哪一层"（ALE = ``parent_of_donor``，其余四工具 =
   ``donor_itself``），后者记录权重到底是 [0,1] 支持度、整数计数还是后验频率；
3. ``AdapterDiagnosticError`` 及其子类（均继承 ``ValueError``，故调用方原有的
   ``except ValueError`` 语义不变）：把"必须让用户知道"的情形从静默返回空集合
   改为抛出带上下文的异常。

设计约束：纯标准库、无第三方依赖、所有结构均可 pickle（``ProcessPoolExecutor``
子进程会返回它）。
"""

import os
import sys
from typing import Any, Dict, Iterable, List, Optional, Set

from maxtic_next.config import DEFAULT_MIN_ENDPOINT_HIT_RATE  # noqa: F401 (再导出)

# ----------------------------------------------------------------------
# 供体端点层级约定（不同工具相差一个物种层级，不可直接混用）
# ----------------------------------------------------------------------
#: ALE：约束 donor 取"供体的父"（``parent(donnor)`` / ``donnor_search``）
CONVENTION_PARENT_OF_DONOR = "parent_of_donor"
#: RANGER-DTLx / ARTra / ecceTERA / AleRax：donor 取"供体本身"
#: （``Mapping -->`` / ``branchingOut speciesLocation`` / 表文件 donor 列）
CONVENTION_DONOR_ITSELF = "donor_itself"

# ----------------------------------------------------------------------
# 权重口径（块数不是采样分母）
# ----------------------------------------------------------------------
#: ``count / 工具自己声明的样本数``（ALE 的 ``N reconciled``、
#: RANGER-DTLx / ARTra 的 ``Total number of optimal solutions: N``）
WEIGHT_DECLARED_SUPPORT = "support_over_declared_samples"
#: ``count / 调和块数``：仅当各块叶子集合一致（确为同一家族的重复调和）时才成立
WEIGHT_BLOCK_SUPPORT = "support_over_replicate_blocks"
#: 分母不可确立：权重为**整数计数**，不是 [0,1] 支持度
WEIGHT_INTEGER_COUNT = "integer_count"
#: AleRax：上游已按后验样本预聚合的频率/期望计数，本适配器不再归一化
WEIGHT_POSTERIOR_FREQUENCY = "posterior_frequency"

#: 适配器文件级缓存键的版本号（键必须含物种树**拓扑**指纹）。
#: 任何使缓存载荷语义改变的改动都要递增此值，使旧条目永不被读到。
CACHE_KEY_VERSION = 2

#: 端点命中率低于该比例即告警；命中率为 0 时直接报错。
#: 数值本身定义在 ``config.DEFAULT_MIN_ENDPOINT_HIT_RATE``（CLI / API / 适配器共用
#: 同一默认值，：默认值不得在多处各写一遍），此处仅再导出。

_INT_KEYS = (
    "files_seen",
    "blocks_seen",
    "trees_seen",
    "transfers_seen",
    "transfers_resolved",
    "transfers_unpaired",
    "numeric_id_resolutions",
    "constraints_kept",
    "dropped_by_label_miss",
    "dropped_by_support",
    "dropped_by_family_size",
    "dropped_by_family_size_probe_failed",
    "self_loop_transfers",
    "malformed_lines",
    "families_skipped",
)
_LIST_KEYS = (
    "declared_sample_counts",
    "unresolved_endpoint_samples",
    "species_label_samples",
    "warnings",
    "errors",
)


def new_diagnostics(
    tool: str, convention: str, weight_semantics: str = WEIGHT_INTEGER_COUNT
) -> Dict:
    """返回一个零值诊断字典（键集合固定，便于机器读取与合并）。"""
    diag: Dict = {
        "tool": tool,
        "donor_endpoint_convention": convention,
        "weight_semantics": weight_semantics,
        "sample_denominator": None,
        "family_size_probe_ok": None,  # None=未知, True/False=探测结论
        "endpoint_hit_rate": None,
        "min_support_is_fractional": weight_semantics
        in (WEIGHT_DECLARED_SUPPORT, WEIGHT_BLOCK_SUPPORT),
    }
    for k in _INT_KEYS:
        diag[k] = 0
    for k in _LIST_KEYS:
        diag[k] = []
    return diag


def merge_diagnostics(dst: Dict, src: Optional[Dict]) -> Dict:
    """把 ``src`` 并入 ``dst``（整数求和、列表拼接去重、标量取后者）。

    ``sample_denominator`` 特殊：跨文件求和（同一批输入里每个家族各自带来一段
    采样），其余标量（口径标签 / 命中率等）以最后一次观测为准。
    """
    if not src:
        return dst
    for k, v in src.items():
        if k in _INT_KEYS:
            dst[k] = int(dst.get(k, 0) or 0) + int(v or 0)
        elif k in _LIST_KEYS:
            cur = dst.setdefault(k, [])
            for item in v or []:
                if item not in cur:
                    cur.append(item)
        elif k == "sample_denominator":
            continue
        elif v is None:
            continue
        else:
            dst[k] = v
    denom = src.get("sample_denominator")
    if denom:
        dst["sample_denominator"] = int(dst.get("sample_denominator") or 0) + int(denom)
    return dst


def reciprocal_subtract(counts: Dict) -> None:
    """互反减法（各适配器共用，原地修改）：若互反对权重更小，从较大者减去较小者。

    等价原版 ``constraints_from_reconciliations.py`` 的减法。**例外**：
    ``"X,X"`` 这种"自互反"键（供体 = 受体的自环事件）不参与减法——
    原版会把 ``counts["X,X"] -= counts["X,X"]`` 归零，使"供体=受体谱系"这一类
    事件从统计里静默消失；此处保留它，并由各适配器的 ``self_loop_transfers``
    计数上报（下游是否按"to itself"处理与原版一致）。
    """
    for k in list(counts.keys()):
        if k not in counts:
            continue
        first, _, second = k.partition(",")
        if first == second:
            continue  # 自环不参与互反抵消

        opp = second + "," + first
        if opp in counts and counts[opp] <= counts[k]:
            counts[k] = counts[k] - counts[opp]
            counts[opp] = counts[opp] * 0


def endpoint_hit_rate(diag: Dict) -> Optional[float]:
    """端点命中率 = 解析成功的转移配对 / 看到的转移配对总数。"""
    seen = int(diag.get("transfers_seen", 0))
    if seen <= 0:
        return None
    return float(diag.get("transfers_resolved", 0)) / float(seen)


def weight_and_support(diag: Dict, count: float) -> tuple:
    """按本批次权重口径，把一个事件计数换算成 ``(weight, support)``。

    * 口径为支持度且分母可确立：``weight = support = count / denominator``（[0,1]）；
    * 口径为整数计数 / 后验频率：``weight = support = count``
      （分母不可确立时**不假装**自己是 [0,1] 支持度）。
    """
    denom = diag.get("sample_denominator")
    if diag.get("weight_semantics") in (WEIGHT_DECLARED_SUPPORT, WEIGHT_BLOCK_SUPPORT) and denom:
        value = float(count) / float(denom)
        return value, value
    return float(count), float(count)


def passes_support(diag: Dict, weight: float, count: float, min_support: float) -> bool:
    """阈值判定：一律按 ``weight > min_support`` 比较，但**口径随权重语义而变**。

    * 支持度口径（``support_over_declared_samples`` /
      ``support_over_replicate_blocks``）：``weight`` 是 [0,1] 支持度，阈值有意义；
    * AleRax 的 ``posterior_frequency``：``weight`` 是后验频率，阈值在同一量纲上比较；
    * ``integer_count``（分母不可确立）：``weight`` 是**整数计数**，此时阈值不再是
      [0,1] 支持度（默认 0.05 于是退化为"至少出现一次"）——因此这里把
      ``min_support_is_fractional`` 置为 ``False``，并由 ``finalize`` 打印警告，
      使"阈值成为空操作"这件事**说出来**而不是静默。
    """
    if diag.get("weight_semantics") == WEIGHT_INTEGER_COUNT:
        diag["min_support_is_fractional"] = False
    return weight > min_support


def warn(message: str, tool: str = "adapter") -> None:
    """向 stderr 产出一行警告（不抛异常、不影响 stdout 逐字节口径）。"""
    print(f"[{tool}] 警告：{message}", file=sys.stderr, flush=True)


def error_notice(message: str, tool: str = "adapter") -> None:
    """向 stderr 产出一行**错误级**提示（已抛异常时无需调用）。"""
    print(f"[{tool}] 错误：{message}", file=sys.stderr, flush=True)


def summary_line(diag: Dict) -> str:
    """把诊断字典压成**一行**可读摘要（stderr 用）。"""
    denom = diag.get("sample_denominator")
    hit = diag.get("endpoint_hit_rate")
    hit_s = "" if hit is None else f" 端点命中率={hit * 100.0:.1f}%"
    return (
        f"[{diag.get('tool', 'adapter')}] donor 端点约定={diag.get('donor_endpoint_convention')}"
        f" 权重口径={diag.get('weight_semantics')}"
        f" 分母={'不可确立' if denom in (None, 0) else denom}"
        f"{hit_s}"
        f" | 文件={diag.get('files_seen', 0)} 块={diag.get('blocks_seen', 0)}"
        f" 转移={diag.get('transfers_seen', 0)}"
        f" 保留={diag.get('constraints_kept', 0)}"
        f" 丢弃[标签不匹配={diag.get('dropped_by_label_miss', 0)},"
        f" 支持度={diag.get('dropped_by_support', 0)},"
        f" 家族规模={diag.get('dropped_by_family_size', 0)},"
        f" 家族规模探测失败={diag.get('dropped_by_family_size_probe_failed', 0)},"
        f" 自环={diag.get('self_loop_transfers', 0)}]"
    )


def emit_summary(diag: Dict, tool: str = "adapter") -> str:
    """产出诊断摘要行（stderr），并返回该行（便于测试断言）。"""
    line = summary_line(diag)
    print(line, file=sys.stderr, flush=True)
    return line


def sample_of(items: Iterable, n: int = 8) -> List[str]:
    """取前 ``n`` 个样本标签（排序后），用于错误信息里"两边命名各长什么样"。"""
    out = sorted({str(x) for x in items})
    return out[:n]


# ----------------------------------------------------------------------
# 诊断异常（全部继承 ``ValueError``：调用方既有 ``except ValueError`` 仍成立）
# ----------------------------------------------------------------------


class AdapterDiagnosticError(ValueError):
    """适配器携带诊断信息的错误基类。"""

    def __init__(self, message: str, diagnostics: Optional[Dict] = None) -> None:
        super().__init__(message)
        self.diagnostics: Dict = diagnostics or {}


class LabelMismatchError(AdapterDiagnosticError):
    """上游端点标签与物种树命名不兼容。

    典型场景：ecceTERA 在 recPhyloXML 里用**数字节点 ID** 标注内部物种节点，
    而用户物种树的内部标签是别的命名（或反之）。
    """


class UpstreamParseError(AdapterDiagnosticError):
    """上游文件解析失败（XML 异常不再被吞成"无约束"）。"""

    def __init__(
        self,
        message: str,
        source: str = "",
        line: int = 0,
        column: int = 0,
        cause: Optional[BaseException] = None,
        diagnostics: Optional[Dict] = None,
    ) -> None:
        loc = f"（文件 {source}"
        if line:
            loc += f"，第 {line} 行"
            if column:
                loc += f"，第 {column} 列"
        loc += "）" if source else "）"
        full = message + (loc if (source or line) else "")
        if cause is not None:
            full += f"；原始异常：{type(cause).__name__}: {cause}"
        super().__init__(full, diagnostics)
        self.source = source
        self.line = line
        self.column = column
        self.cause = cause


def map_payloads(
    worker, args_list: List, max_workers: Optional[int] = None, parallel: str = "process"
) -> List:
    """在各适配器共用的 ``ProcessPoolExecutor`` / ``ThreadPoolExecutor`` 上跑
    ``worker(arg)``，失败自动回退顺序执行；结果**严格按输入顺序**返回。

    并行只影响速度不影响结果（各文件解析互独立，边键聚合与顺序无关）。
    worker 必须是模块级可 pickle 函数，参数只含简单类型。
    """
    n = len(args_list)
    payloads: List = [None] * n
    if n == 0:
        return payloads

    executor_cls: Optional[type] = None
    mp_context: Optional[Any] = None
    if parallel == "process":
        try:
            from concurrent.futures import ProcessPoolExecutor
            import multiprocessing

            executor_cls = ProcessPoolExecutor
            # macOS 默认 spawn 会重导入 __main__（pytest 下会死锁），
            # 故优先 fork；Windows 无 fork 时回退默认上下文。
            try:
                mp_context = multiprocessing.get_context("fork")
            except ValueError:
                mp_context = None
        except Exception:
            executor_cls = None
    elif parallel == "thread":
        try:
            from concurrent.futures import ThreadPoolExecutor

            executor_cls = ThreadPoolExecutor
        except Exception:
            executor_cls = None

    if executor_cls is None:
        return [worker(arg) for arg in args_list]

    try:
        if max_workers is None:
            max_workers = max(1, min(os.cpu_count() or 1, 8, n))
        kwargs: Dict = {"max_workers": max_workers}
        if mp_context is not None:
            kwargs["mp_context"] = mp_context
        with executor_cls(**kwargs) as executor:
            future_to_idx = {executor.submit(worker, arg): i for i, arg in enumerate(args_list)}
            for future in future_to_idx:
                idx = future_to_idx[future]
                try:
                    payloads[idx] = future.result()
                except Exception:
                    payloads[idx] = worker(args_list[idx])
    except Exception:
        payloads = [worker(arg) for arg in args_list]
    return payloads


def check_endpoint_hit_rate(
    diag: Dict, *, species_labels: Set[str], tool: str, threshold: float, raise_on_zero: bool = True
) -> Optional[float]:
    """命中率守卫：0 命中报错、低于阈值告警。

    Args:
        diag: 诊断字典（需含 transfers_seen / transfers_resolved /
            unresolved_endpoint_samples）。
        species_labels: 物种树合法标签样本（写进错误信息里做对照）。
        tool: 工具名（stderr 前缀）。
        threshold: 低于该命中率即告警；``<= 0`` 表示关闭守卫。
        raise_on_zero: 命中率为 0 时是否抛 ``LabelMismatchError``。

    Returns:
        命中率（无转移可判定时为 ``None``）。
    """
    rate = endpoint_hit_rate(diag)
    diag["endpoint_hit_rate"] = rate
    if rate is None:
        return rate
    # ``threshold <= 0`` 只关闭"命中率低于阈值"的**告警**，绝不关闭
    # "端点全部未命中"这一硬错误 —— 后者与阈值无关，且 CLI 帮助承诺的就是
    # "设为 0 关闭告警（0 命中仍报错）"。因此这里不能在 rate <= 0 时提前返回，
    # 否则 ``--min-endpoint-hit-rate 0`` 会把命名不兼容变成静默的 0 约束。
    if rate <= 0.0:
        msg = (
            f"{tool} 输出中的 {diag.get('transfers_seen', 0)} 次转移事件"
            "端点**全部**不在物种树标签集中，"
            "故产出 0 条约束——这是命名不兼容"
            "（如上游用数字节点 ID 标注内部物种节点，而物种树用别的命名），"
            "并不代表“该工具没检出转移”。"
            f" 未命中标本：{diag.get('unresolved_endpoint_samples') or []}；"
            f" 物种树标签样本：{sorted(species_labels)[:8]}"
            f"（内部节点另需与上游 numbering 对齐）。"
        )
        if raise_on_zero:
            raise LabelMismatchError(msg, diag)
        error_notice(msg, tool)
        return rate
    if threshold is None or threshold <= 0:
        # 阈值守卫已按要求关闭；0 命中的硬错误在上方仍会触发。
        return rate
    if rate < threshold:
        msg = (
            f"{tool} 输出的转移事件中仅 {rate * 100.0:.1f}% 端点命中物种树标签"
            f"（{diag.get('transfers_resolved', 0)}/{diag.get('transfers_seen', 0)}，"
            f"低于阈值 {threshold * 100.0:.0f}%）；未命中标本："
            f"{diag.get('unresolved_endpoint_samples') or []}，物种树标签样本："
            f"{sorted(species_labels)[:8]}。约束集可能被静默缩小。"
        )
        diag["warnings"].append(msg)
        warn(msg, tool)
    return rate


def finalize_adapter_diagnostics(
    cset,
    threshold: float,
    species_labels: Set[str],
    tool: str,
    quiet: bool = False,
    files: int = 0,
    notice: Optional[str] = None,
    raise_on_zero_hit: bool = True,
) -> Dict:
    """把逐文件诊断并入 ``cset.diagnostics``、跑守卫、打印 stderr 摘要。

    ：任何"约束变少"都必须可归因（标签不匹配 /
    支持度 / 家族规模 / 自环 / 互反抵消），并且当上游命名与物种树完全
    不兼容时直接报错而不是返回空集合。

    Args:
        cset: 目标 ``ConstraintSet``（鸭子类型：需有 ``diagnostics`` /
            ``constraints`` / ``informative_edges()``）。逐文件诊断应已放在
            ``cset.diagnostics["files"]`` 列表中。
        threshold: 端点命中率告警阈值（``<= 0`` 关闭守卫）。
        species_labels: 物种树标签集合（写进错误/警告信息做对照）。
        tool: 工具名（诊断与 stderr 前缀）。
        quiet: 不打印 stderr（诊断仍随结果返回）。
        files: 本次输入文件数。
        notice: 额外的一行 stderr 提示（如 ALE 的  口径声明）。
        raise_on_zero_hit: 命中率为 0 时是否抛错。

    Returns:
        合并后的诊断字典（即 ``cset.diagnostics``）。
    """
    diag: Dict = cset.diagnostics
    per_file = diag.pop("files", []) or []
    diag["per_file"] = per_file
    diag["tool"] = tool
    for st in per_file:
        merge_diagnostics(diag, st)
    diag["files_seen"] = files or len(per_file)
    denominators = [int(st.get("sample_denominator") or 0) for st in per_file]
    diag["sample_denominator"] = sum(denominators) or None
    semantics = {st.get("weight_semantics") for st in per_file if st.get("weight_semantics")}
    if len(semantics) == 1:
        diag["weight_semantics"] = semantics.pop()
    elif len(semantics) > 1:
        diag["weight_semantics"] = "mixed:" + "+".join(sorted(str(s) for s in semantics))
        diag["warnings"].append(
            f"本批输入含多种权重口径 {sorted(map(str, semantics))}，聚合后的权重不可通约。"
        )
    if per_file:
        diag["family_size_probe_ok"] = all(bool(st.get("family_size_probe_ok")) for st in per_file)
    diag["tool"] = tool
    diag["constraints_kept"] = len(cset.constraints)
    diag["edges_kept"] = len(cset.informative_edges())
    check_endpoint_hit_rate(
        diag,
        species_labels=species_labels,
        tool=tool,
        threshold=threshold or 0.0,
        raise_on_zero=raise_on_zero_hit,
    )
    if files and int(diag.get("transfers_seen", 0)) == 0:
        diag["warnings"].append(
            "上游输出中未看到任何转移事件（0 次）：空约束集是字面结果，"
            "但请确认输入确为该工具的完整输出文件。"
        )
    if str(diag.get("weight_semantics")) == WEIGHT_INTEGER_COUNT and not any(
        "整数计数" in str(w) for w in diag.get("warnings", [])
    ):
        diag["warnings"].append(
            "权重为**整数计数**口径（采样分母不可确立）："
            "--min-support / 支持度阈值本次不具 [0,1] 语义。"
        )
    # 阈值语义以**最终**口径为准（逐文件诊断在合并过程中可能被 integer_count 改写）
    diag["min_support_is_fractional"] = diag.get("weight_semantics") in (
        WEIGHT_DECLARED_SUPPORT,
        WEIGHT_BLOCK_SUPPORT,
    ) and bool(diag.get("sample_denominator"))
    for msg in list(diag.get("warnings", [])):
        if not quiet:
            warn(msg, tool)
    if notice and not quiet:
        print(notice, file=sys.stderr, flush=True)
    if not quiet:
        emit_summary(diag, tool)
    return diag


__all__ = [
    "CONVENTION_PARENT_OF_DONOR",
    "CONVENTION_DONOR_ITSELF",
    "WEIGHT_DECLARED_SUPPORT",
    "WEIGHT_BLOCK_SUPPORT",
    "WEIGHT_INTEGER_COUNT",
    "WEIGHT_POSTERIOR_FREQUENCY",
    "CACHE_KEY_VERSION",
    "DEFAULT_MIN_ENDPOINT_HIT_RATE",
    "new_diagnostics",
    "merge_diagnostics",
    "endpoint_hit_rate",
    "weight_and_support",
    "passes_support",
    "reciprocal_subtract",
    "warn",
    "error_notice",
    "summary_line",
    "emit_summary",
    "sample_of",
    "AdapterDiagnosticError",
    "LabelMismatchError",
    "UpstreamParseError",
    "check_endpoint_hit_rate",
    "finalize_adapter_diagnostics",
    "map_payloads",
]
