"""输出层：三文件写入与 stdout 摘要。

三个输出文件名按 ``output_style`` 选择（默认 ``short``）：

* ``short``（默认，简洁）：``<prefix>.mt.informative.tsv`` /
  ``<prefix>.mt.conflicts.tsv`` / ``<prefix>.mt.partial_order.tsv``；
* ``legacy``（与原版逐字节等价的长名）：
  ``<prefix>_MT_output_filtered_list_of_weighted_informative_constraints`` 等。

两种风格下三文件的**内容逐字节一致**，仅文件名不同；stdout 摘要不含文件名，
逐行复刻原版 ``MaxTiC.py`` 的 print 顺序与措辞（不受命名风格影响）。

写入语义：

* **不静默覆盖**：目标文件已存在且不是本进程本次会话写出的路径时，抛
  ``FileExistsError``；调用方须显式传 ``force=True``（CLI 的 ``--force``）；
* **原子写入**：先写同目录 ``*.tmp.<rand>`` 临时文件，``flush + os.fsync`` 后
  ``os.replace`` 到目标路径，中途崩溃不会留下"半截"产物；
* **缺失父目录自动创建**；创建失败给出带路径的清晰错误。
"""

import os
import tempfile
from typing import Iterable, List, Tuple

from maxtic_next.config import DEFAULT_OUTPUT_STYLE, OUTPUT_SUFFIXES

# 本进程已写出（因而可覆写）的产物路径集合。
# 语义：同一次运行内重跑（例如对同一 Ranker 连续调用 ``run()``）允许覆写自己
# 刚刚产出的文件；跨进程覆盖用户已有产物必须显式 ``force=True``。
_WRITTEN_BY_THIS_PROCESS: set = set()


def reset_output_ownership_for_testing() -> None:
    """清空"本进程已写出"登记（仅供测试使用）。"""
    _WRITTEN_BY_THIS_PROCESS.clear()


def _ensure_parent_dir(path: str) -> None:
    """确保 ``path`` 的父目录存在；失败时给出带路径的清晰错误。"""
    parent = os.path.dirname(os.path.abspath(path))
    if not parent or os.path.isdir(parent):
        return
    try:
        os.makedirs(parent, exist_ok=True)
    except OSError as exc:  # noqa: BLE001 - 转为面向用户的清晰信息
        raise OSError(
            f"无法创建输出目录 {parent!r}（输出前缀 {path!r}）：{exc}。请检查路径拼写与写权限。"
        ) from exc


def _ensure_not_overwritten(path: str, force: bool) -> None:
    """覆盖守卫：目标已存在、且不是本进程写出的，必须显式 ``--force`` 才允许覆盖。

    ``_WRITTEN_BY_THIS_PROCESS`` 存绝对路径，因此查表两侧都必须归一化为绝对路径——
    否则相对前缀在同一进程内重跑时会被误判为"别人的既有结果"。
    """
    key = os.path.abspath(path)
    if os.path.exists(path) and key not in _WRITTEN_BY_THIS_PROCESS and not force:
        raise FileExistsError(
            f"输出文件已存在：{path}\n"
            "为避免静默覆盖既有结果，已中止。请改用 -p/--output-prefix 指定其它前缀，"
            "或加 --force 显式覆盖。"
        )


def _write_atomic(path: str, body: str) -> None:
    """把 ``body`` 原子写入 ``path``（临时文件 + fsync + ``os.replace``）。"""
    _ensure_parent_dir(path)
    parent = os.path.dirname(os.path.abspath(path))
    fd, tmp_path = tempfile.mkstemp(prefix=".maxtic-tmp-", dir=parent)
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as fh:
            fh.write(body)
            fh.flush()
            os.fsync(fh.fileno())
        os.replace(tmp_path, path)
    except BaseException:
        # 失败时清掉临时文件，绝不留半截产物
        try:
            os.unlink(tmp_path)
        except OSError:
            pass
        raise
    _WRITTEN_BY_THIS_PROCESS.add(os.path.abspath(path))


def _atomic_dump(path: str, lines: Iterable[str], force: bool = False) -> None:
    """把 ``lines`` 原子写入 ``path``（临时文件 + ``os.replace``）。

    Args:
        path: 目标路径。
        lines: 逐行内容（不含行尾换行）。
        force: 目标已存在且非本进程写出时是否覆盖。

    Raises:
        FileExistsError: 目标已存在、非本进程写出、且 ``force`` 为假。
        OSError: 目录创建或写入失败。
    """
    _ensure_not_overwritten(path, force)
    _write_atomic(path, "".join(line + "\n" for line in lines))


def write_text_file_atomic(path: str, text: str, force: bool = False) -> str:
    """原子写入一个**整块文本**产物（如自包含 HTML 报告），并复用同一套覆盖守卫。

    HTML 报告与三件套共用同一套 ``--force`` 语义，避免同一条命令下
    "三件套拒绝覆盖、报告却直接覆盖"的两套口径：既有文件需要 ``--force``，
    写入过程原子（不留半截 HTML）。
    """
    _ensure_not_overwritten(path, force)
    _write_atomic(path, text)
    return path


def write_three_files(
    prefix: str,
    informative_lines: List[str],
    conflicting_lines: List[str],
    partial_lines: List[str],
    output_style: str = DEFAULT_OUTPUT_STYLE,
    force: bool = False,
) -> Tuple[str, str, str]:
    """写入三个输出文件，返回其路径三元组。

    Args:
        prefix: 输出文件前缀（通常为约束文件路径）。
        informative_lines: 信息性约束行（``"key weight"``）。
        conflicting_lines: 与最优序冲突的约束行（``"key weight"``）。
        partial_lines: 偏序行（``"a b weight black/green"``）。
        output_style: 命名风格 ``"short"``（默认）或 ``"legacy"``。
        force: 允许覆盖已存在的产物文件（对应 CLI ``--force``）。

    Returns:
        (informative_file, conflicting_file, partial_order_file) 路径三元组。
    """
    suf_filtered, suf_conflicting, suf_partial, _ = OUTPUT_SUFFIXES[output_style]
    informative_file = prefix + suf_filtered
    conflicting_file = prefix + suf_conflicting
    partial_file = prefix + suf_partial

    _atomic_dump(informative_file, informative_lines, force=force)
    _atomic_dump(conflicting_file, conflicting_lines, force=force)
    _atomic_dump(partial_file, partial_lines, force=force)
    return informative_file, conflicting_file, partial_file


def write_aux_file(path: str, lines: Iterable[str], force: bool = False) -> str:
    """原子写入一个辅助产物（如 ``--random-trees`` 的能量/相似度分布）。"""
    _atomic_dump(path, lines, force=force)
    return path


__all__ = ["write_three_files", "write_aux_file", "write_text_file_atomic", "format_summary"]


def _g(x: float) -> str:
    """按原版 ``str(float)`` 的方式格式化浮点数（Python 2/3 均为最短往返表示）。

    原版 ``MaxTiC.py``（Python 2）的 ``print`` 对浮点使用 ``str()``，
    Python 2 与 Python 3 的 ``str(float)`` 行为一致：都采用**最短往返（shortest
    round-trip）**表示。例如 ``str(50.0)`` 在 Py2/Py3 都得到 ``'50.0'``，
    与原版 stdout **逐字节一致**。

    注意 ``"%.12g" % x`` 会偏离原版：``"%.12g" % 50.0`` 得到 ``'50'``、
    大数可能被科学计数法表示，都会破坏与原版的字节级等价。因此这里严格使用
    ``str(x)`` 复刻原版 ``print`` 的浮点显示。
    """
    return str(x)


def format_summary(result) -> str:
    """生成与原版 print 顺序/措辞一致的 stdout 摘要字符串。

    ``result`` 为 ``maxtic_next.ranking.ranker.Result`` 实例（惰性导入以避免循环依赖）。

    与原版 ``MaxTiC.py`` 的**已记录偏离**（，非回归）：

    1. ``uninformative`` 百分比按 ``uninformative / total_weight`` 计（
       原版分母 ``total_transfers + uninformative`` 把无信息权重量双重计入）；
    2. 开启局部搜索时补打原版的 ``after local search ... rejected`` /
       ``best found solution ... (%)`` 两行，使 stdout 不与交付树矛盾；
    3. 空白：原版 ``print a,"uninformative (",b`` 产生**两个**空格，本版与原版
       对齐；
    4. 谱系硬约束 NOTE 的措辞固定为"…并已计入**上方** uninformative 统计（在
       total_weight 分母中只计一次，不重复计权）"。该 NOTE 由 ``Ranker`` 生成、
       经本函数 appended 到摘要**末尾**，故统计行在其**上方**，措辞用"计入上方"；
       "已计入 / 只计一次"是事实口径，
       与 ``tests/baselines/stdout.txt`` 一致。
    """
    from maxtic_next.ranking.ranker import Result  # noqa: F401 (type check only)

    r: Result = result

    total = r.total_weight
    lines: List[str] = []
    lines.append(r.constraint_file)
    # 原版：print "tree with ", N, "internal nodes" -> "tree with  N internal nodes"
    lines.append(f"tree with  {r.internal_node_count} internal nodes")

    lines.append(f"{_g(total)} total weight of constraints from transfers")

    # 分母为 total_weight（其本身已包含 uninformative），不双重计权
    uninf = r.uninformative["total"]
    pct_uninf = (uninf * 100 / total) if total != 0 else 0.0
    lines.append(
        f"{_g(uninf)} uninformative (  {_g(r.uninformative['to_desc'])} to a descendant,"
        f" {_g(r.uninformative['to_leaf'])} to a leaf {_g(r.uninformative['to_anc'])} to an ancestor"
        f" {_g(r.uninformative['to_itself'])} to itself {pct_uninf:.1f}%)"
    )
    if r.from_leaf > 0:
        lines.append(
            "WARNING: there are constraints from leaves, which cannot be met since this"
            " program only orders internal nodes, they are ignored in the scores"
        )
    pct_triv = int(r.trivial_conflict * 100 / total) if total != 0 else 0
    lines.append(
        f"{_g(r.trivial_conflict)} trivially conflicting constraints ({pct_triv}%)"
        " (descendant to ancestor or trivial cycle)"
    )

    def _value_line(label: str, val: float) -> str:
        if total != 0:
            return f"{label} {_g(val)} ({_g(val * 100 / total)}%)"
        return f"{label} {_g(val)}"

    lines.append(_value_line("value of the order given by the input tree", r.values["input"]))
    lines.append(_value_line("value of the greedy heuristic", r.values["greedy"]))
    lines.append(_value_line("value of the mixing heuristic:", r.values["mixing"]))
    lines.append(f"best order is the {r.best_source}")

    # 局部搜索后的真实值与交付序一致地打印（原版 MaxTiC.py:668-669）
    if "local_search" in r.values:
        best_val = r.values.get("best", r.values["local_search"])
        lines.append(f"after local search {_g(r.values['local_search'])} rejected")
        lines.append(_value_line("best found solution", best_val))

    lines.append(r.ranked_newick)
    lines.append(
        f"Similarity of the best order compared with the input order {_g(r.similarity_to_input)}"
    )

    ratio = (r.conflict_with_input / r.partial_total) if r.partial_total != 0 else 0.0
    lines.append(
        f"{_g(ratio)} constraints in agreement with the best tree in conflict with input tree"
    )
    lines.append(
        f"Similarity of the order compared with the input order {_g(r.similarity_to_input)}"
    )

    # --random-trees 统计摘要（等价原版 RANDOM 分支的 4 行打印）
    if r.random_stats is not None:
        rs = r.random_stats
        lines.append(
            f"values from  {rs['n']}  random orders {_g(rs['value_min'])} {_g(rs['value_max'])}"
        )
        # (k+1)/(n+1) 校正，且检验对象为最终交付序（含局部搜索）
        lines.append(
            f"pvalue of the found order: {_g(rs['value_pvalue'])}"
            " [(k+1)/(n+1) corrected; tested order = delivered best order]"
        )
        lines.append(
            f"similarity values from  {rs['n']} random orders"
            f" {_g(rs['sim_min'])} {_g(rs['sim_max'])}"
        )
        lines.append(
            f"pvalue of the similarity with the input order: {_g(rs['sim_pvalue'])}"
            " [(k+1)/(n+1) corrected]"
        )

    # 运行期告警（口径不闭合、选项被忽略、不可复现等）
    for w in getattr(r, "warnings", []) or []:
        lines.append(w)

    return "\n".join(lines)
