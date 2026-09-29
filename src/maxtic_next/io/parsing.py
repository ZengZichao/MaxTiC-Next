"""Newick 与约束文件的读取分发。

本层承担"文件 -> 对象"的健壮性职责：

* 以 ``utf-8-sig`` 读取，带 BOM 的文件不再污染首列标签；
* Newick 允许**折行（多行）**书写：解析前把整份文件的行拼接为一条序列；
* 空文件 / 只有空白行的文件给出清晰错误，而非裸 ``IndexError``；
* 约束文件支持**行内 ``#`` 注释**（``61 67 15.48 # 说明``）；
* 约束权重/距离列在此处即校验：**必须有限（拒绝 ``nan`` / ``inf``）且权重 >= 0**，
  报错信息含**文件名与行号**，非法权重不会一路带进目标函数；
* 约束文件与 Newick 物种树都**透明支持 gzip 输入**（``.gz``、按 ``1f 8b``
  魔数识别的 gzip 流、以及恰好一个成员的 ``.tar.gz`` / ``.tgz`` 归档），多成员归档
  给出含 ``tar -xzOf`` 命令的可执行提示（见 :mod:`maxtic_next.io.compression`）。
"""

import math
from typing import List

from maxtic_next.constraints.constraint import ConstraintSet
from maxtic_next.constraints.parsers import parse_constraints
from maxtic_next.io import compression
from maxtic_next.tree.tree import Tree


def read_text_lines(path: str) -> List[str]:
    """以 ``utf-8-sig`` 读取文本文件为行列表（BOM 自动剥离）。

    ：底层走 :mod:`maxtic_next.io.compression`，因此 ``.gz`` / gzip 流与
    **单成员** ``.tar.gz`` 归档无需手工解压即可读取；多成员归档抛
    :class:`~maxtic_next.io.compression.CompressedArchiveError`（``ValueError``
    子类，消息含 ``tar -xzOf`` 指引）。
    """
    return compression.read_text_lines(path, encoding="utf-8-sig")


# 向后兼容别名（历史内部名）
_read_text_lines = read_text_lines


def strip_inline_comment(line: str) -> str:
    """去掉行内 ``#`` 之后的注释（``#`` 前须为空白或行首，避免切断含 ``#`` 的类群名）。"""
    if "#" not in line:
        return line
    idx = line.find("#")
    while idx >= 0:
        if idx == 0 or line[idx - 1] in " \t":
            return line[:idx]
        idx = line.find("#", idx + 1)
    return line


def _is_data_line(line: str) -> bool:
    """判断一行是否为约束数据行（与 ``parse_constraints`` 的跳过规则一致）。"""
    stripped = line.strip()
    if not stripped:
        return False
    if stripped.startswith("#"):
        return False
    # 锚定到首个字段，而不是整行子串（否则类群名含 "FRQ" 的数据行
    # 会被静默判定为非数据行而丢弃）。与 ``constraints/parsers.py`` 同一规则。
    first_field = stripped.replace(",", " ", 1).split(None, 1)[0]
    if first_field == "FRQ":
        return False
    return True


def _to_float(token: str, path: str, lineno: int, field: str, raw: str) -> float:
    """把 ``token`` 转成 float，失败/非法时给出**含文件名与行号**的清晰错误。"""
    try:
        value = float(token)
    except (TypeError, ValueError):
        raise ValueError(
            f"约束文件 {path} 第 {lineno} 行的{field}无法解析为数值：{token!r}\n"
            f"  原始行：{raw!r}\n"
            "  说明：数值列只接受十进制实数（不接受 '-'、'NA'、空串等记号）。"
        ) from None
    if not math.isfinite(value):
        raise ValueError(
            f"约束文件 {path} 第 {lineno} 行的{field}不是有限实数：{token!r}\n"
            f"  原始行：{raw!r}\n"
            "  说明：'nan' / 'inf' / '-inf' 会污染目标函数与百分比统计，已拒绝。"
        )
    return value


def validate_constraint_lines(lines: List[str], path: str = "<memory>") -> None:
    """校验约束文本行：权重/距离为有限实数、权重非负、端点非空。

    Args:
        lines: 原始行（可含换行符）。
        path: 用于错误信息的文件名。

    Raises:
        ValueError: 任一行含非法数值/空端点，信息含文件名与行号。
    """
    for lineno, raw in enumerate(lines, start=1):
        line = strip_inline_comment(raw.rstrip("\n"))
        if not _is_data_line(line):
            continue
        if "," in line:
            parts = line.split(",")
            if len(parts) < 4:
                raise ValueError(
                    f"约束文件 {path} 第 {lineno} 行的逗号格式列数不足：{raw!r}\n"
                    "  逗号格式必须包含 gene_family 首列，即 "
                    "'family,donor,receptor[,weight[,distance]]'（至少 4 列）；"
                    "若确为空格格式，请改用空格分隔 'donor receptor [weight]'。"
                )
            words = parts[1:]  # 逗号格式：丢弃 gene_family 首列
        else:
            words = line.split()
        if len(words) < 2:
            raise ValueError(
                f"约束文件 {path} 第 {lineno} 行缺少 donor/receptor 两端点：{raw!r}\n"
                "  说明：每行至少需要 'donor receptor'（权重缺省视为 1.0）。"
            )
        if not words[0].strip() or not words[1].strip():
            raise ValueError(f"约束文件 {path} 第 {lineno} 行的端点为空：{raw!r}")
        weight_token = words[2] if len(words) > 2 else "1.0"
        weight = _to_float(weight_token, path, lineno, "约束权重", raw)
        if weight < 0:
            raise ValueError(
                f"约束文件 {path} 第 {lineno} 行的约束权重为负：{weight_token!r}\n"
                f"  原始行：{raw!r}\n"
                "  说明：权重是支持度/概率/计数，必须 >= 0；负权重会使'被违反权重和'"
                "这一目标函数与所有百分比统计失去意义。"
            )
        if "," in line:
            distance_token = words[3] if len(words) > 3 else ""
            if distance_token.strip():
                _to_float(distance_token, path, lineno, "phylogenetic distance", raw)
        elif len(words) > 3:
            _to_float(words[3], path, lineno, "phylogenetic distance", raw)


def validate_constraint_set(cset: ConstraintSet, origin: str = "constraints") -> None:
    """校验一个 ``ConstraintSet``（API/适配器路径）的权重与距离合法性。

    与 :func:`validate_constraint_lines` 同一红线（有限、非负），但不依赖文件行号，
    以 ``donor,receptor`` 边键定位问题约束。
    """
    for c in cset.constraints:
        if not math.isfinite(c.weight) or c.weight < 0:
            raise ValueError(
                f"约束 {c.donor!r}->{c.receptor!r}（来源：{origin}）的权重非法："
                f"{c.weight!r}。权重必须是有限实数且 >= 0。"
            )
        dist = c.metadata.get("distance") if c.metadata else None
        if dist is not None and not math.isfinite(float(dist)):
            raise ValueError(
                f"约束 {c.donor!r}->{c.receptor!r}（来源：{origin}）的"
                f" phylogenetic distance 非法：{dist!r}。"
            )


def read_newick_file(path: str) -> Tree:
    """读取 Newick 文件并解析为 ``Tree``。

    原版为 ``readTree(open(...).readline())``（只看首行）。这是**已记录的偏离**
    ：折行书写的 Newick 在真实数据中常见，故此处把整份文件拼接为一条
    序列后再解析；空文件给出清晰错误而非 ``IndexError``。
    """
    lines = _read_text_lines(path)
    seq = "".join(line.strip() for line in lines).strip()
    if not seq:
        raise ValueError(
            f"物种树文件为空或仅含空白：{path}\n"
            "  说明：需要一个 Newick 字符串（如 '((A:1,B:1)1:1,C:1)2;'）。"
        )
    tree = Tree()
    try:
        tree.read_newick(seq)
    except ValueError as exc:
        raise ValueError(f"物种树文件 {path} 解析失败：{exc}") from None
    return tree


def read_constraints_file(path: str) -> ConstraintSet:
    """读取约束文件并解析为 ``ConstraintSet``。

    注：按距离过滤由 ``ConstraintSet.filter_by_distance`` 在排序阶段执行，
    不在此处发生；解析层不接收任何过滤参数。
    """
    lines = [_strip_inline(ln) for ln in _read_text_lines(path)]
    validate_constraint_lines(lines, path)
    try:
        return parse_constraints(lines)
    except ValueError as exc:
        raise ValueError(f"约束文件 {path} 解析失败：{exc}") from None


def _strip_inline(line: str) -> str:
    """对单行应用行内注释剥离并保留换行符。"""
    stripped = strip_inline_comment(line.rstrip("\n"))
    return stripped + ("\n" if line.endswith("\n") else "")
