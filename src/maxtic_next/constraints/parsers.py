"""约束文件双格式解析。

支持两种输入格式：

1. 空格格式：``donor receptor [weight] [distance]``（``weight`` 与第 4 列
   ``distance`` 均可选；缺失 ``weight`` 默认 ``1.0``，缺失 ``distance`` 视为
   无距离列 ``None``）。``distance`` 供 ``--min-transfer-distance`` 过滤，
   与原版 ``MaxTiC.py:411`` 的 ``len(words)<=3 or float(words[3])>MIN`` 一致。
2. 逗号格式：``family,donor,receptor,[weight]``（**首列 gene_family 被丢弃**，
   ``weight`` 可选）；以及 5 列 ALE 输出 ``family,donor,receptor,weight,distance``，
   其中第 5 列 ``distance`` 存入 ``metadata`` 供 ``--min-transfer-distance`` 过滤。

最高风险约束：**禁止以三列 ``donor,receptor,weight`` 作为逗号输入**，
因为首列丢弃后会把 donor/receptor 与 weight 混淆。故逗号格式少于 4 列（即 2 列或
3 列逗号）一律报错，给出清晰信息而非 ``IndexError``。
"""

from typing import Iterable

from maxtic_next.constraints.constraint import Constraint, ConstraintSet


def parse_constraints(lines: Iterable[str]) -> ConstraintSet:
    """解析约束文件内容为 ``ConstraintSet``。

    Args:
        lines: 约束文件逐行内容（可迭代的字符串）。

    Returns:
        解析得到的 ``ConstraintSet``。

    Raises:
        ValueError: 遇到三列逗号输入（``donor,receptor,weight`` 误用）等不可解析行。
    """
    cset = ConstraintSet()
    for raw in lines:
        line = raw.rstrip("\n")
        stripped = line.strip()
        if not stripped:
            continue
        # 跳过注释与 ALE 的 FRQ 汇总行。FRQ 汇总行的真实形态是"首字段为 FRQ"，
        # 故锚定到首个字段（空格或逗号分隔）——按整行子串匹配会把
        # **类群名恰含 FRQ 子串**的数据行整行静默删除。
        if stripped.startswith("#"):
            continue
        if stripped.replace(",", " ", 1).split(None, 1)[0] == "FRQ":
            continue

        if "," in line:
            parts = line.split(",")
            # 少于 4 列的逗号输入（2 列或 3 列）= family,donor,receptor 结构不完整
            # 或 donor,receptor,weight 误用，必须报错
            if len(parts) < 4:
                raise ValueError(
                    "逗号格式必须包含 gene_family 首列，即 "
                    "'family,donor,receptor[,weight[,distance]]'（至少 4 列）。\n"
                    f"  收到 {len(parts)} 列逗号输入 {parts!r}，首列丢弃规则会"
                    "误丢 donor/receptor 或与其余列混淆，已拒绝。\n"
                    "  若确为空格格式，请改用空格分隔：'donor receptor [weight]'。"
                )
            family = parts[0]
            words = parts[1:]  # 丢弃首列（与原版 words = words[1:] 一致）
            # words: [donor, receptor, (weight), (distance)]
            if len(words) > 3:
                # 5 列 ALE：family,donor,receptor,weight,distance
                weight = float(words[2])
                # 距离列允许为空（适配器写出 distance=None 时的空字段），视同无距离列
                distance = float(words[3]) if words[3].strip() else None
            elif len(words) > 2:
                # 4 列：family,donor,receptor,weight
                weight = float(words[2])
                distance = None
            else:
                # 2 列：family,donor,receptor（无权重）
                weight = 1.0
                distance = None
            donor, receptor = words[0], words[1]
            cset.add(
                Constraint(
                    donor=donor,
                    receptor=receptor,
                    weight=weight,
                    metadata={"family": family, "support": None, "distance": distance},
                )
            )
        else:
            words = line.split()
            if len(words) < 2:
                continue
            donor, receptor = words[0], words[1]
            weight = float(words[2]) if len(words) > 2 else 1.0
            # 与原版 MaxTiC.py:411 一致：空格格式第 4 列为 phylogenetic
            # distance，供 --min-transfer-distance（filter_by_distance）过滤使用；
            # 列数 <= 3 时视为"无距离列"（None，过滤时一律保留）。
            distance = float(words[3]) if len(words) > 3 else None
            cset.add(
                Constraint(
                    donor=donor,
                    receptor=receptor,
                    weight=weight,
                    metadata={"family": None, "support": None, "distance": distance},
                )
            )
    return cset


def parse_constraints_file(path: str) -> ConstraintSet:
    """从文件读取并解析约束（便捷封装）。

    ：读取走 :mod:`maxtic_next.io.compression` 的共享文本读取器，因此
    ``.gz`` / gzip 流与单成员 ``.tar.gz`` 归档可直接作为输入。此处在函数内导入
    （而非模块级），以避免 ``io.parsing ↔ constraints.parsers`` 的循环导入。
    """
    from maxtic_next.io.compression import read_text_lines

    return parse_constraints(read_text_lines(path, encoding="utf-8"))
