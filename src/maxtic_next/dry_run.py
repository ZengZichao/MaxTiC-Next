"""预检与断点续传支持。

``dry_run_check`` 在正式排序**之前**对输入做静态校验，覆盖以下检查项：

* 物种树格式（可解析、存在根）；
* **二叉性**（多歧分支为 error：本项目不支持非二叉物种树）；
* 内部节点标签唯一性 + **叶子名唯一性**；
* 约束文件格式（双格式可解析；拒绝三列逗号误用；4/5 列数值合法）；
* **约束权重必须有限且 >= 0**（``nan`` / ``inf`` / 负值均为 error）；
* 约束端点合法性（donor / receptor 须为物种树中存在的标签）；
* ``d`` 距离列语义（``--min-transfer-distance > 0`` 且无距离列时，该选项**被忽略**
  而非"丢弃全部约束"）；
* ``--from-ale`` / 适配器路径：校验上游输出文件列表可读取。

返回结构化的 ``report``（含 ``issues`` 列表，每项有 ``severity`` / ``message`），
供 CLI 打印并**据 ``report["ok"]`` 决定退出码**（不得以格式化文本的子串
判断成败）。``ok`` 为是否存在致命（error）问题。

注意：本预检**不修改**任何数据，也不消耗随机数，可安全在管线最前端插入。
"""

import os
from typing import Dict, List, Optional

from maxtic_next.constraints.parsers import parse_constraints
from maxtic_next.io.compression import open_text
from maxtic_next.io.parsing import read_text_lines, validate_constraint_lines
from maxtic_next.tree.tree import Tree


def _check_species_tree(tree: Tree, issues: List[Dict]) -> None:
    """校验物种树：可解析、有根、内部标签唯一、二叉、叶子名唯一。"""
    root = tree.get_root()
    if root < 0:
        issues.append({"severity": "error", "message": "物种树无根（解析失败或为空）"})
        return
    internal = tree.internal_node_labels()
    if len(internal) != len(set(internal)):
        dup = [x for x in internal if internal.count(x) > 1]
        issues.append(
            {
                "severity": "error",
                "message": f"物种树内部节点标签不唯一（重复标签：{sorted(map(str, set(dup)))}），"
                "会导致端点映射歧义，必须修正。",
            }
        )
    if not internal:
        issues.append({"severity": "warning", "message": "物种树没有内部节点，排序无意义。"})
    # ：重复叶子名会让节点标识塌缩
    dup_leaves = tree.duplicate_leaf_names()
    if dup_leaves:
        detail = "; ".join(
            f"{name!r} 出现在父节点 {parents}" for name, parents in sorted(dup_leaves.items())[:5]
        )
        issues.append(
            {
                "severity": "error",
                "message": f"物种树存在 {len(dup_leaves)} 个重复叶子名：{detail}。"
                "重复类群名会使排序丢失内部节点并在建树时崩溃"
                "（常见成因：误把基因树/带样本编号的拷贝当作物种树）。",
            }
        )


def _check_binary_tree(tree: Tree, issues: List[Dict]) -> None:
    """校验物种树满足二叉树假设：每个内部节点恰有 2 个子节点。

    原版 ``opt`` / ``mix`` 强假设二叉树；多歧分支（polytomy，内部节点子节点数 != 2）
    会被静默丢弃多余子节点，导致结果不可预期。预检阶段显式报错，避免正式运行
    产生误导性的"最优序"。正式运行路径另有无条件守卫（``Ranker.run()``）。
    """
    for label, n_children in tree.non_binary_internal_nodes():
        issues.append(
            {
                "severity": "error",
                "message": f"物种树存在多歧分支（polytomy）：内部节点 {label!r} "
                f"有 {n_children} 个子节点，而算法强假设二叉树"
                f"（opt/mix 会静默丢弃多余子节点，结果不可靠）。"
                f"请解消多歧，或重新定根使每个内部节点恰有两个子节点"
                f"（未定根 Newick 的根天然三出）。",
            }
        )


def _check_target_clade(tree: Tree, target_clade: Optional[str], issues: List[Dict]) -> None:
    """校验 ``--target-clade`` 标签合法性：若提供且不在物种树中，预检阶段即报错。

    否则正式运行时 ``prune_constraints`` 会抛出 ``KeyError``，错误位置不清晰。
    """
    if target_clade is None:
        return
    try:
        tree.label_to_node_id(target_clade)
    except KeyError:
        issues.append(
            {
                "severity": "error",
                "message": f"--target-clade 指定的节点 {target_clade!r} 不在物种树中"
                f"（预检阶段拦截，避免正式运行 KeyError）。"
                f"请检查标签是否为物种树内部节点的 bootstrap 标签。",
            }
        )


def _check_constraints_file(
    path: str, tree: Tree, min_transfer_distance: float, issues: List[Dict]
) -> None:
    """校验单个约束文件（文本双格式）：格式、权重合法性、端点、d 列语义。

    读取与数值校验统一走 ``io.parsing`` 的helper：带 BOM、
    行内 ``#`` 注释、``nan`` / ``inf`` / 负权重等都在**此层**给出**含文件名与行号**的
    信息，与正式运行的解析层保持同一判据。
    """
    if not os.path.exists(path):
        issues.append({"severity": "error", "message": f"约束文件不存在：{path}"})
        return
    try:
        lines = read_text_lines(path)
    except Exception as exc:  # noqa: BLE001 - 预检需捕获一切读取异常
        issues.append({"severity": "error", "message": f"约束文件 {path} 读取异常：{exc}"})
        return
    try:
        validate_constraint_lines(lines, path)
    except ValueError as exc:
        issues.append({"severity": "error", "message": str(exc)})
        return
    try:
        cset = parse_constraints(lines)
    except ValueError as exc:
        issues.append({"severity": "error", "message": f"约束文件 {path} 解析失败：{exc}"})
        return
    except Exception as exc:  # noqa: BLE001 - 预检需捕获一切解析异常
        issues.append({"severity": "error", "message": f"约束文件 {path} 读取异常：{exc}"})
        return

    valid_labels = set(str(x) for x in tree.internal_node_labels()) | set(tree.get_leaves_names())
    bad_endpoints = 0
    n_with_distance = 0
    for c in cset.constraints:
        if c.donor not in valid_labels or c.receptor not in valid_labels:
            bad_endpoints += 1
        if c.metadata.get("distance") is not None:
            n_with_distance += 1
    if bad_endpoints > 0:
        issues.append(
            {
                "severity": "error",
                "message": f"约束文件 {path} 有 {bad_endpoints} 条约束的端点不在物种树中"
                "（donor/receptor 必须是内部节点 bootstrap 或叶子名）。",
            }
        )
    # d 距离列语义校验
    # 事实：filter_by_distance 对**无距离列**的约束一律保留（与原版 MaxTiC.py:411 一致），
    # 因此该选项是"被忽略"，而不是旧文案所说的"将丢弃全部约束"。
    if min_transfer_distance > 0 and len(cset.constraints) > 0 and n_with_distance == 0:
        issues.append(
            {
                "severity": "error",
                "message": f"--min-transfer-distance={min_transfer_distance} 将被忽略："
                f"约束文件 {path} 不含 phylogenetic distance 列（4 列空格 / 5 列 ALE 格式），"
                "无距离列的约束一律保留，即本次运行不会按距离删除任何约束。"
                "如确需按距离过滤，请补充距离列；如不需要，请去掉该选项。",
            }
        )


def _check_ale_files(rec_files: List[str], issues: List[Dict]) -> None:
    """校验 ALE / 适配器路径：每个上游输出文件可读取。

    ：读取走 :func:`maxtic_next.io.compression.open_text`，故 ``.gz`` /
    gzip 流上游输出在预检里同样透明解压；多成员 ``.tgz`` 归档给出的**不是**裸
    traceback，而是含 ``tar -xzf`` / ``tar -xzOf`` 命令的 error 级条目。
    """
    for f in rec_files:
        if not os.path.exists(f):
            issues.append(
                {"severity": "error", "message": f"上游输出文件（.uml_rec 等）不存在：{f}"}
            )
            continue
        try:
            with open_text(f, encoding="utf-8-sig") as fh:
                fh.readline()
        except Exception as exc:  # noqa: BLE001
            issues.append(
                {"severity": "error", "message": f"上游输出文件（.uml_rec 等）{f} 读取失败：{exc}"}
            )


def dry_run_check(
    tree: Tree,
    constraints_paths: List[str],
    from_ale: bool = False,
    min_transfer_distance: float = 0.0,
    target_clade: Optional[str] = None,
) -> Dict:
    """对输入做静态预检，返回结构化报告。

    Args:
        tree: 已解析的物种树。
        constraints_paths: 约束文件路径列表（文本模式）或上游输出文件列表（适配器模式）。
        from_ale: 是否为适配器模式（校验上游文件可读性而非文本约束格式）。
        min_transfer_distance: 最小转移距离阈值（用于 ``d`` 列语义校验）。
        target_clade: ``--target-clade`` 目标类群根节点标签；若提供但不在物种树中，
            预检阶段即报错（避免正式运行 ``KeyError``）。

    Returns:
        dict，含 ``issues``（列表，每项 ``{severity, message}``）、``ok``（无 error）、
        ``n_errors``、``n_warnings``。``ok`` 是 CLI 退出码的**唯一**依据。
    """
    issues: List[Dict] = []
    _check_species_tree(tree, issues)
    _check_binary_tree(tree, issues)
    _check_target_clade(tree, target_clade, issues)
    if from_ale:
        _check_ale_files(constraints_paths, issues)
    else:
        for p in constraints_paths:
            _check_constraints_file(p, tree, min_transfer_distance, issues)

    n_errors = sum(1 for it in issues if it["severity"] == "error")
    n_warnings = sum(1 for it in issues if it["severity"] == "warning")
    return {
        "issues": issues,
        "ok": n_errors == 0,
        "n_errors": n_errors,
        "n_warnings": n_warnings,
    }


def format_dry_run(report: Dict) -> str:
    """把预检报告格式化为可读字符串。

    注意：本函数**只用于人读**。调用方（CLI / GUI）判定成败必须使用
    ``report["ok"]``，不得对本字符串做子串匹配——旧实现以 ``"错误" in report_str``
    判据，而结论行恒含"错误 N"，导致合法输入也返回退出码 1。
    """
    lines = ["=== 预检（dry-run）报告 ==="]
    if not report["issues"]:
        lines.append("未发现任何问题。")
    for it in report["issues"]:
        prefix = "error" if it["severity"] == "error" else "warning"
        lines.append(f"[{prefix}] {it['message']}")
    lines.append(
        f"结论：{'通过' if report['ok'] else '存在致命问题，已中止'} "
        f"（error {report['n_errors']}, warning {report['n_warnings']}）"
    )
    return "\n".join(lines)
