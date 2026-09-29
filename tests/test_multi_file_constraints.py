# Copyright (C) 2026 MaxTiC-Next rewrite team.
# MaxTiC-Next is free software under the CeCILL 2.1 license (see LICENSE).

"""多约束文件合并与 -o 模式守卫的回归测试。

覆盖：

1. 文本模式下传入多个约束文件时，全部文件必须被合并计权（手册 03 与
   CLI synopsis "CONSTRAINTS [...]" 的约定），任何文件不得被静默忽略；
2. ``-o`` / ``--constraints-out`` 仅在适配器模式下有效，文本模式下
   传入应显式报错而非静默忽略。
"""

import os

import pytest

from maxtic_next.api import rank

# 路径以本文件为锚点，而不是以进程 CWD 为锚点：从仓库根目录之外运行
# ``pytest tests/`` 时，相对路径写法会给出误导性的 FileNotFoundError。
DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")
TREE = os.path.join(DATA, "minitree.tree")
CONS = os.path.join(DATA, "Cyano_CUTConstraints.tsv")


def test_multiple_text_constraint_files_are_merged(tmp_path):
    f1 = tmp_path / "a.tsv"
    f2 = tmp_path / "b.tsv"
    f1.write_text("61 62 10.0\n", encoding="utf-8")
    f2.write_text("62 45 30.0\n", encoding="utf-8")

    r_both = rank(str(TREE), [str(f1), str(f2)], seed=42, print_summary=False)
    r_one = rank(str(TREE), [str(f1)], seed=42, print_summary=False)
    # 两个文件都被计入：总权重 = 10 + 30；仅第一文件 = 10
    assert r_both.total_weight == pytest.approx(40.0)
    assert r_one.total_weight == pytest.approx(10.0)

    # 与"合并为单文件"的结果一致（同一组约束、同一 seed）
    f_merged = tmp_path / "merged.tsv"
    f_merged.write_text(
        f1.read_text(encoding="utf-8") + f2.read_text(encoding="utf-8"), encoding="utf-8"
    )
    r_merged = rank(str(TREE), str(f_merged), seed=42, print_summary=False)
    assert r_merged.total_weight == pytest.approx(r_both.total_weight)
    assert r_merged.best_order == r_both.best_order
    assert r_merged.values == r_both.values


def test_constraints_out_rejected_for_text_mode():
    with pytest.raises(ValueError, match="constraints_out"):
        rank(str(TREE), str(CONS), seed=42, print_summary=False, constraints_out="out.tsv")
