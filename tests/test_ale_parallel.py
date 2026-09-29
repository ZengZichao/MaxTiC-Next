# Copyright (C) 2026 MaxTiC-Next rewrite team.
# This file is part of MaxTiC-Next, a Python 3 rewrite of MaxTiC.
#
# MaxTiC is Copyright (C) Eric Tannier (Inria / CNRS / ENS Lyon).
# Original reference: MaxTiC: Fast ranking of a phylogenetic tree by
# Maximum Time Consistency with lateral gene transfers,
# Biorxiv doi.org/10.1101/127548
#
# MaxTiC-Next is free software: you can redistribute it and/or modify it
# under the terms of the CeCILL 2.1 license as published by the
# Commissariat a l'energie atomique et aux energies alternatives (CEA),
# the Centre National de la Recherche Scientifique (CNRS) and the
# Institut National de Recherche en Informatique et en Automatique (INRIA).
# A copy of the license is included in the LICENSE file.

"""ALE 适配器并行化回归测试。

验证 ``ALEAdapter.convert`` 的线程池并行解析与顺序执行**逐字节一致**：
由于 ``ConstraintSet.add`` 按边键聚合（顺序无关），并行合并后的
"边键 + 权重"集合必须与顺序执行完全相同。

比对待用顺序无关聚合：``ConstraintSet.informative_edges()`` 返回
``边键 -> 总权重`` 字典，对约束列表顺序不敏感 —— 这正是
"顺序无关合并 + 聚合后比对" 的严格判定方式。

注意：测试统一使用 ``cache_dir=None`` 强制重新解析（不命中文件级缓存），
从而真正覆盖"并行 recompute -> 顺序无关合并 = 顺序 recompute" 这条路径。
"""

import os

from maxtic_next.constraints.adapters.ale import convert_from_ale
from maxtic_next.tree.tree import Tree

DATA_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data", "ale_example")
SPECIES = os.path.join(DATA_DIR, "species.tree")
REC = os.path.join(DATA_DIR, "rec.uml_rec")


def _read_tree(path):
    """读入 Newick 物种树。"""
    t = Tree()
    with open(path) as fh:
        t.read_newick(fh.readline())
    return t


def _aggregated(cset):
    """边键 -> 总权重（顺序无关），用于顺序/并行一致性比对。"""
    return sorted((key, round(weight, 9)) for key, weight in cset.informative_edges().items())


def _write_synthetic_rec(path: str, newick: str, family: str = "FAM") -> None:
    """写出一个合规的 ALE .uml_rec 片段（1 个调和家族）。"""
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("1 reconciled\n")
        fh.write(family + "\n")
        fh.write(newick + "\n")


# ---------------------------------------------------------------------------
# 真实单文件：顺序 vs 并行
# ---------------------------------------------------------------------------
def test_sequential_vs_parallel_identical_on_real_file():
    """真实 rec.uml_rec：max_workers=1（顺序执行器）与 max_workers=8 结果一致。"""
    species = _read_tree(SPECIES)
    sequential = convert_from_ale(species, [REC], min_family_size=0, cache_dir=None, max_workers=1)
    parallel = convert_from_ale(species, [REC], min_family_size=0, cache_dir=None, max_workers=8)
    assert _aggregated(sequential) == _aggregated(parallel)
    # 真实样例应至少产出一条约束（测试 test_ale_adapter 已验证 (65,59)）
    assert len(parallel) >= 1


# ---------------------------------------------------------------------------
# 多合成文件：顺序 vs 并行（覆盖 ThreadPoolExecutor 聚合路径）
# ---------------------------------------------------------------------------
def test_sequential_vs_parallel_identical_on_multiple_files(tmp_path):
    """三个互不相同合成家族：顺序 vs 并行，边键+权重集合完全一致。"""
    species = _read_tree(SPECIES)
    f1 = str(tmp_path / "f1.uml_rec")
    f2 = str(tmp_path / "f2.uml_rec")
    f3 = str(tmp_path / "f3.uml_rec")
    # 三个家族的调和事件解析出不同的 donor/receptor，使合并具非空意义
    _write_synthetic_rec(f1, "((CYAP8:1.0,CYAP0:1.0)59.T@61->62:1.0,CYAA5:1.0)X65:1.0;", "FAM1")
    _write_synthetic_rec(f2, "((CYAP8:1.0,CYAA5:1.0)61.T@62->59:1.0,CYAP0:1.0)X65:1.0;", "FAM2")
    _write_synthetic_rec(f3, "((CYAP8:1.0,CYAA5:1.0)62.T@65->61:1.0,CYAP0:1.0)X65:1.0;", "FAM3")
    files = [f1, f2, f3]

    sequential = convert_from_ale(species, files, min_family_size=0, cache_dir=None, max_workers=1)
    parallel = convert_from_ale(species, files, min_family_size=0, cache_dir=None, max_workers=8)

    assert _aggregated(sequential) == _aggregated(parallel)
    assert len(parallel) >= 1


# ---------------------------------------------------------------------------
# 默认并行度（max_workers=None）与显式并行一致
# ---------------------------------------------------------------------------
def test_default_max_workers_matches_explicit_parallel(tmp_path):
    """max_workers=None（自动取 min(8, len, cpu)）应与显式 max_workers=8 一致。"""
    species = _read_tree(SPECIES)
    g1 = str(tmp_path / "g1.uml_rec")
    g2 = str(tmp_path / "g2.uml_rec")
    _write_synthetic_rec(g1, "((CYAP8:1.0,CYAP0:1.0)59.T@61->62:1.0,CYAA5:1.0)X65:1.0;", "FAM1")
    _write_synthetic_rec(g2, "((CYAP8:1.0,CYAA5:1.0)62.T@65->61:1.0,CYAP0:1.0)X65:1.0;", "FAM2")
    files = [g1, g2]

    auto = convert_from_ale(species, files, min_family_size=0, cache_dir=None, max_workers=None)
    parallel = convert_from_ale(species, files, min_family_size=0, cache_dir=None, max_workers=8)

    assert _aggregated(auto) == _aggregated(parallel)


# ---------------------------------------------------------------------------
# trf 来源同样满足顺序无关一致性
# ---------------------------------------------------------------------------
def test_parallel_consistency_for_transfer_source(tmp_path):
    """source='trf'（转移事件）下，顺序 vs 并行仍一致。"""
    species = _read_tree(SPECIES)
    h1 = str(tmp_path / "h1.uml_rec")
    h2 = str(tmp_path / "h2.uml_rec")
    _write_synthetic_rec(h1, "((CYAP8:1.0,CYAP0:1.0)59.T@61->62:1.0,CYAA5:1.0)X65:1.0;", "FAM1")
    _write_synthetic_rec(h2, "((CYAP8:1.0,CYAA5:1.0)61.T@62->59:1.0,CYAP0:1.0)X65:1.0;", "FAM2")
    files = [h1, h2]

    sequential = convert_from_ale(
        species, files, min_family_size=0, source="trf", cache_dir=None, max_workers=1
    )
    parallel = convert_from_ale(
        species, files, min_family_size=0, source="trf", cache_dir=None, max_workers=8
    )

    assert _aggregated(sequential) == _aggregated(parallel)
