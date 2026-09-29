#!/usr/bin/env python3
"""MaxTiC-Next 基准测试脚本。"""

import time
import sys
import os
import random
import statistics
import tempfile

sys.path.insert(0, "src")
from maxtic_next import rank
from maxtic_next.tree.tree import Tree
from maxtic_next.constraints.adapters import registry


def bench(func, n=5, **kwargs):
    times = []
    for _ in range(n):
        t0 = time.perf_counter()
        func(**kwargs)
        t1 = time.perf_counter()
        times.append(t1 - t0)
    return times


print("=" * 70)
print("MaxTiC-Next Benchmark Suite v0.1.0")
print("=" * 70)
print(f"Python {sys.version.split()[0]}")
print()

# ---- Benchmark 1: Core ranking (Cyanobacteria real data) ----
print("=== B1: Core ranking (Cyanobacteria, 13 internal nodes, 123 constraints) ===")
times = bench(
    rank,
    n=5,
    species_tree_path="examples/minitree.tree",
    constraints_path="examples/Cyano_CUTConstraints.tsv",
    seed=42,
    html_report=False,
    force=True,
    print_summary=False,
)
med = statistics.median(times)
print(f" Runs: {[f'{t:.3f}s' for t in times]}")
print(
    f" Median: {med:.3f}s | Mean: {statistics.mean(times):.3f}s | "
    f"Stdev: {statistics.stdev(times):.3f}s"
)
r = rank(
    "examples/minitree.tree",
    "examples/Cyano_CUTConstraints.tsv",
    seed=42,
    html_report=False,
    force=True,
    print_summary=False,
)
print(f" Result: best_source={r.best_source}, similarity={r.similarity_to_input}")
print()

# ---- Benchmark 2: Incremental vs full scoring ----
print("=== B2: Incremental vs Full scoring (local search 500 iters) ===")
times_inc = bench(
    rank,
    n=3,
    species_tree_path="examples/minitree.tree",
    constraints_path="examples/Cyano_CUTConstraints.tsv",
    seed=42,
    local_search=5,
    local_search_max_iterations=500,
    incremental=True,
    html_report=False,
    force=True,
    print_summary=False,
)
times_full = bench(
    rank,
    n=3,
    species_tree_path="examples/minitree.tree",
    constraints_path="examples/Cyano_CUTConstraints.tsv",
    seed=42,
    local_search=5,
    local_search_max_iterations=500,
    incremental=False,
    html_report=False,
    force=True,
    print_summary=False,
)
med_inc = statistics.median(times_inc)
med_full = statistics.median(times_full)
print(f" Incremental: median={med_inc:.3f}s, runs={[f'{t:.3f}s' for t in times_inc]}")
print(f" Full: median={med_full:.3f}s, runs={[f'{t:.3f}s' for t in times_full]}")
print(f" Speedup: {med_full / med_inc:.1f}x")
print()

# ---- Benchmark 3: Adapter parsing speed ----
print("=== B3: Adapter parsing speed ===")
t = Tree()
with open("examples/minitree.tree") as f:
    t.read_newick(f.readline())

adapters = [
    ("ALE", "ale", ["examples/adapters/ale/rec.uml_rec"]),
    ("RANGER-DTLx", "ranger", ["examples/adapters/ranger/FAM1.dtl"]),
    ("ecceTERA", "eccetera", ["examples/adapters/eccetera/FAM1.recphyloxml"]),
    ("ARTra", "artra", ["examples/adapters/artra/FAM1.txt"]),
    ("AleRax", "alerax", ["examples/adapters/alerax/run"]),
]
for name, tool, inputs in adapters:
    times_ad = []
    for _ in range(5):
        t0 = time.perf_counter()
        registry.convert(tool, t, inputs, min_family_size=0, quiet=True)
        t1 = time.perf_counter()
        times_ad.append(t1 - t0)
    med_ad = statistics.median(times_ad)
    print(f" {name:15s}: {med_ad * 1000:.2f}ms (median of 5)")
print()

# ---- Benchmark 4: Two-stage vs one-shot ----
print("=== B4: Two-stage vs One-shot ===")
tmpdir = tempfile.mkdtemp()
cons_out = os.path.join(tmpdir, "cons.tsv")

t0 = time.perf_counter()
rank(
    "examples/minitree.tree",
    "examples/adapters/ranger/FAM1.dtl",
    from_tool="ranger",
    ale_min_family_size=0,
    constraints_out=cons_out,
    html_report=False,
    force=True,
    print_summary=False,
    adapter_quiet=True,
)
t1 = time.perf_counter()
stage1 = t1 - t0

t0 = time.perf_counter()
rank(
    "examples/minitree.tree", cons_out, seed=42, html_report=False, force=True, print_summary=False
)
t1 = time.perf_counter()
stage2 = t1 - t0

t0 = time.perf_counter()
rank(
    "examples/minitree.tree",
    "examples/adapters/ranger/FAM1.dtl",
    from_tool="ranger",
    ale_min_family_size=0,
    seed=42,
    html_report=False,
    force=True,
    print_summary=False,
    adapter_quiet=True,
)
t1 = time.perf_counter()
oneshot = t1 - t0

print(
    f" Two-stage: stage1(adapter)={stage1 * 1000:.1f}ms + "
    f"stage2(rank)={stage2 * 1000:.1f}ms = {(stage1 + stage2) * 1000:.1f}ms"
)
print(f" One-shot: {oneshot * 1000:.1f}ms")
print(f" Overhead: {(oneshot - stage1 - stage2) * 1000:.1f}ms")
print()

# ---- Benchmark 5: Scalability (synthetic data) ----
print("=== B5: Scalability with synthetic data ===")
random.seed(42)


def make_tree(n_leaves):
    labels = [f"L{i}" for i in range(n_leaves)]
    random.shuffle(labels)
    internal_id = 0
    while len(labels) > 1:
        a = labels.pop()
        b = labels.pop()
        internal_id += 1
        labels.append(f"({a}:1,{b}:1){internal_id}")
    return labels[0] + ";"


def make_constraints(nwk, n_constraints):
    import re

    labels = re.findall(r"[A-Za-z_]\w*", nwk)
    labels = [x for x in set(labels) if not x.isdigit()]
    if len(labels) < 2:
        return ""
    lines = []
    for _ in range(n_constraints):
        d, r = random.sample(labels, 2)
        w = random.uniform(0.1, 10.0)
        lines.append(f"{d} {r} {w:.2f}")
    return "\n".join(lines)


for n_nodes in [10, 20, 50, 100]:
    n_leaves = n_nodes + 1
    nwk = make_tree(n_leaves)
    n_cons = n_nodes * 30
    cons = make_constraints(nwk, n_cons)

    tmp_tree = os.path.join(tmpdir, f"tree_{n_nodes}.nwk")
    tmp_cons = os.path.join(tmpdir, f"cons_{n_nodes}.tsv")
    with open(tmp_tree, "w") as f:
        f.write(nwk)
    with open(tmp_cons, "w") as f:
        f.write(cons)

    times_sc = []
    for _ in range(3):
        t0 = time.perf_counter()
        rank(tmp_tree, tmp_cons, seed=42, html_report=False, force=True, print_summary=False)
        t1 = time.perf_counter()
        times_sc.append(t1 - t0)
    med_sc = statistics.median(times_sc)
    print(f" n={n_nodes:3d} nodes, {n_cons:4d} constraints: {med_sc * 1000:.1f}ms")

print()
# ---- Benchmark 6: Reproducibility ----
print("=== B6: Reproducibility (fixed seed) ===")
results = []
for _ in range(5):
    r = rank(
        "examples/minitree.tree",
        "examples/Cyano_CUTConstraints.tsv",
        seed=42,
        html_report=False,
        force=True,
        print_summary=False,
    )
    results.append(r.best_order)
all_same = all(rr == results[0] for rr in results)
print(f" 5 runs with seed=42: {'IDENTICAL' if all_same else 'DIFFERENT'}")
print()

print("=" * 70)
print("Benchmark complete.")
print("=" * 70)
