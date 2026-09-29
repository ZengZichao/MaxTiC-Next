# 00 · Overview

## 0.1 What is MaxTiC-Next

`MaxTiC-Next` is a **faithful Python 3 rewrite and engineering extension** of the
Python 2 phylogenetics tool **MaxTiC** (Eric Tannier, Inria). It solves one well-defined
problem:

> Given a **rooted species tree** and a set of **weighted time constraints** inferred from
> lateral/horizontal gene transfers (HGT) — each of the form "node A is not younger than
> node B" — find an **ordering (relative chronology)** of the tree's internal nodes that
> **maximizes the cumulative weight of satisfied constraints** (equivalently, minimizes the
> weight of the Feedback Arc Set, FAS).

In one sentence: **use the "who is older than whom" signal from HGT to give the internal
nodes of a species tree a maximally self-consistent temporal order.**

## 0.2 Core concepts & glossary

| Term | Meaning |
|------|---------|
| **Species tree** | Rooted Newick tree with **internal node labels stored in the bootstrap field** (e.g. `...)61:2,...`, `61` is the label) |
| **Constraint** | Directed weighted edge `donor → receptor`, meaning "donor no later than receptor" |
| **Donor / receptor** | Species-tree nodes for the source / target branch of a transfer |
| **Weight** | Confidence score of a constraint (e.g. sampling support, posterior frequency) |
| **MTC** | Maximum Time Consistency: maximize weight of non-conflicting constraints |
| **FAS** | Feedback Arc Set: the violated constraints; MTC ≡ minimize FAS weight |
| **Informative constraint** | Neither a tree ancestor-descendant relation nor a trivial cycle; both endpoints internal |
| **Partial order** | Constraints consistent with the best order, forming a relative chronology |
| **Kendall similarity** | Normalized agreement between output order and the input tree order (if ultrametric), ∈[0,1] |

## 0.3 Three ranking heuristics

MaxTiC-Next runs three heuristics on the same problem and keeps the best:

1. **Greedy (`order_from_graph`)**: add constraints by descending weight while avoiding cycles, then topologically sort.
2. **Mixing (`opt`/`mix`)**: bottom-up dynamic programming that merges the in-orders of left/right subtrees.
3. **Local search (`optimisation_locale`)**: Metropolis interval-rotation perturbations on the current best order within a time budget (enabled by `--local-search`).

## 0.4 Capability snapshot

- **Algorithmic equivalence**: core algorithm ported function-by-function. **The default path
  (`--ls 0`) is byte-for-byte reproducible under a fixed `--seed`**; with `--ls > 0` the stop
  condition is wall-clock time, so the iteration count follows machine load and
  `--local-search-max-iters` is also needed for cross-machine reproducibility
  (see 03 §3.8).
- **Dual input formats**: space `donor receptor [weight] [distance]` and comma `family,donor,receptor,[weight],[distance]`.
- **Five upstream adapters**: ALE, RANGER-DTLx, ecceTERA, ARTra, AleRax — one-shot `--from <tool>` or `--from-auto` detection.
- **Interactive HTML report**: on by default; weight distribution, conflict ratio, node-position robustness, ranked-tree view.
- **Robustness / sensitivity analysis**: local search collects near-optimal orders; node-position and pairwise-order statistics.
- **MCMC sampler** (optional, **preliminary**): reversible Metropolis–Hastings sampling over
  topology-compatible linear extensions; convergence diagnostics are **not** validated, so its
  samples must **not** be treated as posterior samples (see 08 §8.2).
- **Engineering**: modern CLI + Python API, Docker/Singularity, Snakemake/Nextflow, dry-run
  preflight, pruning, checkpoint resume, incremental scoring, output overwrite protection
  (`--force`).

## 0.5 Source layout

```
MaxTiC-Next/                  # repository root (https://github.com/ZengZichao/MaxTiC-Next)
├── src/maxtic_next/           # the package (import name: maxtic_next)
│   ├── cli.py  api.py         # CLI / Python API entry points
│   ├── tree/                  # lightweight Tree class
│   ├── constraints/           # data model, dual-format parsing, 5 adapters
│   │   └── adapters/          # ale/ranger_dtl/eccetera/artra/alerax + registry
│   ├── ranking/               # 3 heuristics + main Ranker
│   ├── robustness/            # sensitivity summary + MCMC
│   ├── io/  report/           # 3-file output / HTML report
│   └── benchmarks/            # benchmark suite
├── docs/                      # this manual + adapter specs (the desktop guide lives in the Studio repo)
├── examples/                  # example data (minitree + Cyano + per-tool samples)
├── workflows/                 # Snakemake / Nextflow
├── tests/                     # test suite
├── Dockerfile  Singularity.def
└── pyproject.toml
```

## 0.6 Naming convention (three layers)

| Layer | Value |
|-------|-------|
| Product / display name | `MaxTiC-Next` |
| GitHub repository | [`ZengZichao/MaxTiC-Next`](https://github.com/ZengZichao/MaxTiC-Next) |
| Python import package | `maxtic_next` |
| Console command | `maxtic-next` |

> Next: [01 · Installation](01_installation.md)
