# 02 · Quickstart

This chapter runs one complete ranking in under five minutes on the dataset shipped
with the repository (the small cyanobacterial example).

> Every stdout listing and file list below is **captured from an actual run**
> (Python 3.14, local machine; copy the commands verbatim to reproduce them), and
> matches the program's current output line for line.

## 2.1 Input data

- Species tree: `examples/minitree.tree` (13 internal nodes, Newick, internal labels
  carried in the bootstrap field, **single line**)
- Constraints: `examples/Cyano_CUTConstraints.tsv` (space format `donor receptor weight`,
  **123 lines** — the cyanobacterial constraint set distributed with the original MaxTiC)

Species tree fragment (internal node labels are the numbers `42/59/61/...`):

```
(((((CYAP8:1,CYAP0:1)42:7,(CYAA5:7,UCYNA:7)41:1)59:2,...)69;
```

Constraint fragment (`donor receptor weight`, tab-separated):

```
61	62	15.04
61	67	15.48
61	46	28.73
```

## 2.2 Minimal run

```bash
maxtic-next examples/minitree.tree examples/Cyano_CUTConstraints.tsv --seed 42
```

Without installing (or if you prefer not to), run straight from the source tree at
the repository root:

```bash
PYTHONPATH=src python3 -m maxtic_next examples/minitree.tree \
    examples/Cyano_CUTConstraints.tsv --seed 42
```

The two are exactly equivalent (`maxtic-next` is the console entry point for
`maxtic_next.cli:main`).

## 2.3 Expected stdout summary

Below is the **actual output** of the command above from the repository root (the HTML
report is written at the same time; `--no-html` only removes `*.html`, stdout is
unchanged). Note the **two spaces** after `tree with` (matching the original `print`):

```
examples/Cyano_CUTConstraints.tsv
tree with  13 internal nodes
2218.4 total weight of constraints from transfers
307.45000000000005 uninformative (  307.45000000000005 to a descendant, 0.0 to a leaf 0.0 to an ancestor 0.0 to itself 13.9%)
676.0400000000001 trivially conflicting constraints (30%) (descendant to ancestor or trivial cycle)
value of the order given by the input tree 755.8800000000005 (34.0732059141724%)
value of the greedy heuristic 755.8800000000005 (34.0732059141724%)
value of the mixing heuristic: 755.8800000000005 (34.0732059141724%)
best order is the greedy heuristic
(((((CYAP8:1,CYAP0:1)42:7,(CYAA5:7,UCYNA:7)41:1)59:2,(MICAN:9,(CYAP2:3,CYAP7:3)43:6)45:1)61:2,(SYNY3:11,SYNP2:11)56:1)65:1,(TRIEI:6,(NOSA0:5,(NOSP7:4,(NOSS1:2,ANAVT:2)40:2)46:1)62:1)67:7)69;
Similarity of the best order compared with the input order 1.0
0.0 constraints in agreement with the best tree in conflict with input tree
Similarity of the order compared with the input order 1.0
NOTE: 12 条约束与物种树谱系边同向（权重合计 63.55），已被置为不可违反的拓扑硬约束，并已计入上方 uninformative 统计（在 total_weight 分母中只计一次，不重复计权）。
```

(Translation of the trailing NOTE, which the program emits in Chinese: *"12 constraints
agree with an edge of the species phylogeny (total weight 63.55); they have been turned
into inviolable topological hard constraints and are already counted in the
uninformative summary above (once in the `total_weight` denominator, not twice)."*)

**How to read this summary**:

| Line | Meaning |
|------|---------|
| `examples/Cyano_CUTConstraints.tsv` | Output prefix for this run = path of the **first** constraint file |
| `tree with  13 internal nodes` | Number of internal nodes to rank (two spaces, as in the original) |
| `2218.4 total weight of constraints from transfers` | Sum of all transfer-constraint weights (the denominator) |
| `uninformative (...)` | Weight classified as uninformative (to a descendant / leaf / ancestor / itself) and its share; **the percentage is `uninformative ÷ total_weight`**, which prints `13.9%` here — the original denominator double-counted the uninformative weight and printed `12%` on the same instance |
| `trivially conflicting` | Weight of trivially conflicting constraints (descendant→ancestor or trivial cycle) and its share |
| `value of the ... heuristic` | "Weight of violated constraints" for each of the three rankings (lower is better), plus its share of `total_weight` |
| `best order is the greedy heuristic` | The heuristic finally selected; becomes `... + local search` when `--ls > 0` actually improves things |
| `after local search X rejected` / `best found solution Y (%)` | **Only with `--ls > 0`** — two extra lines reporting the post-search value |
| one Newick line | The **ranked tree**: branch lengths set to the ranking positions |
| `Similarity ...` / `constraints in agreement ...` | Kendall similarity with the input ordering (∈[0,1], 1 = identical), etc. |
| `values from N random orders ...` (4 lines) | **Only with `--random-trees > 0`** (see 08 §8.7) |
| `NOTE: ...` / `WARNING: ...` | Runtime statements of convention and warnings, always appended at the **end** of the summary |

> The trailing NOTE has **no counterpart in the original** (the per-edge print at `MaxTiC.py:444` is
> unreachable on the original's execution path). It reports that 12 constraints point the
> same way as an existing parent→child edge of the species tree, so they are promoted to
> inviolable topological hard constraints.

## 2.4 Output files

Default **short** naming (prefix = path of the **first** constraint file):

| File | Contents |
|------|----------|
| `examples/Cyano_CUTConstraints.tsv.mt.informative.tsv` | Filtered weighted informative constraints |
| `examples/Cyano_CUTConstraints.tsv.mt.conflicts.tsv` | Constraints conflicting with the best order |
| `examples/Cyano_CUTConstraints.tsv.mt.partial_order.tsv` | Partial order (with black/green flags) |
| `examples/Cyano_CUTConstraints.tsv.html` | Interactive HTML report (unless `--no-html`) |

The fourth data file appears only when the randomization test is enabled:

| File | Appears when |
|------|--------------|
| `<prefix>.mt.random_dist.tsv` | `--random-trees > 0` (one `value similarity` pair per line, N lines) |

> - `--output-style legacy` switches back to the original long names (for reproducibility
>   verification); the file **contents are byte-for-byte identical** under both styles.
> - **Overwrites are refused by default**: if the products already exist, a second run
>   aborts with exit code `3` instead of clobbering them. Pass `-f/--force`
>   to override deliberately, or `-p` to pick another prefix.

## 2.5 Opening the HTML report

```bash
# macOS
open examples/Cyano_CUTConstraints.tsv.html
# Linux
xdg-open examples/Cyano_CUTConstraints.tsv.html
```

The report is a **fully static** page (Plotly inlined/CDN), needs no backend and opens
offline. Without `jinja2` / `plotly` installed it does **not** fail — it silently
**degrades** to a basic template (see the FAQ in chapter 09).

## 2.6 A bit more flavor

```bash
# 10 s of local search (adds the two "local search" lines + robustness summary)
maxtic-next examples/minitree.tree examples/Cyano_CUTConstraints.tsv --seed 42 --ls 10

# randomization test (adds 4 p-value lines + the 4th output file)
maxtic-next examples/minitree.tree examples/Cyano_CUTConstraints.tsv --seed 42 --rd 50

# rank directly from a RANGER-DTLx reconciliation report
# (--ale-min-family-size 0 is needed for the tiny bundled fixture)
maxtic-next examples/minitree.tree examples/adapters/ranger/FAM1.dtl \
    --from ranger --ale-min-family-size 0

# auto-detect the upstream tool output format
maxtic-next examples/minitree.tree examples/adapters/eccetera/FAM1.recphyloxml \
    --from-auto --ale-min-family-size 0

# validate the inputs without ranking (exit code 0 when they are legal)
maxtic-next examples/minitree.tree examples/Cyano_CUTConstraints.tsv --dry-run
```

Adapter mode first prints one convention statement and one diagnostic summary to
**stderr** (kept out of stdout so that the original-aligned summary stays
machine-parseable):

```
[ranger] 警告：报告未声明样本数（无 ``Total number of optimal solutions: N`` 行）且只有 1 个调和块：权重按**整数计数**输出，min_support 阈值不具 [0,1] 支持度语义。
[ranger] donor 端点约定=donor_itself 权重口径=integer_count 分母=不可确立 端点命中率=100.0% | 文件=1 块=1 转移=2 保留=2 丢弃[标签不匹配=0, 支持度=0, 家族规模=0, 家族规模探测失败=0, 自环=0]
```

(These lines are emitted in Chinese; the fields read *"donor endpoint convention =
donor_itself, weight semantics = integer_count, sample denominator = undeterminable,
endpoint hit rate = 100.0% | files=1 blocks=1 transfers=2 kept=2 dropped[label mismatch=0,
support=0, family size=0, family-size probe failure=0, self-loops=0]"*.)

The same quantities are also written to `Result.run_metadata["adapter_diagnostics"]`.
Add `--quiet-adapters` to suppress the printing only — errors are still raised.

## 2.7 One-liner in Python

```python
from maxtic_next import rank
result = rank("examples/minitree.tree", "examples/Cyano_CUTConstraints.tsv", seed=42)
print(result.best_source, result.similarity_to_input)
# greedy heuristic 1.0
```

> Prev: [01 · Installation](01_installation.md) ｜ Next: [03 · CLI reference](03_cli_reference.md)
