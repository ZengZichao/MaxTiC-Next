# 03 · CLI reference

```
maxtic-next SPECIES_TREE CONSTRAINTS [CONSTRAINTS ...] [options]
```

- `SPECIES_TREE`: path of the species tree (Newick, internal node labels carried in the
  bootstrap field). **Required**.
- `CONSTRAINTS`: one or more constraint paths. In text mode these are constraint files;
  with `--from <tool>` they are the corresponding upstream tool's output files/directory.
  **Required** (at least one).

> Both positionals **transparently accept** `.gz` / gzip streams and tar archives
> (`.tar.gz` / `.tgz` / `.tar`), so nothing has to be unpacked first; the exact support and
> its limits are in [chapter 05 §5.7](05_io_formats.md).

> The core options mirror the bare `sys.argv` semantics of the original `MaxTiC.py`;
> long and short forms are equivalent (`--local-search` = `--ls`).
> The tables below were transcribed from `build_parser()` in `src/maxtic_next/cli.py`
> (defaults from `src/maxtic_next/config.py`) and agree with `maxtic-next --help`.

## 3.0 Meta

| Option | Notes |
|--------|-------|
| `--version` | Prints `MaxTiC-Next <version>` and exits 0 (currently `MaxTiC-Next 0.1.1`) |
| `-h`, `--help` | Full option list |

## 3.1 Core ranking options (equivalent to the original)

| Option | Alias | Type / domain | Default | Notes |
|--------|-------|---------------|---------|-------|
| `--seed` | — | int ≥ 0 | `42` | Random seed driving tie-breaking in `mix` and the local search |
| `--local-search` | `--ls` | float ≥ 0 | `0.0` | Local-search **wall-clock budget** in seconds, 0 = off; turning it on also triggers the robustness/sensitivity summary |
| `--temperature` | `--t` | float > 0 | `0.001` | Metropolis temperature of the local search (scale for accepting worse solutions) |
| `--random-type` | `--r` | int ∈ {0,1,2} | `0` | Randomization type (`0` keeps data as-is; `2` fully randomizes nodes) |
| `--min-transfer-distance` | `--d` | float ≥ 0 | `0` | Minimum transfer distance, filters on the **phylogenetic distance** column, keeps `distance > threshold` |
| `--threshold-constraints` | `--ts` | float ∈ **[0, 1]** | `0.0` | Weight-threshold fraction; removes the lowest-weight edges up to that cumulative share |
| `--random-trees` | `--rd` | int ≥ 0 | `0` | Number of random trees; >0 computes empirical p-values and writes a distribution file |

> **Domains are enforced in argparse**: `--r 3`,
> `--ts 5`, `--ls -5`, `--t -1`, `--rd -1` are not accepted (they would quietly
> change the semantics); each of them exits with code **2** and an actionable message.
> `--ts` is a **fraction**, not an absolute value, because it means "delete the constraints
> making up the lowest X% of cumulative weight": `--ts 5` is illegal, while
> `--ts 0.5` deletes low-weight constraints down to a 50% share.

> ⚠️ **What `--d` actually does**: it filters on
> phylogenetic distance **only**, never on weight; and constraints **without** a distance
> column are **all kept**, not dropped. So for 3-column space input or 4-column comma
> input, `--d 5` removes nothing. This case is **error-level** (not a soft
> warning) and appends to stdout:
>
> ```
> ERROR: --min-transfer-distance=5.0 被忽略（123/123 条约束不含 phylogenetic distance 列，……）
> ```
>
> (message emitted in Chinese: *"ignored — 123/123 constraints carry no phylogenetic
> distance column"*). The run still exits **0**, because this is a statement of convention
> and does not alter the computation; `--dry-run` reports it as an error and exits 1.
> To filter by distance, supply 4-column (space) or 5-column (ALE comma) input; otherwise
> drop the option.

## 3.2 Upstream-tool options

| Option | Values | Default | Notes |
|--------|--------|---------|-------|
| `--from` | `ale`/`ranger`/`eccetera`/`artra`/`alerax`/`auto` | none | Treat CONSTRAINTS as that tool's output and convert via the registry |
| `--from-auto` | — | off | Equivalent to `--from auto`: auto-detect format from content/extension |
| `--from-ale` | — | off | Backwards-compatible switch, equivalent to `--from ale` (CONSTRAINTS are `.uml_rec` files) |
| `--ale-min-support` | float ∈ [0,1] | `0.05` | Minimum within-family support/frequency, keeps `> threshold` (applies to all reconciliation adapters) |
| `--ale-min-family-size` | int ≥ 0 | `5` | Minimum gene family size (reconciled-tree leaf count); families at or below the threshold are skipped (not applicable to alerax) |
| `--ale-cache-dir` | path | none | Per-family pickle cache directory (resumable) |
| `--ale-source` | `rec`/`trf` | **`trf`** | ALE constraint source — see the dedicated note below |
| `--ale-parallel` | `process`/`thread` | `process` | Parallel parsing: process pool (default, sidesteps the GIL) / thread pool |
| `--artra-transfer-kind` | `all`/`replacing`/`additive` | `all` | Which ARTra transfer classes to count |
| `--min-endpoint-hit-rate` | float ∈ [0,1] | `0.5` | Adapter **endpoint hit-rate** threshold: below it a stderr warning is emitted; a 0% hit rate is a hard error. Set to `0` to disable the warning (0 hits still errors) |
| `--quiet-adapters` | — | off | Suppress the adapters' routine stderr convention statement and diagnostic summary; **printing only** — diagnostics remain in `Result.run_metadata["adapter_diagnostics"]` and errors are still raised |

> The `--ale-` prefix is historical: `--ale-min-support`, `--ale-min-family-size`,
> `--ale-cache-dir`, `--min-endpoint-hit-rate` and `--quiet-adapters` also apply to
> `ranger/eccetera/artra/alerax`.

### 3.2.1 `--ale-source`: the default is `trf`

| Value | Generated from | Semantics |
|-------|----------------|-----------|
| **`trf` (default)** | **transfer** events; equivalent to `constraints_from_transfers` in ALE's own MaxTiC integration | "the **parent** of the donor branch is older than the **child** of the receptor branch" |
| `rec` | **reconciliation** events (the reconciliation records in `.uml_rec`) | produces a constraint set of **different** size and weight scale |

- `trf` is the default because **ALE's own** MaxTiC integration documents
  `constraints_from_transfers` (trf) as the direct input to MaxTiC. These are not two
  interchangeable spellings: **switching replaces the entire input constraint set** (edge
  count and total weight differ by a factor of several on the official example), so every
  downstream number changes.
- The effective convention is printed to stderr, e.g.
  `[ale] 约束口径：source=trf（……ALE 官方 MaxTiC 集成口径），供体端点约定=parent_of_donor，产出边数=…，总权重=…`
- To reproduce results computed from **reconciliation** events, pass `--ale-source rec` explicitly.

### 3.2.2 The adapter diagnostic line (stderr)

Each conversion prints one line with a fixed set of fields:

```
[tool] donor 端点约定=<parent_of_donor|donor_itself> 权重口径=<support_over_declared_samples|support_over_replicate_blocks|integer_count|posterior_frequency> 分母=<N|不可确立> 端点命中率=NN.N% | 文件=N 块=N 树=N 转移=N 保留=N 丢弃[标签不匹配=N, 支持度=N, 家族规模=N, 家族规模探测失败=N, 自环=N]
```

(The line is emitted in Chinese; the fields read *"donor endpoint convention / weight
semantics / sample denominator / endpoint hit rate | files / blocks / trees / transfers /
kept / dropped [label mismatch, support, family size, family-size probe failure,
self-loops]"*.)

- `donor endpoint convention`: ALE uses `parent_of_donor` (the donor's parent), while
  RANGER-DTLx / ARTra / ecceTERA / AleRax use `donor_itself` (the donor itself) — they
  differ by **one species level**. Mixed `--from-auto` input warns separately per tool instead of
  silently summing the two.
- `weight semantics`: ALE yields [0,1] support; RANGER/ARTra yield **integer counts** when
  the report declares no sample count, in which case `--ale-min-support` has no [0,1]
  meaning. Weights of different kinds are **not comparable**.
- `自环=N` (*self-loops*): the number of events with `donor == receptor` (donor and receptor
  are the same lineage). This is a **counter**, not "N constraints were discarded": self-loop
  constraints are **kept** and contribute their weight to the original
  **`to itself`** statistic in the stdout summary (see 3.2.3). All five adapters,
  ALE included, behave this way.
- The same dictionary is available from Python: `result.run_metadata["adapter_diagnostics"]`.

### 3.2.3 Self-loop transfers and `--d`

When donor and receptor are the same lineage, their species-tree topological distance is
always 0. Had the adapter written that 0 into the distance column,
`--min-transfer-distance` would have removed the event as "too close", and it would have
**silently vanished** from the `to itself` statistic. **All five adapters — RANGER-DTLx,
ecceTERA, ARTra, AleRax and ALE** — write **`distance = None`** (no distance information)
for self-loops, so under the existing rule "constraints without a distance column are always
kept" they are **not eaten by `--d`**; the regression test is
`tests/test_ale_selfloop.py`.

Measured (input A: `examples/adapters/ale/species.tree` plus a RANGER-DTLx report in which one
transfer line has its `Recipient -->` set to its own `Mapping` label (62 → 62)):

```
$ maxtic-next species.tree selfloop.dtl --from ranger --d 1
[ranger] donor 端点约定=donor_itself … 端点命中率=100.0% | 文件=1 块=1 转移=2 保留=2
         丢弃[标签不匹配=0, 支持度=0, 家族规模=0, 家族规模探测失败=0, 自环=1]
1.0 uninformative (  0.0 to a descendant, 0.0 to a leaf 0.0 to an ancestor 1.0 to itself 50.0%)
```

The `to itself` figure is identical for `--d 1` and `--d 0` (both `1.0`), i.e. the self-loop
survived the distance filter.

Measured (input B: **the same species tree through `--from ale`** — a `.uml_rec` containing
`T@62->65`, which under the default `trf` convention means donor = `parent(62)` = `65` = the
receptor, i.e. a genuine "transfer to its own lineage". The file is exactly
`1 reconciled` / `FAM_SELFLOOP` /
`(((CYAP8:1.0,CYAP0:1.0)59:1.0,(CYAA5:1.0,CYAP2:1.0)61.T@61->59:1.0)65:1.0,(NOSP7:1.0,ANAVT:1.0)62.T@62->65:1.0)69:1.0;`):

```
$ maxtic-next species.tree selfloop.uml_rec --from ale
[ale] 约束口径：source=trf（trf=转移事件 parent(donor)->receptor，ALE 官方 MaxTiC 集成口径），供体端点约定=parent_of_donor，产出边数=2，总权重=2.0
[ale] donor 端点约定=parent_of_donor 权重口径=support_over_declared_samples 分母=1 端点命中率=100.0% | 文件=1 块=1 转移=2 保留=2 丢弃[标签不匹配=0, 支持度=0, 家族规模=0, 家族规模探测失败=0, 自环=1]
selfloop.uml_rec
tree with  5 internal nodes
2.0 total weight of constraints from transfers
1.0 uninformative (  0.0 to a descendant, 0.0 to a leaf 0.0 to an ancestor 1.0 to itself 50.0%)
```

(`约束口径` = "constraint convention", `供体端点约定` = "donor endpoint convention",
`权重口径` = "weight semantics", `端点命中率` = "endpoint hit rate", `转移` = "transfers",
`保留` = "kept", `丢弃` = "dropped", `自环` = "self-loops".) ALE agrees with the other four:
`自环=1` is only a counter, and that constraint still contributes the `1.0` seen under
`to itself`; nothing was dropped at the default `--d 0`, so no NOTE was printed. The two-stage
product says the same thing — the 5th column of the self-loop row is left **empty** instead of
carrying `0.0`:

```
$ maxtic-next species.tree selfloop.uml_rec --from ale -o unified.tsv
[ale] … (the two stderr diagnostic lines above, elided) …
[约束生成] ale 约束已写出至：unified.tsv        (*"ale constraints written to unified.tsv"*)
$ cat unified.tsv
#family,donor,receptor,weight,distance
selfloop,65,59,1.0,2.0
selfloop,65,65,1.0,
```

**The text-input path behaves the same way**: 3-column space or 4-column comma input
carries no distance column, so self-loops are always kept and always counted under
`to itself`. Only when **you** write an explicit 0 distance in the text input does the event
get dropped, and then the program says so:

```
NOTE: --min-transfer-distance=1.0 丢弃 1/2 条约束（保留 1 条）。
NOTE: 上述丢弃中有 1 条是**自环约束**（donor == receptor，其拓扑距离恒为 0）；
      它们不会进入原版口径的 "to itself" 统计。
```

(Chinese: *"dropped 1 of 2 constraints (kept 1)"* and *"1 of those was a self-loop
constraint (donor == receptor, whose topological distance is always 0); it will not reach the
original 'to itself' statistic"*.) In other words, whether `to itself` is non-zero depends on
"does the input carry a distance column with value 0", not on the event being a self-loop.

## 3.3 Two-stage workflow option

| Option | Alias | Notes |
|--------|-------|-------|
| `-o` | `--constraints-out` | **`--from <tool>` / `--from-ale` modes only**: write the converted unified constraint set to this path (5-column ALE comma format, re-readable by the text parser), then **return early without ranking**, for Snakemake/Nextflow reuse. Combined with `--dry-run` the file **is still written**, and `-o` is never silently ignored. Passing `-o` in text mode is a **hard error** (the input is already unified, so ignoring it silently would mislead users) |

## 3.4 Pruning and preflight

| Option | Notes |
|--------|-------|
| `--target-clade <label>` | Root label of the target clade; by default only constraints with both endpoints inside the clade are kept  |
| `--target-clade-ancestor-map` | Enable the conservative "external ancestor mapped to the root" rule (requires `--target-clade`; off by default to avoid spurious partial orders) |
| `--dry-run` | Static preflight: tree format / binarity / internal-label uniqueness / leaf-name uniqueness / constraint format / weight validity / endpoint validity / `d`-column semantics. In `--from` mode it **actually parses the upstream output** and reports zero-constraint and low-hit-rate problems. **Any adapter parse failure makes the preflight FAIL** and aborts with **exit code 1** (worked example in 08 §8.4b). Does not rank |

## 3.5 Reporting and output

| Option | Values | Default | Notes |
|--------|--------|---------|-------|
| `--no-html` | — | off (HTML **is** produced by default) | Disable the interactive HTML report |
| `--output-style` | `short`/`legacy` | `short` | Output naming style; `legacy` gives the original long names. **Contents are byte-for-byte identical** under both styles |
| `-p` | `--output-prefix` | none | Output prefix; defaults to the path of the **first** constraint file (see below) |
| `-f` | `--force` | off | Allow overwriting existing output files. By default overwrites are **refused** with exit code **3** |

> **Prefix and multiple files**: the output prefix comes from the **first**
> file in the `CONSTRAINTS` list; the remaining files take part in the computation but do
> not appear in the file names. With several files and no `-p`, one declaration is written
> to **stderr**:
>
> ```
> NOTE: 输入了 2 个约束文件，三输出文件与 HTML 报告均以第一个文件为前缀：<file1>
> （其余 1 个文件已参与计算：['<file2>']）。如需为不同来源分别留存产物，
> 请对每次运行显式指定 -p/--output-prefix。
> ```
>
> (Chinese: *"2 constraint files given; the three output files and the HTML report all use
> the first file as prefix … the other 1 file did take part in the computation … pass
> -p/--output-prefix explicitly to keep products per source."*)
>
> It is also stored in `result.run_metadata["multi_file_prefix_note"]`.
>
> **Write semantics**: products are written to a temporary file in the target directory,
> `fsync`-ed and then moved with `os.replace`, so a crash mid-write cannot leave a
> half-written file; missing parent directories are created automatically (failures report
> the offending path).

## 3.6 Robustness / MCMC (advanced)

| Option | Type / domain | Default | Notes |
|--------|---------------|---------|-------|
| `--mcmc` | — | off | Reversible Metropolis–Hastings sampling over linear extensions compatible with the species-tree topology. **Note: this chain has not passed convergence diagnostics; its samples must not be treated as posterior samples** |
| `--mcmc-iters` | int > 0 | `1000` | Number of MH chain steps |
| `--mcmc-temperature` | `auto` or float > 0 | **`auto`** | Boltzmann temperature. `auto` = `max(total weight, 1) / 100`; the effective value is printed and recorded in `run_metadata["params"]`. A positive number fixes the temperature |
| `--mcmc-burn-in` | int ≥ 0 | **50%** of `iters` | Burn-in steps: advance the chain without collecting |
| `--mcmc-thin` | int > 0 | `max(1, iters // 100)` | Thinning interval, to reduce autocorrelation between retained samples |
| `--near-optimal-top-k` | int ≥ 1 | `50` | **Near-optimal collector capacity**: how many **deduplicated** near-optimal orders the robustness/sensitivity summary may consider (API `top_k`, see 3.6b). Meaningful only with `--ls > 0` |

`--mcmc` appends a convention statement, the effective parameters, sample statistics and
**convergence diagnostics** to stdout (measured: 1000 steps, `--seed 42`, cyanobacterial
example):

```
[MCMC] Metropolis–Hastings over linear extensions (preliminary; convergence diagnostics not validated)
[MCMC] 口径声明：这是一条**初步实现**的 Metropolis–Hastings 链，收敛诊断未经验证；链上样本彼此自相关（非独立观测），未通过下列诊断前不得用作不确定性度量
[MCMC] 温度 T=22.184（auto：max(总权重,1)/100）；iters=1000、burn-in=500、thin=10
[MCMC] 保留 100 个样本，其中唯一排序 51 个；能量 mean=835.7672, min=795.6300（唯一样本数远低于样本数说明链未混合，均值不可用作不确定性度量）
[MCMC] 收敛诊断（convergence diagnostics，按抽稀前的记录链计）：唯一样本占比=0.074、相邻重复率=0.907、最长同态连续段占比=0.056、有效样本量 ESS=7.6（0.008 倍样本量）、状态覆盖=n/a（状态空间过大，未精确枚举）；能量 sd=21.2451
```

> The default temperature is **not** a fixed `0.01`: on the example
> with total weight 2218.4, `auto` yields T=22.184 — three orders of magnitude apart.
> At 0.01 the chain is effectively frozen and explores very few states.
> `--mcmc-temperature 0.01` fixes it explicitly when you want that value.

### 3.6b `--near-optimal-top-k`: the support-set cap of the robustness summary

`--near-optimal-top-k` (API `api.rank(top_k=...)`, `Ranker.run(top_k=...)`) **caps** how many
**distinct** (deduplicated) near-optimal orderings the "robustness/sensitivity summary over
solutions visited by the local search" may consider:

- default **50** (`config.DEFAULT_NEAR_OPTIMAL_TOP_K`): the collector holds at most K
  deduplicated orders, so with the default the 51st order can never take part in the
  counts;
- values are validated by `ranking/ranker.py:check_near_optimal_top_k()`: it must be an
  **integer ≥ 1**. Out-of-range on the CLI (e.g. `--near-optimal-top-k 0`) is rejected by
  argparse with **exit code 2**
  (`argument --near-optimal-top-k: 近优解收集容量 必须 > 0，收到 0`); on the Python side
  `0` / `-1` / `2.5` / `'50'` / `True` / `inf` all raise `ValueError`, while the integral
  float `100.0` normalizes to `100`;
- the **effective** value is recorded in `Result.run_metadata["near_optimal_top_k"]`, echoed
  in the `run_metadata["params"]` string (`--near-optimal-top-k 12`) and in the HTML report
  (`top_k：12（本次报告使用 12）`); the summary dict also carries `top_k` / `top_k_reported`.

**Why the cap matters**: with `--ls > 0` what the program reports is **not** robustness under
resampling of the input data (`is_resampling_robustness` is always `False`) — it is a
**deduplicated view of the set of solutions the local search visited**. K therefore bounds the
**support** of every frequency in it: orders outside K take no part in the counts. Measured
(`examples/minitree.tree` + `examples/Cyano_CUTConstraints.tsv`, `--ls 1`,
`--local-search-max-iters 60000`, `--t 60`, `--seed 42`):

| `--near-optimal-top-k` | deduplicated orders retained and counted | distinct orders the search actually visited |
|------------------------|-------------------------------------------|---------------------------------------------|
| `3` | 3 | 9141 |
| `10` | 10 | 9139 |
| `50` (default) | 50 | 8982 |
| `400` | 400 | 7940 |

The chain visits far more distinct orders than K; the summary only sees the K best by
objective. Raising K widens the near-optimal neighbourhood the summary can speak about, at
the cost of memory and `summary()` work growing like **K·n²** (pairwise order frequencies are
accumulated over K orders of length n).

> ⚠️ **The terminology rule still applies**: raising K does **not** turn this into a
> posterior. The module remains a "robustness/sensitivity summary" and
> "posterior" / "confidence" / "confidence interval" stay **forbidden** (see 08 §8.1).
> Whenever you report those frequencies, state the effective `--near-optimal-top-k` too —
> without K, a number like "node 61 sits at rank 3 with frequency …" has no interpretable
> denominator.

## 3.7 Performance and resuming (advanced)

| Option | Type / domain | Default | Notes |
|--------|---------------|---------|-------|
| `--incremental` / `--no-incremental` | — | **on** | Incremental scoring: each interval rotation in the local search costs O(b−a) instead of a full O(\|E\|) recompute. On by default (byte-identical to the full path, measured 3.1–8.2× faster); `--no-incremental` falls back to the original full recompute |
| `--checkpoint <path>` | str | none | Checkpoint file; a search interrupted midway can be resumed (only effective with `--local-search>0`) |
| `--checkpoint-interval` | float > 0 | `60.0` | Automatic checkpoint save interval in seconds |
| `--local-search-max-iters` | int ≥ 0 | `0` | **Iteration cap for the local search** (0 = bounded by `--ls` only). Setting a positive number makes the fixed-`--seed` local-search path **reproducible across machines** |

## 3.8 The exact boundary of reproducibility

| Path | Reproducible under a fixed `--seed`? | Notes |
|------|--------------------------------------|-------|
| Default (`--ls 0`) | **Yes**, byte-for-byte | All three heuristics and the outputs are driven by the seed-determined random stream; constraint line order and multi-file order do not change the products |
| `--ls > 0`, no iteration cap | **Not guaranteed across machines** | The stop condition is wall-clock time, so the iteration count follows machine load. The program says so itself: a WARNING is appended to stdout |
| `--ls > 0` + `--local-search-max-iters N` | **Yes** (same Python/platform) | The iteration count and hence the random-stream consumption are deterministic; this is the option to use when you need a deterministic search |
| `--r 1` / `--r 2` | Determined by `--seed` | Randomization alters the data itself, so "agrees with the default path" is meaningless, but the same seed still reproduces |

> In other words, the blanket promise that "a fixed `--seed` is reproducible" holds only in
> this precise form: **the default path (`--ls 0`) is reproducible;
> `--ls > 0` needs `--local-search-max-iters` to be reproducible across machines.**

## 3.9 Command cheatsheet

```bash
# basic
maxtic-next species.tree constraints.tsv --seed 42

# 180 s of local search + no HTML (add an iteration cap for cross-machine determinism)
maxtic-next species.tree constraints.tsv --ls 180 --local-search-max-iters 2000000 --no-html

# filter by distance (needs a distance column, otherwise the option is ignored + ERROR)
maxtic-next species.tree ale_constraints.tsv --d 3

# weight-threshold filtering + random-tree p-values
maxtic-next species.tree constraints.tsv --ts 0.1 --rd 1000

# upstream tools (explicit / auto)
maxtic-next species.tree FAM1.dtl --from ranger
maxtic-next species.tree out.recphyloxml --from-auto

# ALE: select the `rec` (reconciliation-event) convention explicitly (the default is `trf`)
maxtic-next species.tree fam.uml_rec --from ale --ale-source rec

# target-clade pruning
maxtic-next species.tree constraints.tsv --target-clade 61

# preflight
maxtic-next species.tree constraints.tsv --dry-run

# MCMC (auto temperature + default burn-in/thin; --mcmc-temperature 0.01 to fix it)
maxtic-next species.tree constraints.tsv --mcmc --mcmc-iters 2000

# let the robustness/sensitivity summary consider more near-optimal topologies (default K=50)
maxtic-next species.tree constraints.tsv --ls 180 --near-optimal-top-k 400

# compressed input, fed directly (no unpacking first; see 05 §5.7)
maxtic-next species.tree.gz constraints.tsv.gz --seed 42
maxtic-next species.tree examples/reconciliations.tgz --from ale

# suspicious upstream output: the preflight really parses it, a parse failure exits 1
maxtic-next species.tree out.recphyloxml --from eccetera --dry-run; echo "exit=$?"

# long search + incremental + resumable
maxtic-next species.tree constraints.tsv --ls 600 --incremental \
    --checkpoint run.ckpt --checkpoint-interval 30

# re-run into the same directory, deliberately allowing overwrite of previous products
maxtic-next species.tree constraints.tsv --seed 42 --force
```

## 3.10 Exit codes

| Code | Meaning |
|------|---------|
| `0` | Success (including a passed `--dry-run`) |
| `1` | `--dry-run` found **error-level** problems (decided from the structured `run_metadata["dry_run_ok"]`, never by substring matching). In `--from` mode "the upstream output cannot be parsed at all" lands here too |
| `2` | argparse domain error (illegal `--r`, out-of-range `--ts`, negative `--ls`, unknown `--ale-source`, `--near-optimal-top-k < 1`, …), returned by argparse itself |
| `3` | Output file already exists and `--force` was not given (overwrite protection) |
| `4` | Runtime exception: an input fails an algorithmic precondition (non-binary tree, duplicate leaf names, duplicate internal labels, constraint endpoints absent from the tree, illegal/negative weights, 0% adapter hit rate, …) |

> In a container this code propagates directly (`ENTRYPOINT ["maxtic-next"]`).

> **Front-end exit codes**: the desktop app `maxtic-studio` lives in the separate repository
> `MaxTiC-Next-Studio` (see "Where is the graphical front end?" in the README). It reuses this
> same table, where **`3`** additionally means "GUI dependencies missing — an install hint is
> printed". This repository does not register that entry point.

> Prev: [02 · Quickstart](02_quickstart.md) ｜ Next: [04 · Python API](04_python_api.md)
