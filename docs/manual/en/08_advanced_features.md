# 08 · Advanced features

## 8.1 Robustness / sensitivity summary (local search)

With local search enabled (`--local-search T > 0`), MaxTiC-Next collects a **set of
near-optimal solutions** during the search and summarizes:

- **node position distribution** (`node_position_distribution`): for each internal node, how
  often it appears at each rank among the near-optimal orders;
- **pairwise order frequency** (`pairwise_order_frequency`): for each ordered pair (a, b),
  how often a precedes b.

```bash
maxtic-next species.tree constraints.tsv --ls 30
```

```python
r = rank("species.tree", "constraints.tsv", local_search=30, print_summary=False)
summ = r.sensitivity_summary
print(summ["title"])                 # robustness/sensitivity summary over solutions visited by the local search
print(summ["n_orders_collected"])    # number of deduplicated near-optimal orders collected
print(summ["node_position_distribution"]["61"])   # position distribution of node 61
```

> **Terminology rule**: the local search is a Metropolis neighbourhood search and the
> solutions it visits are **not** a statistical posterior. This module therefore speaks only
> of a "robustness/sensitivity summary"; "confidence interval / posterior probability /
> posterior distribution" are **forbidden**. Raising the capacity below does not loosen this
> rule one bit.

### 8.1b `--near-optimal-top-k`: the hard cap on the summary's support

What the summary tallies is not resampling of the input but the **set of solutions the local
search visited**, deduplicated (`is_resampling_robustness` is always `False`). The collector
keeps only **K** deduplicated near-optimal orders, ordered by objective value, so K is the
**support cap** of every frequency reported below: orders outside K take no part in the
counts. **K defaults to 50** and is overridable from three layers —
CLI `--near-optimal-top-k`, `api.rank(top_k=...)`, `Ranker.run(top_k=...)` — validated by
`ranking/ranker.py:check_near_optimal_top_k` (**integer ≥ 1**; exit code 2 on the CLI,
`ValueError` in Python), recorded in `Result.run_metadata["near_optimal_top_k"]`, echoed in
`run_metadata["params"]` and in the HTML report, with `top_k` / `top_k_reported` inside the
summary dict.

Measured (`examples/minitree.tree` + `examples/Cyano_CUTConstraints.tsv`, `--ls 1`,
`--local-search-max-iters 60000`, `--t 60`, `--seed 42`):

| `--near-optimal-top-k` | deduplicated orders retained and counted | distinct orders the search actually visited |
|------------------------|-------------------------------------------|---------------------------------------------|
| `3` | 3 | 9141 |
| `10` | 10 | 9139 |
| `50` (default) | 50 | 8982 |
| `400` | 400 | 7940 |

```bash
maxtic-next species.tree constraints.tsv --ls 180 --near-optimal-top-k 400
```

```python
r = rank("species.tree", "constraints.tsv", local_search=180, top_k=400)
r.run_metadata["near_optimal_top_k"]        # 400 — the machine-readable entry point
r.sensitivity_summary["n_orders_retained"]  # deduplicated orders that took part
```

- A **larger K** lets the summary speak about a wider near-optimal neighbourhood; the price is
  memory and `summary()` work growing like **K·n²** (pairwise order frequencies are accumulated
  over K orders of length n).
- Raising K does **not** make it a posterior and does **not** improve chain mixing; it only
  changes where the denominator comes from.
- Always publish the effective K together with these frequencies, otherwise a statement like
  "node 61 lands at rank 3 with frequency …" has no interpretable denominator.

> **Reproducibility**: with `--ls > 0` the stop condition is **wall-clock time**, so the
> iteration count follows machine load and a fixed `--seed` is **not guaranteed to reproduce
> across machines** (the program emits a WARNING saying exactly this). For determinism also
> set `--local-search-max-iters N` (see 8.5b). The default path (`--ls 0`) is always
> byte-for-byte reproducible.

## 8.2 MCMC sampler (optional, preliminary)

`--mcmc` enables a reversible Metropolis–Hastings sampler over the state space of **linear
extensions compatible with the species-tree topology** (symmetric proposal +
`min(1, exp((E_old−E_new)/T))`, so detailed balance holds).

```bash
maxtic-next species.tree constraints.tsv --mcmc --mcmc-iters 2000
```

| Option | Default | Notes |
|--------|---------|-------|
| `--mcmc-iters` | `1000` | Chain steps (must be > 0) |
| `--mcmc-temperature` | **`auto`** | Boltzmann temperature. `auto` = `max(total weight, 1) / 100`; a positive number fixes it |
| `--mcmc-burn-in` | **50%** of `iters` | Burn-in: advance the chain without recording samples |
| `--mcmc-thin` | `max(1, iters // 100)` | Thinning interval, reducing autocorrelation between retained samples |

- **The default temperature is `auto`, not a fixed `0.01`.** At real total
  weights (hundreds to thousands) a fixed `0.01` nearly freezes the chain, so it explores very
  few states. `auto` scales with the energy, and the effective
  value is printed to stdout:

  ```
  [MCMC] 温度 T=22.184（auto：max(总权重,1)/100）；iters=1000、burn-in=500、thin=10
  ```

- Results land in `Result.mcmc_samples` (list of linear-extension samples), and sample
  statistics plus **convergence diagnostics** are printed:

  ```
  [MCMC] 保留 100 个样本，其中唯一排序 51 个；能量 mean=835.7672, min=795.6300（唯一样本数远低于样本数说明链未混合，均值不可用作不确定性度量）
  [MCMC] 收敛诊断（convergence diagnostics，按抽稀前的记录链计）：唯一样本占比=0.074、相邻重复率=0.907、最长同态连续段占比=0.056、有效样本量 ESS=7.6（0.008 倍样本量）、状态覆盖=n/a（状态空间过大，未精确枚举）；能量 sd=21.2451
  ```

  (Chinese fields read: *"kept 100 samples, of which 51 unique orders; energy mean/min (the
  unique-sample count being far below the sample count shows the chain has not mixed, so the
  mean must not be used as an uncertainty measure)"* and *"convergence diagnostics (counted
  over the recorded chain before thinning): unique-sample fraction, adjacent-duplicate rate,
  longest homogeneous-run fraction, effective sample size ESS, state coverage, energy sd"*.)

> ⚠️ **Honest capability boundary**: this is a **preliminary implementation that
> has not passed convergence diagnostics**. The first line it prints says
> `preliminary; convergence diagnostics not validated`, and the measurements show the chain
> is **unmixed** (ESS far below the sample count, high adjacent-duplicate rate). Therefore
> **until the diagnostics pass, its samples must not be used as an uncertainty measure**.
> MCMC is merely the one module where the vocabulary of "distribution" is *permitted* — that
> does not mean it has earned it. **Off by default**; when off, the P0/P1 paths are unchanged
> (zero drift).
> The sampler also **honours cooperative cancellation** (`MCMCSampler.sample(stop_check=...)`,
> threaded through by `api.rank` / `Ranker.run`): the predicate is checked between chain steps,
> the samples recorded so far are returned, and `run_metadata["mcmc_stopped_early"]` reports
> whether the chain really was interrupted. A chain cancelled **before producing any sample**
> raises `RuntimeError` instead of silently returning empty statistics (see 8.9).

## 8.3 Target-clade pruning

When you only care about a subclade of the species tree, `--target-clade <label>` focuses the
ranking on that clade:

```bash
# default: keep only constraints whose both endpoints are inside clade 61
maxtic-next species.tree constraints.tsv --target-clade 61

# enable the conservative "external ancestor mapped to the root" rule (use with care)
maxtic-next species.tree constraints.tsv --target-clade 61 --target-clade-ancestor-map
```

- Default policy: cross-clade constraints are dropped, and **no spurious partial order is
  ever introduced silently**.
- `--target-clade-ancestor-map`: an external endpoint is mapped to the root only if it can be
  determined to be an ancestor of the clade and lies on the donor side; side branches are
  always dropped. Off by default.
- This step **does not modify tree topology or re-root**; it only filters / remaps constraint
  endpoints.
- A `--target-clade` label absent from the tree is an **error** in both `--dry-run` and a
  real run, not a deep `KeyError`.
- Caveat: a constraint mapped to the root holds in every legal order (the root's rank is 0),
  so it only dilutes the weight denominator without carrying information.

## 8.4 Preflight (dry-run)

Statically checks the inputs before ranking — **it consumes no random numbers and modifies
no data**, so it is safe to put at the front of a pipeline:

```bash
maxtic-next species.tree constraints.tsv --dry-run
maxtic-next species.tree fam1.uml_rec --from ale --dry-run   # actually parses upstream output
```

Checks performed (transcribed from `src/maxtic_next/dry_run.py`):

| Check | Severity |
|-------|----------|
| Species tree parses and is rooted | error |
| Internal node labels unique | error |
| **Binary-tree assumption** (each internal node has exactly 2 children) | error |
| **Leaf names unique** | error |
| `--target-clade` label exists in the species tree | error |
| Constraint file exists, readable, both formats parse (3-column comma misuse rejected) | error |
| Constraint weight/distance are finite reals and weight ≥ 0 (with file name and line number) | error |
| All constraint endpoints present in the species tree | error |
| `--d > 0` but constraints carry **no** distance column (option would be ignored) | **error** |
| Species tree has no internal nodes | warning |
| Adapter mode: upstream output fails to parse / zero constraints / no endpoint hits at all | error |
| Adapter mode: endpoint hit rate below `--min-endpoint-hit-rate` | warning |

The output is a structured report ending with `结论：通过 / 存在致命问题，已中止
（error N, warning N）` ("conclusion: passed / fatal problems found, aborted"). On legal
input (`examples/minitree.tree` + `examples/Cyano_CUTConstraints.tsv`) the measured output is:

```
=== 预检（dry-run）报告 ===
未发现任何问题。
结论：通过 （error 0, warning 0）
```

> **The exit code comes from the structured field** `run_metadata["dry_run_ok"]`,
> never from substring matching on the report text. Passing preflight exits `0`; any
> error-level finding exits `1`.
> Preflight in adapter mode **actually parses the upstream output**: two
> lines of garbage posing as a `.dtl` file must not be judged "passed", because dropping
> `--dry-run` would then produce a complete, plausible-looking summary on 0 total weight.
> Zero-constraint and no-endpoint-hit cases are error-level. When `--dry-run` and `-o` are both
> given, `-o` still writes the constraint file.

### 8.4b A parse failure always fails the preflight

In adapter mode the preflight really runs the upstream conversion once. **Any adapter parse
failure** — corrupt/truncated XML, an unreadable file, a multi-member archive — is converted
into an **error-level** report entry, after which `n_errors` / `ok` are **recomputed**, so the
CLI aborts with **exit code 1**. Appending the error without recomputing `ok` would print
`[error]` while concluding *passed* and exiting 0, which is exactly what this recomputation
prevents.

Worked example: `examples/adapters/eccetera/FAM1.recphyloxml` truncated to 400 bytes
(`head -c 400`), leaving the XML unclosed:

```bash
$ head -c 400 examples/adapters/eccetera/FAM1.recphyloxml > truncated.recphyloxml
$ maxtic-next examples/adapters/ale/species.tree truncated.recphyloxml \
      --from eccetera --dry-run
=== 预检（dry-run）报告 ===
[error] 上游输出（格式 eccetera）解析失败：UpstreamParseError: /tmp/truncated.recphyloxml：
        recPhyloXML 解析失败（文件损坏 / 被截断 / 非 XML）（文件 /tmp/truncated.recphyloxml，
        第 13 行，第 11 列）；原始异常：ParseError: unclosed token: line 13, column 11
结论：存在致命问题，已中止 （error 1, warning 0）
$ echo $?
1
```

(The message is emitted in Chinese: *"upstream output (format eccetera) failed to parse:
recPhyloXML parse failed (file corrupt / truncated / not XML) … original exception: unclosed
token: line 13, column 11"*, and the closing line reads *"conclusion: fatal problems found,
aborted"*. Captured with
`PYTHONPATH=src python3 -m maxtic_next … --from eccetera --dry-run`; wrapped here for
legibility.) For contrast, the **untruncated** file passes preflight with exit `0`
(`未发现任何问题。/ 结论：通过 （error 0, warning 0）`). From Python, decide with
`result.run_metadata["dry_run_ok"]` / `dry_run_n_errors`, never by matching text.

## 8.5 Incremental scoring (performance)

The local search's neighbourhood moves are **interval rotations**. `--incremental` computes only
the increment over the affected interval (O(b−a)) instead of **recomputing `value()` in full**
(O(|E|)) after each move, which is much faster for large |E| / long searches, with identical
results:

```bash
maxtic-next species.tree constraints.tsv --ls 600 --incremental
```

> **On by default**: the incremental and full paths are verified byte-for-byte identical
> over 4 inputs × 5 seeds × 7 artifact fields (guarded by `tests/test_output_equivalence.py`),
> and the incremental path is 3.1–8.2× faster. It affects **performance only** and changes
> neither the search trajectory nor the result; use `--no-incremental` to fall back to the
> original full recompute.

## 8.5b Making the local search reproducible across machines

```bash
maxtic-next species.tree constraints.tsv --ls 600 --local-search-max-iters 2000000
```

- `--local-search-max-iters 0` (default) = bounded only by the `--ls` wall-clock budget →
  the **iteration count varies with machine load**, so the same seed can yield different
  results on different machines (the program appends a WARNING to stdout saying so).
- A positive integer N → stop after at most N iterations. Random-stream consumption becomes
  deterministic, so a **fixed `--seed` is deterministic**; `--ls` then degrades to an
  additional upper bound on runtime.
- The default path (`--ls 0`) is unaffected and always reproducible.

## 8.6 Checkpointing / resuming

Long local searches can be interrupted. `--checkpoint` periodically saves the search state,
and re-running the same command resumes from it:

```bash
maxtic-next species.tree constraints.tsv --ls 3600 \
    --checkpoint run.ckpt --checkpoint-interval 30
# if interrupted, repeat the command above to resume from run.ckpt
```

- Only effective with `--local-search > 0`.
- `--checkpoint-interval`: automatic save interval in seconds (default 60).
- Checkpointing makes a **long search resumable**, but it does **not** by itself provide
  cross-machine determinism — that is what `--local-search-max-iters` is for (8.5b). The two
  are complementary: checkpoints solve "it will not finish", caps solve "different machines
  give different answers".

> Note: the ALE **parsing stage** has its own per-gene-family file-level cache
> (`--ale-cache-dir`), an independent resume mechanism. Its cache key includes a
> **species-tree topology fingerprint**, so two trees with the same label set
> but different topology never contaminate each other's results.

## 8.7 Random-tree significance (p-values)

`--random-trees N` generates N random orders that are **legal topological orders**, computes
empirical p-values for the value and similarity of the **final delivered order**, and writes a
distribution file:

```bash
maxtic-next species.tree constraints.tsv --rd 1000
```

Measured output with `--rd 50` on the cyanobacterial example:

```
values from  50  random orders 781.3100000000004 1058.9899999999996
pvalue of the found order: 0.0196078431372549 [(k+1)/(n+1) corrected; tested order = delivered best order]
similarity values from  50 random orders 0.4 0.9111111111111111
pvalue of the similarity with the input order: 0.0196078431372549 [(k+1)/(n+1) corrected]
```

Two statistical conventions define these p-values:

1. **The `(k+1)/(n+1)` correction** (Phipson & Smyth 2010). Reporting `k/n`, as the original
   does, would at N=50 squeeze the smallest non-zero p-value down to `0.0` and manufacture
   an illusion of absolute significance; with 50 samples the minimum p reported here is
   `1/51 ≈ 0.0196`. The printout states `[(k+1)/(n+1) corrected]` explicitly.
2. **The tested order is the final delivered order** (after local search), not the
   pre-search heuristic order — testing an order that is not the tree handed to
   the user would describe a different tree. The printout states
   `tested order = delivered best order`.

The distribution file is the **4th output product**: `<prefix>.mt.random_dist.tsv`
(`_distribution_random` under the legacy style), one `value similarity` pair per line, N lines.
Random orders are driven by `--seed` and are byte-for-byte reproducible under a fixed seed.

> A p-value only says "the observed order is better than N random legal orders". It is **not**
> a posterior probability and does not constitute proof that the ranking is correct.

## 8.8 Randomization type `--random-type`

| Value | Meaning |
|-------|---------|
| `0` (default) | Keep data as-is |
| `1` | Keep nodes, randomize directions |
| `2` | Fully randomize nodes |

The domain is enforced as `{0,1,2}` in argparse: `3`, `-1` and the like are **rejected outright**
with exit code 2, because silently treating them as `0` would amount to running a control
experiment on un-randomized data. `--random-type 1/2` changes the data itself, so "agreement
with the default path"
is meaningless — but the same seed still reproduces. Use them only when a control is needed.

## 8.9 Cooperative cancellation (`stop_check`)

A long local search and a long MCMC chain can be **cancelled cooperatively**. This is an
**API / GUI** capability (there is no CLI switch for it); the parameter is
`stop_check: () -> bool`:

| Insertion point | Signature | Where it is checked |
|-----------------|-----------|---------------------|
| `api.rank(...)` | `stop_check=None` | passed through to `Ranker.run` |
| `Ranker.run(...)` | `stop_check=None` | passed through to the local search and the MCMC chain |
| `ranking.local_search.optimisation_locale(...)` | `stop_check=None` | **first statement of the search loop body**, i.e. at every iteration boundary |
| `robustness.mcmc.MCMCSampler.sample(...)` | `stop_check=None` | between chain steps (in burn-in and in the recording stage alike) |

**What it guarantees**:

- it stops at an **iteration boundary** — never in the middle of an iteration, so the delivered
  solution is always self-consistent;
- it **returns** the best order found as of the cancellation point and writes the three files +
  HTML report as usual: it **does not raise** and does not throw away computed work. On
  cancellation stdout ends with `NOTE: 本次运行按请求**取消**：局部搜索在第 N 次迭代边界停止
  （计划的搜索时长为 T 秒，未跑满）……` (*"this run was cancelled as requested: the local search
  stopped at iteration boundary N, the planned T-second budget was not used up"*);
- whether it really interrupted is machine-readable: `Result.run_metadata["cancelled"]`
  (= the local search's `stats_out["cancelled"]` **or** the chain's `mcmc_stopped_early`),
  `run_metadata["stop_check_requested"]` (was a predicate passed at all) and
  `run_metadata["mcmc_stopped_early"]`;
- **the one exception**: an MCMC chain cancelled **before it recorded any sample** raises
  `RuntimeError` (*"the chain was cancelled (stop_check) before producing its first sample, so
  this run has no approximate-posterior sample to deliver"*) rather than silently returning
  statistics over zero samples — with no samples there is nothing to report, and the program
  refuses to invent it.

**With no `stop_check` (default `None`)** every byte of output is what it always was and
`cancelled` is `False`.

```python
import itertools
from maxtic_next import rank

calls = itertools.count()
def stop_after_50():                    # predicate: turns true after the 50th check
    return next(calls) > 50

r = rank("species.tree", "constraints.tsv",
         local_search=30.0, temperature=60.0, stop_check=stop_after_50)
r.run_metadata["cancelled"]                     # True
r.run_metadata["local_search_iterations"]       # 51 — planned 30 s, stopped at a boundary
r.values["best"]                                # best as of that point, delivered normally
```

Measured (same cyanobacterial instance, `--ls 30`, `--t 60`, predicate turning true at the
51st check): `cancelled=True`, 51 iterations, elapsed far below 30 s, results written.

**Studio's "Cancel ■" button is wired to exactly this path** (Studio's `maxtic_studio/cancellation.py`): button →
`RunEngine.cancel()` → `StopToken.cancel()`; `cancellation_hook` inspects the callee's parameter
names and injects `stop_check=token.should_stop` into `api.rank`, so cancellation really
interrupts the computation at an iteration boundary instead of the previous behaviour of "wait
for it to finish, then throw the result away". If a future core stopped accepting the
parameter, the log panel prints the verbatim `gui/i18n.py` string
`log_cancel_unsupported` — "当前内核版本不接受取消回调，取消只能在阶段边界生效（长 local-search
任务会跑完当前阶段后丢弃结果）。" (*"this kernel version does not accept a cancellation callback;
cancellation only takes effect at stage boundaries — a long local-search task will finish its
current stage and then discard the results"*) — rather than pretending it can interrupt.

> Prev: [07 · Workflows & deployment](07_workflows_deployment.md) ｜ Next: [09 · FAQ & troubleshooting](09_faq_troubleshooting.md)
