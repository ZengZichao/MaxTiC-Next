# 09 · FAQ & troubleshooting

## 9.1 FAQ

**Q1: How is MaxTiC-Next related to the original MaxTiC? Are the results the same?**
A: MaxTiC-Next is a Python 3 faithful rewrite of the original (Python 2), ported
function by function; the reachability-matrix acceleration was checked case by case against
the original `path()` greedy and matched 400/400. But it is **not a blanket "byte-for-byte
equivalent"**: this version deliberately fixes several accounting defects of the original
(the uninformative percentage denominator, the partial-order sentinel threshold, the total
weight after threshold filtering, the write-back of the value after local search, the random-tree
p-value correction, the ALE constraint-source default, …) and adds precondition checks the
original lacked. See the "Intentional deviations from the original MaxTiC" section of the
README for the full list. `--output-style legacy` changes **file names only** (the original
long names); contents are byte-for-byte identical under both styles.

**Q2: What identifies donor/receptor in the constraints?**
A: Internal nodes use **the label in the species tree's bootstrap field** (e.g. `61`); leaves
use leaf names. Constraint endpoints must lie in the species-tree label set.

**Q3: Does `--seed` guarantee identical results across environments?**
A: There are three cases — stop summarizing them as "a fixed seed is reproducible":
- **Default path (`--ls 0`)**: reproducible, byte-for-byte. Constraint line order and
  multi-file order do not affect the products.
- **`--ls > 0`**: the stop condition is **wall-clock time**, so the iteration count follows
  machine load → the same seed **may differ** across machines (the program appends a WARNING to
  stdout saying exactly this).
- **`--ls > 0` + `--local-search-max-iters N`**: the iteration count is fixed → determinism
  restored.

In addition, the greedy heuristic depends on the insertion order of the `edge` dict when
weights **tie**, and Newick parse order feeds into that, so "reproducible" is always qualified
as "fixed seed + same parser + same Python version". Pin versions and containerize.

**Q4: Why did my `--d` do nothing?**
A: `--d` filters on the **phylogenetic distance column only**. Space format needs a 4th
distance column, comma format needs the 5-column ALE layout. If the constraints carry **no**
distance column, `--d>0` **removes nothing at all** (distance-less constraints are always kept,
matching the original) — it is "ignored", not "drops everything". This is stated at
**error level**:

```
ERROR: --min-transfer-distance=5.0 被忽略（123/123 条约束不含 phylogenetic distance 列，……）
```

(the message is emitted in Chinese: *"ignored — 123/123 constraints carry no phylogenetic
distance column"*). Ranking still exits 0; `--dry-run` reports an error and exits 1. Add a
distance column to filter for real, or drop the option.

**Q5: Can local-search output be read as confidence intervals / posterior probabilities?**
A: **No.** The local search is a Metropolis neighbourhood search and yields only a
"robustness/sensitivity summary". `--mcmc` is the one module where the word "distribution" is
permitted — **but it is itself a preliminary implementation that has not passed convergence
diagnostics**: measured ESS is far below the sample count and the adjacent-duplicate rate is
high, i.e. the chain has not mixed. Until those diagnostics pass, **neither** set of numbers
should be published as an uncertainty measure.

**Q6: My small upstream demo families get filtered away — what now?**
A: For small demo samples add `--ale-min-family-size 0`; on real data keep the default `5` or
tune it.

**Q7: Few / zero constraints after an ecceTERA conversion?**
A: ecceTERA labels internal species nodes with its **own numeric node IDs** (leaves by name),
which need not match your species tree's bootstrap labels. Mismatched endpoints are dropped —
and **counted and reported**. Read the stderr diagnostic line (`端点命中率=…`,
`丢弃[标签不匹配=N,…]`) or take
`result.run_metadata["adapter_diagnostics"]`. A hit rate below `--min-endpoint-hit-rate`
(default 0.5) warns; a **0% hit rate is a hard error**. Fix it by aligning your species-tree
internal labels with ecceTERA's naming (or the reverse, i.e. switch to numeric IDs).

**Q8: No HTML report?**
A: The full report needs `jinja2 + plotly` (`pip install "MaxTiC-Next[report]"`). Without them
it does **not** fail: it degrades to a basic template (no interactive plots, and the fallback
template does not HTML-escape). If generation raises, one line is printed to stdout
(`[警告] HTML 报告生成失败，已跳过：…` — "HTML report generation failed, skipped") and ranking
and file output are unaffected. You can also disable it explicitly with `--no-html`.

**Q9: Can I generate constraints only, without ranking?**
A: Yes. `--from <tool> -o constraints.tsv` writes the unified constraints and returns early
(Stage 1 of the two-stage flow). `-o` is **adapter-mode only**: passing it with text constraints
is a hard error instead of a silent no-op, and `--dry-run -o …` does write the file.

**Q10: Can I feed a polytomy (non-binary) species tree?**
A: **No.** The original's `opt`/`mix` hard-assume binarity, and this version enforces that as a
**precondition**: a polytomy raises `ValueError` (exit code 4) instead of silently dropping
nodes and distorting the objective. Resolve the polytomy or re-root first (an unrooted Newick's
root is naturally trifurcating). `--dry-run` reports it as an error too.

**Q11: Why does re-running the same command complain that the output file exists?**
A: That is **overwrite protection**: existing products are not clobbered by default,
and the run aborts with exit code `3` so a previous result cannot be silently destroyed. Pass
`-f/--force` to confirm, or `-p` for another prefix. Products are also written **atomically**
(temporary file + `os.replace`), so a crash cannot leave a half-written file.

**Q12: With several constraint files, why do the products only name the first one?**
A: The prefix comes from the **first** entry of the `CONSTRAINTS` list; the others take part in
the computation but not in the file names. This is declared explicitly on
**stderr** (`NOTE: 输入了 N 个约束文件…`, Chinese for *"N constraint files given; products use
the first as prefix"*) and recorded in `run_metadata["multi_file_prefix_note"]`. To keep products
per source, pass an explicit `-p` for each run.

**Q13: Should `--ale-source` be `rec` or `trf`?**
A: The default is **`trf`** (= `constraints_from_transfers`, the convention documented by ALE's
own MaxTiC integration). The two produce constraint sets that **differ in size and weight
scale**, so switching replaces the entire input. Use `--ale-source rec` only to reproduce results
from versions that defaulted to `rec`. The effective convention is printed to stderr.

**Q14: The upstream example is a `.tgz` — how do I feed it to MaxTiC-Next?**
A: **You do not unpack it first.** `io/compression.py` decides by the gzip magic number
(`1f 8b`) and the `ustar` magic, so `.gz`, extension-less gzip streams and archives with
**exactly one** file member (`.tar.gz`·`.tgz`·`.tar`, renamed included) are decompressed
transparently at **every** entry point (species tree, text constraints, `--dry-run`, all five
adapters, and `--from-auto`'s format sniffing). Family names are taken from the stem after the
compression suffix is stripped, so `.gz` changes no per-family accounting.

The **official example bundles** of the upstream tools can be fed directly: ALE's
`examples/reconciliations.tgz`, for instance — which contains a thousand `*.uml_rec` files —
works **as a whole** in the constraint / upstream-input position, because `api.rank` / the CLI
unpack the archive safely into a process temporary directory and treat **each member** as an
input:

```bash
maxtic-next species.tree examples/reconciliations.tgz --from ale
# stderr: NOTE: 输入 … 是 tar 归档，已自动展开为 N 项输入（解包于 …）
# mapping kept in run_metadata["archives_expanded"]; the prefix comes from the first member,
# so pass -p explicitly if you want per-member products
```

So what does the error mean when a **multi-member** archive is read? It means the archive was
handed over where it would have to be treated as **one** input (typically the species-tree
position, or calling `read_constraints_file` directly). MaxTiC-Next refuses to concatenate the
members on purpose: each member is a separate upstream output (its own gene family / sample
set), and gluing them together would count as a single family, mis-applying every per-family
threshold (`--ale-min-family-size`, `--ale-min-support`, …). That raises
`CompressedArchiveError` (CLI exit code 4; an error-level `--dry-run` entry) whose message
lists the members and gives **directly runnable** commands: `tar -xzf archive -C ./unpacked`
(then feed `./unpacked/*`), or `tar -xzOf archive member > file` to take one member.
Windows 10+ ships `tar`. Details in 05 §5.7.

**Q15: What is `--near-optimal-top-k` for? Should I tune it?**
A: It is the **support-set cap of the robustness/sensitivity summary** — how many **deduplicated**
near-optimal orders the local search keeps for it (API `top_k`, `Ranker.run(top_k=...)`; default
`50`). The summary is not resampling but a deduplicated view of the solutions the search
visited, so orders outside K simply take no part in the counts; measured runs show the chain
walking through nearly ten thousand distinct orders while the summary sees 50. Raise it for a
wider neighbourhood (memory and `summary()` grow like K·n²). The effective value is in
`run_metadata["near_optimal_top_k"]` and in the HTML report, and **must** be published together
with the frequencies. Values must be integers ≥ 1: `--near-optimal-top-k 0` exits 2, and Python
raises `ValueError`. Raising K does **not** make it a posterior (03 §3.6b, 08 §8.1b).

**Q16: I want to stop a long run halfway — is the work so far lost?**
A: No. `api.rank(..., stop_check=...)` / `Ranker.run(..., stop_check=...)` take an
`()`→`bool` predicate that is checked in the local-search main loop and between MCMC chain
steps, so cancellation happens at an **iteration boundary**: the program **returns** the best
order as of that point and writes the products, raising nothing. Whether it really interrupted
is in `run_metadata["cancelled"]` (the local search's `stats_out["cancelled"]` or
`mcmc_stopped_early`). **Studio's "Cancel ■" button is bound to this path**
(the `StopToken` in Studio's `maxtic_studio/cancellation.py`) — the semantics are "deliver the
best order as of the cancellation point", not "wait for the run to finish,
then discard the result". **The one exception**: an MCMC chain cancelled before producing any
sample raises `RuntimeError`, because zero samples is nothing to report. Details in 08 §8.9.

**Q17: Does `自环=N` (*self-loops*) in the adapter diagnostic line mean "N constraints dropped"?**
A: **No** — it is a **counter**. Events with `donor == receptor` have topological distance 0 by
construction; had the adapters written that 0 into the distance column,
`--min-transfer-distance` would have removed them and they would have vanished from the
original `to itself` statistic. **All five adapters — RANGER-DTLx, ecceTERA, ARTra,
AleRax and ALE** — write `distance = None` for self-loops, so under the existing "`--d` keeps
distance-less constraints" rule they survive the filter and `to itself` is not eaten by
`--d` — **which is exactly what the text-input path does** (3-column space / 4-column
comma input has no distance column at all). ALE follows the same rule: a self-loop
arriving through `--from ale` is kept even at the default `--d 0` and shows up in the
summary — measured on a `.uml_rec` whose `T@62->65` is a self-loop under the default `trf`
convention: `1.0 uninformative ( … 1.0 to itself 50.0%)` (03 §3.2.3, measurement B; regression
test `tests/test_ale_selfloop.py`). An event only gets dropped when you yourself write distance 0
— in the text input, or by overwriting the `None` the writer left in the 5th column — and then
stdout says so with `NOTE: --min-transfer-distance=… 丢弃 …` plus
`NOTE: 上述丢弃中有 N 条是**自环约束**…`. See 03 §3.2.3 and 05 §5.2.5.

## 9.2 Typical error messages (verbatim from the implementation)

Messages the program emits in Chinese are annotated with a translation.

| Error / symptom | Exit | Cause | Fix |
|-----------------|------|-------|-----|
| `错误：约束文件 <f> 第 N 行的逗号格式列数不足：'61,62,15.04'` ("comma format has too few columns") | 4 | 3-column comma `donor,receptor,weight` misuse | Use space format, or add the family column to make 4 |
| `错误：约束端点不在物种树中：donor='999', receptor='888'（约束文件/来源：<f>）` ("endpoints not in the species tree") | 4 | donor/receptor labels do not match the tree | Check labels; ecceTERA needs its numeric internal IDs aligned |
| `错误：物种树不是二叉树：内部节点 '1' 有 3 个子节点。` ("not a binary tree") | 4 | polytomy | Resolve it or re-root |
| `错误：物种树存在重复的叶子名（共 1 个名字重复）：'A' 出现在父节点 ['9', '8']` ("duplicate leaf names") | 4 | Non-unique leaf names collapse node identity | Check you did not feed a gene tree / numbered copies as the species tree |
| `错误：物种树内部节点标签不唯一（重复标签：[…]）` ("internal labels not unique") | 4 | Duplicate bootstrap labels | Make internal labels unique |
| `约束文件 <f> 第 N 行的约束权重为负：'-3.0'` ("weight is negative") | 4 | Negative weights void the objective and all percentages | Correct the weights (must be finite and ≥ 0) |
| `约束文件 <f> 第 N 行的约束权重不是有限实数：'nan'` ("not a finite real") | 4 | `nan`/`inf` pollute the statistics | Use real values; delete the line if missing |
| `输出文件已存在：<path>` + "请改用 -p… 或加 --force" ("output file exists") | 3 | Overwrite protection (refused by default) | `-f/--force`, or another `-p` prefix |
| `argument --threshold-constraints/--ts: … 取值域 [0.0, 1.0]，收到 5.0` | 2 | `--ts` out of range (silently accepted before) | Use a fraction in 0–1 |
| `argument --random-type/--r: invalid choice: '3' (choose from '0', '1', '2')` | 2 | Illegal randomization type | Use 0/1/2 |
| `ERROR: --min-transfer-distance=… 被忽略（N/N 条约束不含 phylogenetic distance 列…）` | 0 (1 with `--dry-run`) | `--d>0` without a distance column | Drop `--d` or provide 4/5-column input |
| `[tool] 错误：…` / `LabelMismatchError` (0% adapter endpoint hit) | 4 | Upstream naming inconsistent with the species tree | Align naming, or lower `--min-endpoint-hit-rate` (0 hits still errors) |
| `[error] 上游输出（格式 ranger）未产出任何约束：…` ("upstream output produced no constraints") | 1 (dry-run) | Wrong `--from` tool / family threshold too high / naming mismatch | Read the stderr diagnostics; locate it with `--dry-run` |
| `[error] 上游输出（格式 eccetera）解析失败：UpstreamParseError: …recPhyloXML 解析失败（文件损坏 / 被截断 / 非 XML）…` ("upstream output failed to parse") | 1 (`--dry-run`) | The upstream output itself cannot be read (truncated / not XML / wrong encoding) | The preflight **always** fails and exits 1; re-export or unpack the upstream output |
| `错误：输入文件是一个含 N 个文件成员的 tar 归档（gzip）：<path>` plus `tar -xzf` / `tar -xzOf` commands ("input is a tar archive with N file members") | 4 (error-level entry in `--dry-run`) | A multi-member archive read as **one** input (typically in the species-tree position) | Follow the hint: unpack and feed `./unpacked/*`, or `tar -xzOf archive member > file` for one member. As a constraint/upstream input the archive is expanded automatically (05 §5.7) |
| `argument --near-optimal-top-k: 近优解收集容量 必须 > 0，收到 0` | 2 | Summary support-set cap out of range (must be an integer ≥ 1) | Use an integer ≥ 1; from Python the same rule raises `ValueError` |
| `MaxTiC-Next Studio 需要 GUI 依赖 PySide6（以及 matplotlib）。请在 Studio 仓库根目录执行：pip install -e .` | 3 | `maxtic-studio` run without its GUI dependencies | Install the Studio package (`pip install -e .` in the Studio repo), or use `maxtic-next` (the CLI is unaffected) |
| `WARNING: 输入约束集为空（0 条时间约束）。…` ("input constraint set is empty") | 0 | Data or thresholds left no constraint | **Do not read the 0.0 values as "perfect consistency"** |
| `WARNING: --local-search > 0 以墙钟时长为停止条件…` (wall-clock stop condition) | 0 | Not reproducible across machines | Add `--local-search-max-iters` |
| `--target-clade 指定的节点 … 不在物种树中` ("node not in the species tree") | 4 (error in dry-run) | Wrong clade label | Use a species-tree internal bootstrap label |
| `maxtic-next: command not found` | — | Not pip-installed | `pip install -e .`, or `PYTHONPATH=src python3 -m maxtic_next` |
| HTML report "ran, but looks plain" | 0 | jinja2/plotly absent → silent degradation | `pip install "MaxTiC-Next[report]"` |

## 9.3 Reproducibility checklist

- [ ] Fix `--seed` (default 42)
- [ ] Stay on the default path `--ls 0`; if `--ls > 0` is required, also pin
      `--local-search-max-iters`
- [ ] Pin the Python version (3.11 recommended) and dependency versions
- [ ] Fix `--ale-source` (default `trf`; switching to `rec` replaces the whole constraint set)
- [ ] Use the same Newick parsing path (the bundled `Tree`)
- [ ] Record the full command and parameters (both the HTML report and
      `run_metadata["params"]` carry the run metadata)
- [ ] Avoid `--random-type 1/2` (they change the data itself)
- [ ] Containerize and **pin the image tag** (`maxtic-next:0.1.1`, not `:latest`)
- [ ] When reporting p-values, state the `--random-trees` N and the "(k+1)/(n+1) correction"
- [ ] When reporting a robustness/sensitivity summary, state the effective
      `--near-optimal-top-k` (default 50) — it is the cap on the support set those frequencies
      are computed over (03 §3.6b)
- [ ] If the inputs are compressed: `.gz` and archives can be fed directly, but an expanded
      archive changes the **output prefix** and the members live in a temporary directory
      (see `run_metadata["archives_expanded"]`); pass `-p` for per-member products (05 §5.7)

## 9.4 Performance advice

| Goal | Means |
|------|-------|
| Speed up the local search | incremental scoring is on by default (`--no-incremental` restores the full recompute; same results) |
| Survive interruption of a long search | `--checkpoint` + `--checkpoint-interval` |
| Speed up bulk upstream parsing | `--ale-parallel process` + `--ale-cache-dir` (resumable; key includes the topology fingerprint) |
| Avoid re-parsing | Two stages: freeze `constraints.tsv`, then rerun only the ranking |
| Bound the cost of the randomization test | Choose `--random-trees` deliberately (each draw computes value + similarity) |

## 9.5 Getting help

- Manual overview: [../README.md](../README.md)
- Adapter framework: [../../adapters/README.md](../../adapters/README.md)
- Original tool author: Eric Tannier (eric.tannier@inria.fr); cite MaxTiC, Biorxiv doi.org/10.1101/127548

> Prev: [08 · Advanced features](08_advanced_features.md) ｜ Back: [Manual overview](../README.md)
