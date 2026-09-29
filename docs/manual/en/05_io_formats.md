# 05 · Input / output formats

## 5.1 Input 1: the species tree (Newick)

- A **rooted** Newick tree.
- **All internal node labels must be written in the bootstrap field** (`)label:branchlength`).
  This is MaxTiC's core convention; `donor`/`receptor` in the constraints refer to internal
  nodes through these labels.
- Leaves are referenced by leaf name.
- **It must be binary**: every internal node has exactly 2 children. The original `opt`/`mix`
  hard-assume binarity, and this version enforces that assumption as a **precondition** —
  a polytomy is a hard error (exit code 4) rather than silently dropping nodes.
  `--dry-run` reports it at error level too (see 08 §8.4).
- **Internal node labels must be unique**, and **leaf names must be unique too**:
  a duplicated leaf name collapses node identity.
- Ultrametry is **not required**, but if the tree is ultrametric one can compute a
  Kendall similarity between the input and output rankings from it.
- The file may be **wrapped over several lines**: before parsing, all lines of the file are
  concatenated into one sequence (the original only read the first line; this is a documented
  deviation). Files are read as `utf-8-sig`, so a BOM cannot corrupt the first
  column.
- It may be `.gz` / a gzip stream / a **single-member** `.tar.gz`·`.tgz`·`.tar` archive; no
  unpacking is needed (see 5.7).

Example (internal labels `42/41/59/61/...`; this three-line indented spelling parses as-is):

```
(((((CYAP8:1,CYAP0:1)42:7,(CYAA5:7,UCYNA:7)41:1)59:2,
  (MICAN:9,(CYAP2:3,CYAP7:3)43:6)45:1)61:2,(SYNY3:11,SYNP2:11)56:1)65:1,
 (TRIEI:6,(NOSA0:5,(NOSP7:4,(NOSS1:2,ANAVT:2)40:2)46:1)62:1)67:7)69;
```

## 5.2 Input 2: the constraint file (two formats)

One constraint per line; both formats are **auto-detected** (a line containing a comma is
treated as comma format):

### 5.2.1 Space format

```
donor receptor [weight] [distance]
```

- `weight` is optional, defaults to `1.0` when absent;
- the optional 4th column `distance` is the **phylogenetic distance**, used by
  `--min-transfer-distance`. **Constraints without a distance column are all kept** — i.e.
  `--d > 0` is a complete no-op on pure 3-column input, and an ERROR-level statement of
  convention tells you so ("the option was ignored"; see 03 §3.1).

```
61 62 15.04
61 67 15.48 3.0     # 3.0 is the distance; after the inline # is a comment
```

### 5.2.2 Comma format

```
gene_family,donor,receptor,[weight],[distance]
```

- **the leading `gene_family` column is discarded** (traceability only);
- 4 columns = `family,donor,receptor,weight`;
- 5 columns (ALE) = `family,donor,receptor,weight,distance` (5th column is the distance).

```
fam001,147,149,0.09
fam002,197,187,0.12,2.5
```

### 5.2.3 ⚠️ The forbidden 3-column comma misuse

```
donor,receptor,weight     # wrong! rejected with an error
```

Because the comma format **always discards the first column**, a 3-column line would drop
`donor`, misalign donor/receptor and lose the weight. The parser raises `ValueError` when it
detects a 3-column comma line (exit code 4). If you meant space format, use spaces.

### 5.2.4 Comments, special lines and value validation

- lines starting with `#` are skipped; **inline `#` comments are supported**
  (`61 67 15.48 # note`) — the `#` must be at line start or preceded by whitespace, so that
  lineage names containing `#` are not truncated;
- lines containing `FRQ` are skipped (compatibility with ALE summary lines);
- **weights and distances must be finite reals, and weights must be ≥ 0**: `nan` / `inf` /
  negative weights are rejected **at the parsing layer**, with the **file name and line
  number** in the message (the original carried them straight into the
  objective). Numeric columns accept decimal reals only — not `-`, `NA` or empty strings.
- these files are also read transparently from `.gz` / gzip streams and archives
  (see 5.7).

### 5.2.5 Self-loop constraints and the distance column

A constraint with `donor == receptor` (donor and receptor are the same lineage) has
topological distance 0 by construction, so `--min-transfer-distance` would delete it and it
would disappear from the `to itself` statistic. **All five adapters — RANGER-DTLx, ecceTERA,
ARTra, AleRax and ALE — never emit that 0**: a self-loop is written with `distance = None` (no
distance information) and therefore passes the existing "`--d` keeps distance-less constraints"
rule (regression test `tests/test_ale_selfloop.py`,
see 03 §3.2.3). Consequently `--from ale -o` leaves the 5th column of a self-loop row
**empty** (`selfloop,65,65,1.0,`), and reading that file back yields the same endpoints and
the same statistic — the parser treats an **empty distance column** as "no distance column",
which is the one exception to the "numeric columns reject an empty string" rule of 5.2.4, and
it applies to the distance column only.
**The text-input path always behaved this way** — 3-column space / 4-column comma input carries
no distance column, so self-loops were always kept. Only if **you** write `0.0` in the 4th
(space) or 5th (comma) column does the event get dropped, and then stdout appends
`NOTE: 上述丢弃中有 N 条是**自环约束**…` (*"N of those were self-loop constraints…"*)
saying what `to itself` lost. See 03 §3.2.3.

## 5.3 Output: files

The prefix defaults to the path of the **first constraint file** (CLI `-p/--output-prefix`,
API `output_prefix`). When several constraint files are given, the others take part in the
computation but do not appear in the names, and this is declared once on **stderr**
(see 03 §3.5). `--output-style` selects the naming style:

| Contents | short (default) | legacy (original long names) |
|----------|-----------------|------------------------------|
| Filtered weighted informative constraints | `<prefix>.mt.informative.tsv` | `<prefix>_MT_output_filtered_list_of_weighted_informative_constraints` |
| Constraints conflicting with the best order | `<prefix>.mt.conflicts.tsv` | `<prefix>_MT_output_list_of_constraints_conflicting_with_best_order` |
| Partial order | `<prefix>.mt.partial_order.tsv` | `<prefix>_MT_output_partial_order` |
| Random-order distribution (**only `--random-trees > 0`**) | `<prefix>.mt.random_dist.tsv` | `<prefix>_distribution_random` |

> - Under both styles the file **contents are byte-for-byte identical**; only names differ.
> - **Overwrite protection**: if a target file exists and was not written by this run, the
>   CLI aborts with exit code **3**; only an explicit `-f/--force` overwrites.
> - **Atomic writes**: a temporary file in the target directory is written, `fsync`-ed and
>   then moved into place with `os.replace`, so a crash cannot leave a half-written product.
>   Missing parent directories are created automatically.

### File contents

- **informative**: one `donor,receptor weight` per line — the weighted informative
  constraints after uninformative ones are filtered out. Hard constraints that agreed with a
  phylogenetic edge and were set to `MAX_NUMBER` (1e10), and 0-weight edges, do not appear.
- **conflicts**: `donor,receptor weight` — constraints **violated** by the best ranking.
- **partial_order**: `a b weight color` — constraints **consistent** with the best ranking
  (they form the partial order);
  - `black`: same direction as the input tree's ranking;
  - `green`: opposite direction (i.e. flipped relative to the input).
  - The sentinel-edge test uses `MAX_NUMBER` (1e10) rather than the hard-coded `100000` of
    the original, so informative constraints with weight ≥ 1e5 are not
    silently removed.
- **random_dist**: one `value similarity` pair per line, `--random-trees` lines in total —
  the empirical null distribution behind the p-values.

## 5.4 Output: the stdout summary

The print order and wording of the original `MaxTiC.py` are reproduced line by line (full
worked example in [02 §2.3](02_quickstart.md)): internal node count, total constraint
weight, uninformative breakdown and share, trivially-conflicting share, the value and
percentage of the input/greedy/mixing rankings, the best-order identifier, the ranked-tree
Newick, and the Kendall similarity with the input tree.

Lines appended depending on options:

| Trigger | Appended |
|---------|----------|
| `--ls > 0` | an `attempting a local search …` notice; and after `best order is the …` the two lines `after local search X rejected` and `best found solution Y (%)`, so stdout cannot contradict the delivered tree |
| `--rd N > 0` | 4 lines: random value range, value p-value, random similarity range, similarity p-value. p-values use the **(k+1)/(n+1)** correction and test the **final delivered order** |
| `--mcmc` | 5 `[MCMC]` lines: status statement, convention statement, effective temperature with burn-in/thin, sample statistics, convergence diagnostics |
| whenever applicable | `NOTE:` / `WARNING:` / `ERROR:` statements of convention and warnings appended at the **end** of the summary (phylogenetic hard constraints, empty constraint set, `--d` ignored, self-loops removed by the distance filter, `--ls` not reproducible across machines, …) |
| **cancelled** run (API/GUI `stop_check`) | `NOTE: 本次运行按请求**取消**：局部搜索在第 N 次迭代边界停止（计划的搜索时长为 T 秒，未跑满）……` (*"this run was cancelled as requested: the local search stopped at iteration boundary N (the planned budget of T seconds was not used up)"*) — what is delivered is the best order **as of the cancellation point** (see 08 §8.9) |

Percentage conventions (deviations from the original, documented deliberately):

- the `uninformative` percentage = `uninformative ÷ total_weight`. The original denominator
  `total_transfers + uninformative` **double-counts** the uninformative weight and prints
  `12%` on the same instance; this version prints `13.9%`.
- when `--threshold-constraints` is active, the total weight is computed **after** edges are
  removed (the original sums **before** removal, so deleted weight stays in the denominator).

## 5.5 Output: the interactive HTML report (default)

- A single fully static file `<prefix>.html`, no backend, opens offline.
- Contents: run metadata (seed/parameters/version/adapter diagnostics), constraint weight
  distribution, conflict share, node-position robustness summary (when local search is on),
  and a visualization of the ranked tree.
- Disable with `--no-html`.
- Without `jinja2` / `plotly` installed (`pip install "MaxTiC-Next[report]"`) it does **not**
  fail: it degrades to a basic template — no interactive plots, and the fallback template
  does **not** HTML-escape. Install the extras for the full report.

## 5.6 Output: the robustness/sensitivity summary (when local search is on)

Presented in `Result.sensitivity_summary` (Python) and in the HTML report: top-K near-optimal
orders, per-node position distribution, pairwise order frequencies. It is called **only** a
"robustness/sensitivity summary" and never uses "confidence interval / posterior probability"
(that vocabulary belongs to MCMC, which itself has not passed convergence diagnostics —
see [chapter 08](08_advanced_features.md)).

It is a **deduplicated view of the set of solutions the local search visited**, not
robustness under resampling (`is_resampling_robustness` is always `False`), and the support of
those frequencies is capped by `--near-optimal-top-k` (API `top_k`): by default **50**, i.e.
the summary considers at most the 50 best deduplicated near-optimal orders, while the search
typically walked through far more distinct ones (`n_unique_orders_total`). The effective value
is recorded in `run_metadata["near_optimal_top_k"]`, and the summary dict carries `top_k` /
`top_k_reported`; report K alongside any of these frequencies. See 03 §3.6b and 08 §8.1b.

## 5.7 Compressed input: `.gz` / gzip streams and tar archives (transparent)

Upstream tools ship their examples and results **compressed** (MaxTiC's
`examples/reconciliations.tgz`, ALE/Ranger/AleRax outputs bundled as `.tar.gz`, single-file
examples as `xxx.gz`). `io/compression.py` is the project's **only** text-input reading layer,
so **nothing has to be unpacked first**. Detection uses *both* the content magic number and the
file extension, which avoids the failure mode each one alone has:

| Input | Decided by | Behaviour |
|-------|-----------|-----------|
| `x.gz` / any file whose first two bytes are `1f 8b` | gzip magic (RFC 1952) | decompressed transparently; **works with no extension at all** |
| extension says `.gz` but the content is not gzip | magic mismatch | read as **plain text** (a misnamed uncompressed file still works; no more `BadGzipFile`) |
| `.tar.gz` / `.tgz` / `.tar` with **exactly one** file member (renamed ones included) | `ustar` magic at offset 257 + `tarfile` | that member's content is used, exactly like a plain text file |
| tar archive with **several** file members | as above | **never concatenated**: raises `CompressedArchiveError` (a `ValueError` subclass; CLI exit code 4) whose message lists the member count, some member names and a **directly runnable** `tar -xzf archive -C dir` / `tar -xzOf archive member > file` command |

### 5.7.1 Which entry points actually take it

- **Species tree (5.1)**: `.gz` / gzip streams / **single-member** archives all work (through
  `io.parsing.read_newick_file`). A multi-member archive **is** an error here — one tree can
  only be one member.
- **Constraint files (5.2)**: same, via `io.parsing.read_constraints_file` and the convenience
  wrapper `constraints.parsers.parse_constraints_file`.
- **`--dry-run` preflight**: the same reading layer, so compressed inputs are transparent in
  the preflight too; a multi-member archive becomes an **error-level** report entry rather
  than a raw traceback.
- **Upstream adapters**: the per-file parsers of ALE / RANGER-DTLx / ecceTERA / ARTra / AleRax
  all go through this layer, so `FAM1.dtl.gz` and `fam.uml_rec.gz` can be handed to
  `--from <tool>` directly; `--from-auto`'s format sniffing reads the **decompressed** text as
  well (raw compressed bytes would otherwise be reported as "unrecognized").
  Family names come from the stem **after** the compression suffix is stripped
  (`gene_1.modif.ale.uml_rec.gz` and `gene_1.modif.ale.uml_rec` yield the **same**
  `metadata["family"]`), so `.gz` does not change any per-family accounting.
- **Archives given as input** (`CONSTRAINTS` positional / the constraint list passed to
  `api.rank`): `api.rank` first **safely unpacks** `.tgz` / `.tar.gz` / `.tar` into a
  process-local temporary directory and then treats the **members** as the inputs. A
  multi-member archive — e.g. ALE's official `examples/reconciliations.tgz`, which holds a
  thousand `*.uml_rec` files — can therefore be fed **whole** in the constraint position.
  One `NOTE: 输入 … 是 tar 归档，已自动展开为 N 项输入（解包于 …）` line goes to **stderr**
  and the mapping is stored in `run_metadata["archives_expanded"]` (the members live in a
  temporary directory, so this must stay traceable). An AleRax output-tree archive is
  expanded to its **root directory**, which the AleRax adapter then reads by its own layout
  convention. Unpacking is sandboxed: only regular files are extracted, and absolute paths,
  members containing `..`, and link/device/FIFO members are skipped.

> **Do not overstate it**: concatenating a multi-member archive into **one** text input is
> never done — each member is a separate upstream output (its own gene family / sample set),
> and gluing them together would be treated as a single family, which would miscompute every
> per-family threshold (`--ale-min-family-size`, `--ale-min-support`, …). A multi-member
> archive is only auto-expanded where "each member is its own input" makes sense, i.e. the
> constraint / upstream-output positions.

### 5.7.2 Measured

```bash
# .gz constraints — products byte-for-byte identical to the uncompressed input
maxtic-next minitree.tree constraints.tsv.gz --seed 42          # exit 0
diff constraints.mt.informative.tsv constraints_gz.mt.informative.tsv   # no difference

# single-member .tar.gz (constraints) and .tgz (species tree) both work
maxtic-next minitree.tree single.tar.gz --seed 42               # exit 0
maxtic-next single.tgz constraints.tsv --seed 42                # exit 0

# an official multi-member .tgz handed over as upstream input
maxtic-next species.tree reconciliations.tgz --from ale         # expanded into N .uml_rec files

# a multi-member archive read as ONE input (here: the species-tree position) → exit 4 + tar commands
maxtic-next multi.tar.gz constraints.tsv --seed 42; echo "exit=$?"
错误：输入文件是一个含 2 个文件成员的 tar 归档（gzip）：/tmp/multi.tar.gz
  成员：Cyano_CUTConstraints.tsv, minitree.tree
  说明：MaxTiC-Next 透明解压**单文件** gzip（.gz / 无扩展名的 gzip 流）与**恰好一个成员**的
  tar 归档；多成员归档不会被自动拼接成一个输入，……
  请任选其一，把成员变成真正的输入：
    1) 解包整个归档，再把解出的文件（或目录，AleRax 支持目录）作为输入：
       tar -xzf /tmp/multi.tar.gz -C ./unpacked && ls ./unpacked
       maxtic-next species_tree ./unpacked/* [--from auto]
    2) 只取其中一个成员：
       tar -xzOf /tmp/multi.tar.gz Cyano_CUTConstraints.tsv > Cyano_CUTConstraints.tsv
  （Windows 10+ 自带 tar；--dry-run 会原样复述本提示）
exit=4
```

This module depends only on the standard library (`gzip` / `tarfile`); no third-party
decompression dependency is introduced.

> Prev: [04 · Python API](04_python_api.md) ｜ Next: [06 · Upstream integration](06_upstream_integration.md)
