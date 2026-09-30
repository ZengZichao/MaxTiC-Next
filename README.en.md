# MaxTiC-Next

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23053051.svg)](https://doi.org/10.5281/zenodo.23053051)

> 🇨🇳 中文版：[README.md](README.md) | 📖 Full bilingual manual: [docs/manual/README.md](docs/manual/README.md)

`MaxTiC-Next` is a Python 3 rewrite and engineering extension of **MaxTiC**
(Eric Tannier, Inria), a phylogenetics tool from the Python 2 era. On top of a
**functionally equivalent** port of the original ranking algorithm under controlled conditions
(fixed seed, fixed Python version, fixed parser), it provides a modern CLI/Python API, a
pluggable adapter framework for **5 upstream tools** (ALE / RANGER-DTLx / ecceTERA / ARTra /
AleRax), robustness/sensitivity summaries plus an optional MCMC sampler (preliminary), an
interactive HTML report, and complete packaging and workflow integration
(Docker/Singularity, Snakemake/Nextflow).

> **Documentation navigation**: the full bilingual manual — installation / usage /
> input-output / Python API / upstream integration / workflows & deployment / advanced
> features / FAQ — lives in [`docs/manual/`](docs/manual/README.md) (10 chapters × 2 languages).

## Algorithmic equivalence

A core design goal is a **functionally equivalent** port of the original algorithm under
controlled conditions (fixed seed, fixed Python version, fixed parser):

- `path()` / `value()` / `edgeweights()` / `mix()` / `opt()` / `order_from_graph()` /
  `optimisation_locale()` are ported line by line from `MaxTiC.py`;
- all `dict.keys()/values()/items()` views are wrapped in `list()` before `.sort()` / `del`;
  integer division `/` → `//`; `print` becomes a function; `cmp` → `functools.cmp_to_key`;
- all randomness flows through a single `random.Random(seed)` instance
  (`maxtic_next.random_.RandomWrapper`), default `--seed 42`;
- `--min-transfer-distance` (originally `d`) filters on the **phylogenetic distance** column
  (keeping `distance > threshold`) and is **never** repurposed as a weight filter; constraints
  without a distance column are always kept (the option is then a no-op, and the program says
  so explicitly);
- the stdout summary and file contents reproduce the line order and wording of `MaxTiC.py`.

### The **exact boundary** of equivalence (do not quote this as "byte-for-byte")

1. **Accounting fixes**: this version deliberately fixes several statistical/accounting
   defects of the original (see "Intentional deviations" below), so stdout and the products are
   **not** identical to the original's. Any "byte-for-byte equivalent" claim must be qualified
   as *"byte-for-byte equivalent except for the documented deviations"*.
2. **Reproducibility**: **the default path (`--local-search 0`) is byte-for-byte reproducible
   under a fixed `--seed`**. With `--ls > 0` the stop condition is **wall-clock time**, so the
   iteration count follows machine load and the same seed is **not guaranteed to reproduce
   across machines** (the program appends a WARNING to stdout saying exactly this); for
   determinism also set `--local-search-max-iters`. In addition the greedy heuristic depends on
   `edge` dict insertion order **when weights tie**, so "reproducible" is always qualified as
   "fixed seed + same parser + same Python version".
3. **Two different kinds of evidence** (do not conflate them):
   - `tests/baselines/` (a `stdout.txt` snapshot plus the three product files) are the
     **self-produced deterministic output snapshots of this port at `seed=42`**. Their role is
     **regression prevention** — the same code must produce the same snapshot. They are **not**
     original MaxTiC output and they do **not** prove cross-implementation equivalence;
   - **the actual cross-implementation gate** is `tests/test_reference_equivalence.py`: an
     **independently implemented** Python 3 reference implementation following the control flow and
     formulas of the original `MaxTiC.py` (sharing no code with the tested
     `ReachabilityMatrix` / `EdgeBuilder` / `order_from_graph` / `mix` / `opt` / `value` /
     `RandomWrapper`), diffed **field by field** against `maxtic_next` on the real cyanobacterial
     example. Any drift in the core algorithm or in the statistics makes it fail. It also
     **explicitly asserts** the documented deviations (the uninformative denominator, the
     partial-order sentinel, the post-threshold total weight, the non-deletable sentinel edges;
     and it records that the exclusion of zero-weight edges from the informative file matches the
     original and is therefore *not* a deviation). Newick parsing
     is out of its scope (that file reuses this project's `Tree`, as noted in it).
   - The "run the original and compare" case in `tests/test_equivalence_stats.py` needs both
     `python2` and the original `MaxTiC.py` in the environment; it skips when they are absent and
     the differential test above stays active as the fallback.

## Installation

```bash
pip install -e .            # basic install (core functionality, standard library only)
pip install -e ".[report]"  # full install (adds jinja2/plotly for the interactive HTML report)
```

> The console command is lowercase `maxtic-next` (see `[project.scripts]` in pyproject).
> This repository ships **only** the command line and the Python API — no GUI code, no GUI
> dependencies.

## Where is the graphical front end?

The desktop app **MaxTiC-Next Studio** (PySide6, bilingual zh/en + light/dark themes) lives in
a separate repository, **[`MaxTiC-Next-Studio`](https://github.com/ZengZichao/MaxTiC-Next-Studio)**, with its own installation, packaging
and release cycle:

```bash
pip install .                # this repo (core algorithm)
pip install "git+https://github.com/ZengZichao/MaxTiC-Next-Studio.git"  # Studio, which provides maxtic-studio
maxtic-studio
```

The dependency runs one way: Studio declares `MaxTiC-Next>=0.1.0` in its `pyproject.toml`, while
this repository never references Studio — so where you clone either one, and what you name the
local folders, is irrelevant.

An install-free macOS bundle, `MaxTiC-Next-Studio.app` (frozen with PyInstaller — just
double-click), is produced by that repository's `packaging/build_studio_app.sh` into its
`release/` directory; its usage manual is `docs/studio.en.md` there. Both repositories share the
same algorithm: Studio calls this repository's `maxtic_next.api.rank` in-process and rewrites
no computation.

Offline / without installing, from the repository root:

```bash
export PYTHONPATH=src
python -m maxtic_next --help
python -m maxtic_next --version   # MaxTiC-Next 0.1.0
```

## Usage

```bash
# run with the standard library only (no network install needed)
PYTHONPATH=src python3 -m maxtic_next examples/minitree.tree examples/Cyano_CUTConstraints.tsv --seed 42

# or after installing
maxtic-next examples/minitree.tree examples/Cyano_CUTConstraints.tsv --seed 42
```

### Options (core)

| Option | Description | Default |
|--------|-------------|---------|
| `species_tree` | Species tree (Newick, internal labels in the bootstrap field; must be **binary**, with unique internal labels and unique leaf names; `.gz` / gzip streams / archives are read transparently) | required |
| `constraints` | One or more constraint files (space `donor receptor [weight] [distance]` or comma `family,donor,receptor,[weight],[distance]`); likewise transparent about `.gz` / archives, and official upstream `.tgz` example bundles can be fed whole | required |
| `--seed` | Random seed (drives `mix` tie-breaking and the local search) | `42` |
| `--local-search` / `--ls` | Local-search **wall-clock budget** in seconds, 0 = off | `0` |
| `--local-search-max-iters` | Local-search iteration cap (0 = bounded by time only); set positive for cross-machine determinism | `0` |
| `--temperature` / `--t` | Metropolis temperature | `0.001` |
| `--random-type` / `--r` | Randomization type; only 0/1/2 accepted | `0` |
| `--min-transfer-distance` / `--d` | Minimum transfer distance threshold (filters on the phylogenetic distance column) | `0` |
| `--threshold-constraints` / `--ts` | Constraint weight threshold fraction, domain [0, 1] | `0.0` |
| `--random-trees` / `--rd` | Random-tree sample count (p-values use the (k+1)/(n+1) correction) | `0` |
| `--near-optimal-top-k` | Near-optimal collector capacity: how many **deduplicated** near-optimal orders the robustness/sensitivity summary may consider (API `top_k`; must be ≥ 1; only meaningful with `--ls > 0`) | `50` |
| `-f` / `--force` | Allow overwriting existing products (default: refuse, exit code 3) | off |

> For the complete set (upstream adaptation via `--from` and `--min-endpoint-hit-rate`, clade
> pruning, `--dry-run`, `--mcmc*` including `--mcmc-burn-in` / `--mcmc-thin` / `auto`
> temperature, `--incremental`, checkpointing, `--output-style`, `-p/--output-prefix`, …) see the
> [CLI reference manual](docs/manual/en/03_cli_reference.md).

### Python API

```python
from maxtic_next import rank, build_constraints

result = rank("examples/minitree.tree", "examples/Cyano_CUTConstraints.tsv", seed=42)
print(result.best_source, result.similarity_to_input)   # greedy heuristic 1.0
```

## Upstream tool integration

MaxTiC-Next ships adapters for the five upstream reconciliation / transfer-inference tools —
**ALE / RANGER-DTLx / ecceTERA / ARTra / AleRax** — unified into weighted constraints through a
plugin registry plus auto-detection:

```bash
maxtic-next species.tree FAM1.dtl        --from ranger
maxtic-next species.tree out.recphyloxml --from-auto
maxtic-next species.tree alerax_run/     --from alerax
```

- **ALE's `--ale-source` defaults to `trf`** (transfer events, `parent(donor) -> receptor`), i.e.
  the convention of ALE's own MaxTiC integration (`constraints_from_transfers`). `rec`
  (reconciliation events) remains selectable, but the two produce constraint sets of different
  size and scale, so switching replaces the entire input.
- Every conversion prints one diagnostic line to **stderr** (endpoint convention, weight
  semantics, sample denominator, endpoint hit rate and **per-reason drop counters**) and records
  it in `result.run_metadata["adapter_diagnostics"]`; a 0% endpoint hit rate is a hard error, and
  anything below `--min-endpoint-hit-rate` (default 0.5) warns. `--quiet-adapters` suppresses the
  printing only.

See the [upstream integration manual](docs/manual/en/06_upstream_integration.md) and the
[`docs/adapters/`](docs/adapters/README.md) specifications.

## Output

The output prefix defaults to the path of the **first constraint file** (multi-file runs declare
this on stderr; override with `-p`). Default **short** naming (concise, with `.tsv` extensions):

- `<constraints>.mt.informative.tsv` — filtered weighted informative constraints
- `<constraints>.mt.conflicts.tsv` — constraints conflicting with the best order
- `<constraints>.mt.partial_order.tsv` — partial order (with black/green flags)
- `<constraints>.mt.random_dist.tsv` — the `value similarity` distribution of the random orders,
  written **only with `--random-trees > 0`** (the 4th data product)

`--output-style legacy` switches back to original MaxTiC long names (handy when comparing by
filename):

- `<constraints>_MT_output_filtered_list_of_weighted_informative_constraints`
- `<constraints>_MT_output_list_of_constraints_conflicting_with_best_order`
- `<constraints>_MT_output_partial_order`
- `<constraints>_distribution_random`

> Under both styles the **contents are byte-for-byte identical**; only names differ, and stdout is
> unaffected by the style.
> Products are written **atomically** (temporary file + `os.replace`) and overwriting existing
> files is **refused by default** (exit code 3; needs `-f/--force`).

An interactive HTML report is produced by default (`--no-html` disables it; without
`jinja2`/`plotly` it **degrades** to a basic template rather than failing). The stdout summary
covers the internal node count, total constraint weight, uninformative breakdown and share,
trivially conflicting share, each heuristic's value, the best order, the ranked-tree Newick and
the Kendall similarity with the input tree, followed by runtime `NOTE:` / `WARNING:` / `ERROR:`
statements of convention.

### Exit codes

| Code | Meaning |
|------|---------|
| `0` | Success (including a passed `--dry-run`) |
| `1` | `--dry-run` found error-level problems (in `--from` mode, "the upstream output cannot be parsed at all" lands here too) |
| `2` | argparse domain error (including `--near-optimal-top-k < 1`) |
| `3` | Output exists and `--force` was not given |
| `4` | Input fails an algorithmic precondition (non-binary tree, duplicate labels, missing endpoints, illegal weights, 0% adapter hit rate, a multi-member archive read as one input, …) |

## Intentional deviations from the original MaxTiC

Every deviation below is **deliberate** and answers a specific defect of the original or an
accounting decision of this version. They are also spelled out in the declarations of
`tests/test_reference_equivalence.py` and `tests/test_equivalence.py`, where they are asserted
item by item rather than hidden inside an "everything is equal" assertion.

### Fixes of original defects (improvements)

| Item | Original behaviour | This version |
|------|--------------------|--------------|
| Self-loop constraint `X X w` | **infinite loop** (`MaxTiC.py:201-204`) | `v != current` guard in `greedy.py` + counted under `to_itself`; verified not to hang |
| Kendall on a single-internal-node tree | `ZeroDivisionError` (`MaxTiC.py:138`) | returns `1.0`, with an explanatory comment |
| Label uniqueness / endpoint existence | bare `KeyError` | contextual `ValueError` (file name, edge key, suggested fix) |
| `d=MIN_TRANSFER_DIST` | `int(words[1])` **truncates** the distance (`MaxTiC.py:45`) | `float` |
| Threshold removing everything | `IndexError` (`MaxTiC.py:483`) | guard + explicit `WARNING`: a 0.0 value means "nothing to test", not "perfect consistency" |
| 2/3-column comma misuse | `IndexError` / misaligned read | `parsers.py` rejects it with an actionable message |
| Random source | global `random` (unseeded) | all randomness driven by a single `RandomWrapper(seed)` |
| Cross-source auto-detection | not available | added; mixed input is grouped per tool and warned about, never silently summed |
| `uninformative` percentage denominator | `total_transfers + uninformative` (**double-counts** the uninformative weight) → prints `12%` on the example | `total_weight` → prints `13.9%` |
| Partial-order sentinel test | hard-coded `edge[e] < 100000` (its `MAX_NUMBER` is `1e10`) → silently deletes real informative constraints with weight ≥ 1e5 | tests against `MAX_NUMBER` |
| `--threshold-constraints` total-weight basis | summed **before** removing edges, so deleted weight stayed in the denominator | computed **after** removal |
| Threshold deletion set | could delete the `MAX_NUMBER` phylogenetic hard constraints too → produced a "legal" ranking no longer constrained by the species tree | sentinel edges are excluded from the deletable set |
| stdout after local search | printed the **pre-search** value and source, contradicting the delivered tree | writes back `values` / `best_source` and adds the `after local search … rejected` / `best found solution …` lines |
| `--random-trees` p-value | `k/n`, squeezing the smallest non-zero p to `0.0` at n=50; and it tested the **pre-search** order | `(k+1)/(n+1)` correction (Phipson & Smyth 2010), testing the **final delivered order** |
| Unvalidated option domains | `--r 3`, `--ts 5`, `--ls -5` silently accepted | domains enforced in argparse; violations exit 2 |
| Wrapped Newick, BOM, inline `#` comments | first line only / BOM leaks / `float('#')` crashes | whole file concatenated, read as `utf-8-sig`, inline comments supported; errors carry file name and line number |
| Illegal weights | `nan`/`inf`/negative values carried into the objective | rejected at the parsing layer |
| Non-binary species tree | unchecked (polytomies silently dropped) | enforced as a **precondition**; error, exit code 4 |
| Duplicate leaf names | unchecked (node identity collapses) | error naming the duplicate and its parents |
| ALE constraint-source default | ALE's own MaxTiC integration specifies `constraints_from_transfers` (trf) | default **`trf`** (a `rec` default silently replaces the whole constraint set) |
| ALE file-level cache key | — | key includes a **species-tree topology fingerprint** (two trees with the same labels but different topology cannot contaminate each other) |
| Silent adapter drops | — | every drop path is **counted and reported** (one stderr line + `run_metadata`); a 0% endpoint hit rate is a hard error |
| Output overwrite and atomicity | `open(...,"w")` clobbers, can leave half-written files | **refuses to overwrite by default** (`--force` opts in) + temporary file + `fsync` + `os.replace` |
| Self-loop events (donor == receptor) removed by the distance filter | their topological distance is always 0, so `d` ate them and `to itself` silently under-counted | **all five adapters — RANGER-DTLx, ecceTERA, ARTra, AleRax and ALE — write `distance = None`** for self-loops, so the "distance-less constraints are kept" rule lets them through and `to itself` is not eaten by `--d`; the text-input path has always behaved this way. All five adapters share this convention: regression test `tests/test_ale_selfloop.py`, measurement in manual 03 §3.2.3 |
| Compressed input | plain text only | `.gz` / gzip streams decompressed by magic number and single-member tar archives by the `ustar` magic; official upstream `.tgz` example bundles can be fed whole (expanded into their members). A multi-member archive is **never** concatenated into "one" input — it raises with runnable `tar -xzf` / `tar -xzOf` commands |
| `--dry-run` on upstream output | — (the original has no preflight) | the preflight **actually parses** it: corrupt/truncated XML, zero constraints and 0% endpoint hits are all error-level, so the CLI **exits 1** |
| Cancelling the GUI front end | — (the original has no GUI) | the Cancel button interrupts at an **iteration boundary** via `stop_check` and delivers the best order so far; it does not mean "wait for the run to finish, then discard it" |

### Terminology discipline (not a deviation, but you must know it)

- Local-search output is called **only** a "robustness/sensitivity summary over solutions visited
  by the local search"; "confidence interval / posterior probability" is never used.
- That summary is a **deduplicated view of the local search's visit set**
  (`is_resampling_robustness` is always `False`), and its support is capped by
  `--near-optimal-top-k` (API `top_k`, default **50**): it considers at most the 50 best
  deduplicated near-optimal orders even though the chain visited far more distinct ones. Raising
  K widens the support, it does **not** make this a posterior; whenever you report those
  frequencies you **must** state the effective K
  (`run_metadata["near_optimal_top_k"]`).
- `--mcmc` is a **preliminary implementation**: convergence diagnostics are not validated
  (measured ESS far below the sample count, high adjacent-duplicate rate), so its samples **must
  not** be treated as posterior samples. Every line it prints carries that qualification.

## License

Inherits the original **CeCILL 2.1**, retaining author credit to Eric Tannier and the citation
(Biorxiv doi.org/10.1101/127548). See `LICENSE`.

## Current capabilities

- **Core rewrite**: faithful port of the ranking algorithm (greedy / mixing / local search), with
  the accounting defects fixed and a differential test as the gate.
- **Interfaces**: the `maxtic-next` command line and the `maxtic_next` Python API (the desktop
  front end is the separate [`MaxTiC-Next-Studio`](https://github.com/ZengZichao/MaxTiC-Next-Studio) repository, built on the same API).
- **Upstream adapters**: ALE / RANGER-DTLx / ecceTERA / ARTra / AleRax via plugin registry +
  auto-detection, with endpoint-hit-rate and drop-counter diagnostics.
- **Uncertainty**: local-search robustness/sensitivity summary (its support set capped by
  `--near-optimal-top-k`, default 50) + an optional MH sampler (preliminary; terminology kept
  strictly separated).
- **Reporting**: interactive HTML report with no backend dependency (degrades rather than fails).
- **Input**: constraints and species trees are read transparently from `.gz` / gzip streams /
  tar archives, and official upstream `.tgz` example bundles can be fed whole.
- **Long-run task control**: cooperative cancellation (`api.rank(stop_check=...)`, which is what
  Studio's Cancel button drives), a local-search iteration cap, and checkpoint resumption.
- **Analysis aids**: clade pruning, `--dry-run` preflight, random-tree p-values with the
  (k+1)/(n+1) correction, incremental scoring, checkpointing.
- **Deployment**: Docker/Singularity images and Snakemake/Nextflow wrappers, with parallel
  constraint generation.

## Performance notes (what is actually wired)

- **Acceleration always on**: `ReachabilityMatrix` in `ranking/reachability.py` (incremental
  transitive closure, O(1) lookups instead of O(V) traversal) is active **in every ranking**, and
  was shown to agree with the original `path()` greedy on 400/400 cases. It is **on by default**,
  not an optional switch.
- **Incremental scoring on by default**: `--incremental` (each interval rotation in the
  local search costs O(b−a) instead of a full O(|E|) recompute). It is the default because the
  incremental and full paths are verified **byte-for-byte identical** over 4 inputs × 5 seeds ×
  {stdout summary, ranked newick, best_order, values, the three TSVs}, while being 3.1–8.2×
  faster; `--no-incremental` keeps the original full-recompute semantics as an escape hatch and
  the equivalence is guarded long-term by `tests/test_output_equivalence.py`.
- **Faster default path** (no change to summation order): prefix membership tests in the `mix`
  dynamic program use incrementally maintained sets instead of list slices, cutting the
  synthetic benchmark point (n=3199, |E|=12796) from a 21.7–32.9 s median to **3.7–7.6 s** per run
  (medians of two rounds on this host; round-to-round noise on the same input is up to ~2x. The
  input is a deterministic synthetic tree, not a real biological dataset).
- **Optional switches**: `--checkpoint` (resumption covering the `--local-search` stage only), and
  ALE's `--ale-cache-dir` (per-family cache for the parsing stage).
- `PathCache` in `tree/cache.py` is a legacy implementation **never wired into the ranking hot
  path**, kept exported only so external imports do not break. It is **not** "enabled by
  `--incremental` / `--checkpoint`" — those two options refer to the real paths above. This is
  stated in that module's own docstring.

## Future work

- Under strict equivalence-preserving regression control, evaluate whether wiring tree-topology
  caching of the `PathCache` kind into the hot path is worth it at all.
- Broader real-sample coverage of upstream tools; containerized cross-version equivalence
  verification (Python 2 original vs this version) in CI — currently the independent
  reference implementation in `tests/test_reference_equivalence.py` serves as that gate.
- Convergence diagnostics and adaptive temperature for MCMC (until those pass, `--mcmc` keeps its
  "preliminary" qualification).

## Known limitations

- **`--ls > 0` is not reproducible across machines** (unless `--local-search-max-iters` is also
  set): its stop condition is wall-clock time.
- **`--checkpoint` covers only the local-search stage**: checkpointing applies to the Metropolis
  search driven by `--local-search`; the ALE parsing stage has its own per-gene-family file-level
  cache and is not subject to this limitation.
- **The species tree must be binary and uniquely labelled**: this version does not generalize
  `opt`/`maximum_distance`/`random_order` to multifurcating trees; it enforces binarity as a
  precondition instead (see "Intentional deviations").
- **Weight scales across adapters are not commensurable**: ALE gives [0,1] support; RANGER-DTLx /
  ARTra give **integer counts** when the report declares no sample count; AleRax gives
  upstream-preaggregated posterior frequencies. When mixing them, `--ale-min-support` changes
  meaning accordingly, and the diagnostic line tells you which case you are in.
