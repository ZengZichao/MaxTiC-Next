# 04 · Python API

MaxTiC-Next exposes a compact Python API for embedding in scripts and pipelines. The core
exports live at package top level:

```python
from maxtic_next import rank, build_constraints
```

## 4.1 `rank(...)` — end-to-end ranking

The signature below is transcribed verbatim from `src/maxtic_next/api.py` (default constants
resolve through `config.py`):

```python
def rank(species_tree_path, constraints_path,
         seed=42, local_search=0.0, temperature=0.001,
         random_type=0, min_transfer_distance=0,
         threshold_constraints=0.0, random_trees=0,
         output_prefix=None, print_summary=True,
         *, from_ale=False, from_tool=None,
         target_clade=None, ancestor_map=False, dry_run=False,
         html_report=True, ale_min_support=0.05,
         ale_min_family_size=5, ale_cache_dir=None,
         ale_source="trf", ale_parallel="process",
         artra_transfer_kind="all",
         adapter_min_endpoint_hit_rate=None, adapter_quiet=False,
         mcmc=False, mcmc_iters=1000, mcmc_temperature=0.0,
         mcmc_burn_in=None, mcmc_thin=None,
         incremental=True, checkpoint_path=None, checkpoint_interval=60.0,
         constraints_out=None, output_style="short",
         force=False, local_search_max_iterations=0,
         top_k=50, stop_check=None) -> Result
```

- `constraints_path` may be a **str or a list**; several files are all parsed and merged
  (edge-key aggregation is independent of file order), and the output prefix is taken from the
  **first** file. In the multi-file case stdout is unchanged, **stderr** carries one NOTE, and
  `result.run_metadata["multi_file_prefix_note"]` records it.
- Everything after the 10th positional parameter (from `from_ale` on) is **keyword-only**.
- `mcmc_temperature=0.0` is the **`auto` sentinel** (`max(total weight,1)/100`), not "zero
  temperature"; passing `0.0` equals CLI `--mcmc-temperature auto`. Non-positive explicit
  temperatures are rejected.
- `top_k=50` is the **near-optimal collector capacity** (CLI `--near-optimal-top-k`): how many
  **deduplicated** near-optimal orders the robustness/sensitivity summary may consider. It must
  be an integer ≥ 1 (`ranking/ranker.py:check_near_optimal_top_k`), it only has an effect when
  `local_search > 0`, and the effective value is in `run_metadata["near_optimal_top_k"]`
  (03 §3.6b).
- `stop_check=None` is a **cooperative-cancellation predicate** (`()`→`bool`): it interrupts at
  the **iteration boundary** of the local search and of the MCMC chain and delivers the result
  as of that point. Studio's Cancel button uses exactly this (08 §8.9).
- `species_tree_path` / `constraints_path` likewise **transparently accept** `.gz` / gzip
  streams and tar archives (05 §5.7): `api.rank` first unpacks `.tgz` / `.tar.gz` / `.tar`
  safely into a process temporary directory and treats the members as the inputs
  (`run_metadata["archives_expanded"]` keeps it traceable), while single-file gzip is
  decompressed in place by the reading layer.
- `threshold_constraints` has domain **[0, 1]**; `random_type` ∈ {0,1,2}; `temperature` > 0;
  `local_search` ≥ 0; `random_trees` ≥ 0; `mcmc_iters` > 0; `top_k` ≥ 1. All enforced by
  `Ranker._validate_params`, which raises `ValueError` — **including on the `dry_run=True`
  path**.
- `adapter_min_endpoint_hit_rate=None` means "use the adapter-layer default" (`0.5`).
- With `force=False`, an existing product raises `FileExistsError` (the CLI turns it into exit
  code 3).

### Mapping to CLI options

| Python parameter | CLI option |
|------------------|------------|
| `species_tree_path` | positional `SPECIES_TREE` |
| `constraints_path` (str or list) | positional `CONSTRAINTS` |
| `seed` | `--seed` |
| `local_search` | `--local-search/--ls` |
| `temperature` | `--temperature/--t` |
| `random_type` | `--random-type/--r` |
| `min_transfer_distance` | `--min-transfer-distance/--d` |
| `threshold_constraints` | `--threshold-constraints/--ts` |
| `random_trees` | `--random-trees/--rd` |
| `from_tool` | `--from` |
| `from_ale` | `--from-ale` |
| `ale_min_support` | `--ale-min-support` |
| `ale_min_family_size` | `--ale-min-family-size` |
| `ale_cache_dir` | `--ale-cache-dir` |
| `ale_source` (default `"trf"`) | `--ale-source` |
| `ale_parallel` | `--ale-parallel` |
| `artra_transfer_kind` | `--artra-transfer-kind` |
| `adapter_min_endpoint_hit_rate` | `--min-endpoint-hit-rate` |
| `adapter_quiet` | `--quiet-adapters` |
| `target_clade` / `ancestor_map` | `--target-clade` / `--target-clade-ancestor-map` |
| `dry_run` | `--dry-run` |
| `html_report` (default True) | `--no-html` (inverted) |
| `mcmc` / `mcmc_iters` / `mcmc_temperature` | `--mcmc` / `--mcmc-iters` / `--mcmc-temperature` |
| `mcmc_burn_in` / `mcmc_thin` | `--mcmc-burn-in` / `--mcmc-thin` |
| `incremental` / `checkpoint_path` / `checkpoint_interval` | `--incremental` / `--checkpoint` / `--checkpoint-interval` |
| `local_search_max_iterations` | `--local-search-max-iters` |
| `top_k` | `--near-optimal-top-k` |
| `constraints_out` | `-o/--constraints-out` |
| `output_style` | `--output-style` |
| `force` | `-f/--force` |
| `output_prefix` | (API-only, equivalent to `-p`) output prefix, defaults to the first constraint file path |
| `print_summary` | (API-only) whether to print the stdout summary |
| `stop_check` | (API/GUI-only, **no CLI option**) cooperative-cancellation predicate, see 4.4.6 |

> CLI-only: `--from-auto` (equivalent to `from_tool="auto"`), `--version`, `--help`.
> Every CLI option reaches `api.rank`; there are no dead options. The converse does not hold:
> `stop_check` exists only on the Python / Studio side — the CLI has no interruption switch
> beyond `Ctrl-C`.

## 4.2 The `Result` object

`rank()` returns `maxtic_next.ranking.ranker.Result` (a dataclass). Key fields:

| Field | Type | Meaning |
|-------|------|---------|
| `constraint_file` | str | The file this run's prefix came from (first constraint file) |
| `input_order` | list[str] | Internal-node ranking implied by the input tree |
| `greedy_order` | list[str] | Greedy heuristic ranking |
| `mixing_order` | list[str] | Mixing heuristic ranking |
| `best_order` | list[str] | Final best ranking |
| `best_source` | str | `"greedy heuristic"` / `"mixing heuristic"`; becomes `"… + local search"` when `--ls > 0` improves the solution |
| `values` | dict | `{"input":..,"greedy":..,"mixing":..}`; with local search also `"local_search"` and `"best"` |
| `uninformative` | dict | `{"total","to_desc","to_leaf","to_anc","to_itself"}` |
| `uninformative_percent` | float/None | **= total ÷ `total_weight` × 100**, not the original's double-counting denominator |
| `from_leaf` | float | Weight of constraints whose donor is a leaf (unsatisfiable, ignored in scoring) |
| `trivial_conflict` | float | Trivially conflicting weight |
| `total_weight` | float | Total constraint weight (computed **after** `--ts` removal when that option is active) |
| `removed_by_threshold_weight` / `removed_by_threshold_count` | float / int | Weight and count removed by `--ts` |
| `informative_count` | int | Number of informative constraints |
| `internal_node_count` | int | Number of internal nodes |
| `ranked_newick` | str | Ranked-tree Newick |
| `similarity_to_input` | float | Kendall similarity with the input tree |
| `conflict_with_input` / `partial_total` | float | Numerator/denominator of the "agree with best tree yet conflict with input tree" ratio in the summary |
| `informative_lines` / `conflicting_lines` / `partial_lines` | list[str] | Raw lines of the three files |
| `informative_file` / `conflicting_file` / `partial_order_file` | str | Paths of the three output files |
| `html_report_file` | str | HTML report path (when produced) |
| `sensitivity_summary` | dict/None | Robustness/sensitivity summary (when local search is on) |
| `mcmc_samples` | list/None | MCMC linear-extension samples (`--mcmc`) |
| `random_stats` | dict/None | Random-tree p-value statistics (`--random-trees>0`), including `value_pvalue`, `sim_pvalue`, `correction="(k+1)/(n+1)"`, `tested_order="delivered best order (after local search)"`, `distribution_file` (the 4th product's path) |
| `warnings` | list[str] | Runtime convention statements and warnings (same source as the trailing stdout lines) |
| `run_metadata` | dict | Run metadata: `params`, `seed`, `temperature`, `time_for_search`, `min_transfer_distance`, `threshold_constraints`, `random_type`, `distance_filter`, `informative_count`, `removed_by_threshold_*`, `adapter_diagnostics`, `constraint_files_merged`, `mcmc*` (incl. `mcmc_status`, `mcmc_temperature_is_auto`, `mcmc_chain_diagnostics`), `local_search_*`, `near_optimal_top_k` (the effective `--near-optimal-top-k`), `cancelled` / `stop_check_requested` / `mcmc_stopped_early` (cooperative cancellation), `archives_expanded` / `archive_notes` (which members an archive input expanded to), and — on the relevant paths — `dry_run`/`dry_run_ok`/`dry_run_n_errors`/`dry_run_n_warnings`, `constraints_out`, `multi_file_prefix_note` |
| `dry_run_report` | str | Preflight report (`dry_run=True`) |

> **There is no per-node "MTC score"** in MaxTiC: the objective is the total violated weight of
> the whole ranking (`ranking/value.py`). Any per-node figure is a property of the input
> constraint set and does not explain the ranking.

## 4.3 `build_constraints(path)` — parsing only

```python
from maxtic_next import build_constraints
cset = build_constraints("constraints.tsv")   # returns a ConstraintSet
for c in cset.constraints:
    print(c.donor, c.receptor, c.weight, c.metadata)
```

Distance filtering does **not** happen here (it is applied by
`ConstraintSet.filter_by_distance` during ranking). Illegal weights (`nan`/`inf`/negative)
raise `ValueError` at the parsing layer, with file name and line number.

## 4.4 Typical usage

### 4.4.1 Basic ranking

```python
from maxtic_next import rank

r = rank("examples/minitree.tree", "examples/Cyano_CUTConstraints.tsv",
         seed=42, print_summary=False, html_report=False)
print("best heuristic:", r.best_source)           # greedy heuristic
print("similarity:", r.similarity_to_input)       # 1.0
print("three files:", r.informative_file, r.conflicting_file, r.partial_order_file)
```

### 4.4.2 Ranking straight from upstream tool output

```python
# RANGER-DTLx reconciliation report
r = rank("species.tree", ["FAM1.dtl", "FAM2.dtl"], from_tool="ranger",
         ale_min_family_size=5, print_summary=False, html_report=False)

# auto-detection
r = rank("species.tree", ["out.recphyloxml"], from_tool="auto", html_report=False)

# adapter conventions and drop counters
print(r.run_metadata["adapter_diagnostics"])
```

ALE defaults to `ale_source="trf"` (transfer events, the convention of ALE's own MaxTiC
integration). To line up with reconciliation-event results pass `ale_source="rec"` explicitly — these are
**different constraint sets**, not two spellings of one thing.

### 4.4.3 Two stages: emit unified constraints first, then rank

```python
# stage 1: convert ALE .uml_rec files into a unified constraint file (no ranking)
rank("species.tree", ["fam1.uml_rec", "fam2.uml_rec"], from_tool="ale",
     ale_cache_dir="cache/", constraints_out="constraints.tsv")

# stage 2: rank the constraint file
rank("species.tree", "constraints.tsv", local_search=180)
```

### 4.4.4 Local search + robustness summary

```python
r = rank("species.tree", "constraints.tsv", local_search=30,
         print_summary=False)
summ = r.sensitivity_summary   # dict with node_position_distribution, etc.
print(summ["title"])           # robustness/sensitivity summary over solutions visited by the local search
print("near-optimal orders:", summ["n_orders_collected"])
print("capacity K:", summ["top_k"], r.run_metadata["near_optimal_top_k"])
print("distinct orders the search actually visited:", summ["n_unique_orders_total"])
```

To reproduce this path across machines, add `local_search_max_iterations=2000000`.
The summary only counts the `top_k` best **deduplicated** near-optimal orders (default 50), so
reporting those frequencies **requires** reporting K as well; raising it (`top_k=400`) merely
widens the support set — it does **not** turn the summary into a posterior (03 §3.6b,
08 §8.1b).

### 4.4.5 Preflight only, no ranking

```python
r = rank("species.tree", "constraints.tsv", dry_run=True)
ok = r.run_metadata["dry_run_ok"]        # bool — the sole basis of the CLI exit code
print(r.dry_run_report)                  # human-readable report
print(r.run_metadata["dry_run_n_errors"])
```

This holds in adapter mode too: when the upstream output **cannot be parsed**, `dry_run_ok` is
`False` and `dry_run_n_errors >= 1`, so the CLI exits 1 rather than reporting a pass:

```python
r = rank("examples/adapters/ale/species.tree", "truncated.recphyloxml",
         from_tool="eccetera", dry_run=True, print_summary=False)
r.run_metadata["dry_run_ok"]          # False
r.run_metadata["dry_run_n_errors"]    # 1
```

### 4.4.6 Cooperative cancellation (`stop_check`)

```python
import itertools
from maxtic_next import rank

calls = itertools.count()
r = rank("species.tree", "constraints.tsv", local_search=30.0, temperature=60.0,
         stop_check=lambda: next(calls) > 50)      # turns true after the 50th check

r.run_metadata["cancelled"]                     # True
r.run_metadata["local_search_iterations"]       # 51 — planned 30 s, stopped at an iteration boundary
r.values["best"]                                # best value as of that point, delivered normally
```

Guarantees: it stops only at an **iteration boundary**, returns the best order so far, writes
the products as usual, and raises nothing; with no `stop_check` the behaviour is byte-for-byte
what it always was. `MCMCSampler.sample(stop_check=...)` accepts it too; the **one exception**
is an MCMC chain cancelled before it produced any sample — that raises `RuntimeError` rather
than returning statistics over zero samples. The underlying switches are
`Ranker.run(stop_check=...)` → `optimisation_locale(stop_check=...)`, the latter writing
`stats_out["cancelled"]` back into `run_metadata`. See 08 §8.9.

### 4.4.7 Feeding compressed input directly

```python
rank("species.tree.gz", "constraints.tsv.gz")                      # the reading layer decompresses
rank("species.tree", "examples/reconciliations.tgz", from_tool="ale")
#  ↑ the archive is unpacked safely into a process temp dir; each member is its own input
r.run_metadata["archives_expanded"]    # {archive: [member paths…]}, for traceability
```

When a multi-member archive is read as **one** text input, `CompressedArchiveError` is raised,
and its message contains the runnable `tar -xzf` / `tar -xzOf` commands. Scope and measurements
are in 05 §5.7.

## 4.5 Low-level adapter API (advanced)

```python
from maxtic_next.tree.tree import Tree
from maxtic_next.constraints.adapters import registry

species = Tree(); species.read_newick(open("species.tree").read())

# explicit tool
cset = registry.convert("ranger", species, ["FAM1.dtl"], min_family_size=0)
# auto-detection
cset = registry.convert_auto(species, ["out.recphyloxml"])
# available tool names
print(registry.REGISTRY.names())   # ['ale','alerax','artra','eccetera','ranger']
```

The `ConstraintSet` returned by `convert` / `convert_auto` carries `.diagnostics` (endpoint hit
rate, per-reason drop counters, `donor_endpoint_convention`, `weight_semantics`). `--from auto`
handles **mixed formats**: files are grouped per tool, each parsed by its own
adapter and combined with `registry.merge_constraint_sets(...)`, while the fact that "donor
endpoint levels differ / weight scales are not commensurable" is reported on stderr and recorded
in `diagnostics`.

Per-tool convenience functions:

```python
from maxtic_next.constraints.adapters import (
    convert_from_ale, convert_from_ranger_dtl, convert_from_eccetera,
    convert_from_artra, convert_from_alerax,
)
```

> Prev: [03 · CLI reference](03_cli_reference.md) ｜ Next: [05 · Input/output formats](05_io_formats.md)
