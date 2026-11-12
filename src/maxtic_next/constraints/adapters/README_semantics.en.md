# Adapter Semantics (Donor Endpoint Hierarchy / Weight Convention / Observability)

> 🌐 中文版：[README_semantics.md](README_semantics.md)

This file is the **normative semantic annotation** for the adapters in
`src/maxtic_next/constraints/adapters/` (the ALE / AleRax / ARTra / ecceTERA /
RANGER-DTLx adapters and the shared registry code).
It answers three questions: which species-tree level each tool's donor endpoint
actually refers to, whether the weight is really a [0,1] support value, and how the
user can see it when "constraints become fewer". Wherever the CLI / manual
(`docs/adapters/*.md`) does not repeat this content, this file prevails.

## 1. Donor Endpoints Differ by One Species-Tree Level

A transfer event `donor -> receptor` can be written as two different kinds of MaxTiC
constraints, and the two are **not equivalent**:

| Route | What the donor endpoint refers to | Diagnostic value `donor_endpoint_convention` | Upstream evidence |
|---|---|---|---|
| ALE (`.uml_rec`) | **The parent of the donor**: `parent(donnor)` (trf) or the lineage reached by `donnor_search()` (rec) | `parent_of_donor` | original `constraints_from_reconciliations.py:160-161`; `ale.py` replicates this faithfully |
| RANGER-DTLx / ARTra | **The donor itself**: `Mapping --> <donor>` | `donor_itself` | `DTL-algorithm.h:1840` |
| ecceTERA | **The donor itself**: `<branchingOut speciesLocation="D">` of the parent clade | `donor_itself` | `DTLGraph.cpp:2738-2869` |
| AleRax | **The donor itself**: first column of `*_transfers.txt` | `donor_itself` | `scripts/extract_families_transfer.py` |

The original MaxTiC ALE integration additionally expresses "the child of the receptor"
as the lower bound (ALE `maxtic/README.md:131`:
"the **father of the donor** branch should be older than the **child of the
receptor** branch"), whereas `Mapping -->` / `branchingOut` give the lineage in which
the event itself resides. **This package does not convert between levels**: each
adapter keeps its tool's native convention (converting it would break equivalence with
that tool's official output) and only records the convention **in the diagnostics**.

When merging the two conventions, `registry.merge_constraint_sets(sets)` and
`registry.convert_auto(..., allow_mixed=True)`:

* print a warning to stderr stating that the constraint sets being merged come from
  **different donor endpoint level conventions** (translation of the warning the program
  emits in Chinese: "待合并的约束集来自**不同的供体端点层级约定**");
* set `ConstraintSet.diagnostics["convention_conflict"] = ["donor_itself",
  "parent_of_donor"]` (likewise, `semantics_conflict` is recorded when the weight
  conventions are incommensurable).

As a result, mixing ALE and RANGER-DTLx-family output under `--from-auto` no longer
silently adds endpoints taken from two different species levels into the same edge —
it says so.

## 2. Weight Convention: the Block Count Is Not the Sampling Denominator

`diagnostics["weight_semantics"]` takes one of the following four values (or
`mixed:a+b` when sets with several conventions are merged), and is returned together
with each constraint's `metadata["weight_semantics"]`:

| Convention | Meaning | When used |
|---|---|---|
| `support_over_declared_samples` | `count / N`, where N is the sample count **declared by the tool itself** | ALE's `<N> reconciled`; RANGER-DTLx/ARTra's `Total number of optimal solutions: N` (`DTL-algorithm.h:1785/:2224`, `output.txt:206`) |
| `support_over_replicate_blocks` | `count / number of blocks` | **ALE**: when there is no `<N reconciled>` header line, the number of reconciliations actually parsed in this file is used as the denominator (`ale.py:640-651`, see `docs/adapters/ale.en.md` §3); **RANGER-DTLx / ARTra**: no declared value, and more than one block whose leaf sets are exactly identical (genuine replicate reconciliations of the same family); a warning is attached. **This leaf-set check exists only in the RANGER-DTLx / ARTra shared core (`_recon_report._leaf_sets_consistent`)**; ecceTERA does not check leaf sets and, when there is no declared value, uses the block count as the denominator whenever the number of blocks is greater than 0 (see `docs/adapters/eccetera.en.md` §2) |
| `integer_count` | **Integer counts**, not [0,1] support values | No declared value and the block count is not trustworthy (blocks with different leaf sets = multiple distinct gene trees, or only a single block), or neither a sample count nor a block count can be obtained at all |
| `posterior_frequency` | Posterior frequencies pre-aggregated upstream | both versions of AleRax's summary files, `*_transfers.txt` / `*_meanTransfers.txt` (`alerax._TRANSFER_SUFFIXES`; no longer divided by any denominator) |

Key points:

1. **Real ARTra output** prints only 1 optimal solution yet declares
   `Total number of optimal solutions: 24` ⇒ the denominator is 24. The old
   implementation normalized by "number of blocks = 1", so all weights became `1.0`
   and thresholds such as `--ale-min-support 0.95` became no-ops.
2. Under the `integer_count` convention, `min_support` is still compared against the
   **same numeric values** (it does not pretend to pass), but the `[0,1]` semantics no
   longer hold: `diagnostics["min_support_is_fractional"] = False`, with an explicit
   warning on stderr.
3. Weights under different conventions **must not be added or compared**: ALE's 0.3 and
   ARTra's count of 3 are not the same thing.

## 3. Drops Must Be Visible

The `ConstraintSet.diagnostics` returned by every adapter's `convert()` contains at
least:

```
tool, donor_endpoint_convention, weight_semantics, sample_denominator,
files_seen, blocks_seen, trees_seen, transfers_seen, transfers_resolved,
transfers_unpaired, numeric_id_resolutions, constraints_kept, edges_kept,
dropped_by_label_miss, dropped_by_support, dropped_by_family_size,
dropped_by_family_size_probe_failed, malformed_lines, families_skipped,
self_loop_transfers, family_size_probe_ok, endpoint_hit_rate,
unresolved_endpoint_samples, species_label_samples, warnings[], per_file[]
```

and prints a one-line stderr summary. The line below is the program's **literal**
output: `_diagnostics.summary_line()` always uses the Chinese labels shown, regardless
of the locale, so it is reproduced verbatim rather than translated (English gloss:
"[artra] donor endpoint convention=donor_itself, weight convention=
support_over_declared_samples, denominator=24, endpoint hit rate=100.0% | files=1
blocks=1 transfers=13 kept=0 dropped[label mismatch=0, support=13, family size=0,
family-size probe failed=0, self-loop=0]"):

```
[artra] donor 端点约定=donor_itself 权重口径=support_over_declared_samples 分母=24
端点命中率=100.0% | 文件=1 块=1 转移=13 保留=0 丢弃[标签不匹配=0, 支持度=13,
家族规模=0, 家族规模探测失败=0, 自环=0]
```

Three hard rules:

* **Family-size probe failure ⇒ skip filtering**: when `Leaf Node` lines
  (RANGER-DTLx/ARTra) or `<leaf>` elements (ecceTERA) cannot be probed, `family_size` would
  be 0, and the old implementation's
  `if family_size <= min_family_size: transfers = {}` would then **wipe out all
  transfers** (`--ale-min-family-size 0` triggered this invariably). It now skips the
  filter and warns instead.
* **Endpoint hit-rate guard**: adapters such as ecceTERA write the numeric
  `node->getId()` for internal species nodes (`DTLGraph.cpp:2818-2821`; the numbering
  scheme is the bottom-up, level-by-level breadth-first traversal in
  `MySpeciesTree.cpp:112-136`). `_species.numeric_species_id_map` resolves them back;
  if after resolution **none** of them match ⇒ raise `LabelMismatchError` (the error
  message also includes samples of both the unmatched specimen labels and the
  species-tree labels); if a partial match falls below the threshold (default 50%,
  configurable via `min_endpoint_hit_rate`) ⇒ warn. It never silently returns `{}`.
  **ALE is the exception**: `convert_from_ale` has no `min_endpoint_hit_rate`
  parameter and calls `finalize_adapter_diagnostics(..., raise_on_zero_hit=False)`
  with a 0.0 threshold, so even a 0% hit rate does not raise — it only prints an
  error-level notice; passing the parameter through the registry / CLI is ignored and
  recorded in `adapter_params_not_applicable` (see `docs/adapters/ale.en.md` §3, §4).
* **XML parse failure ⇒ raise**: `UpstreamParseError` carries the file, the
  line/column numbers and the original exception; the pure function
  `parse_recphyloxml(..., strict=False)`, kept with its old signature for backward
  compatibility, at least prints a loud warning containing the file and the
  line/column, and additionally supports `strict=True`. Corrupt files can henceforth
  be distinguished from "the tool detected no transfers".

## 4. Distances and Self-Loops

* `_species.distance_from` returns `None` for labels **not present in the species
  tree** (the old implementation returned `0.0`, which was indistinguishable from the
  true distance of 0 for "donor = receptor", and that fake 0 would let
  `--min-transfer-distance` treat the constraint as having a valid distance and delete
  it silently). In `ConstraintSet.filter_by_distance`, `None` is synonymous with "the
  text input has no 4th column": the constraint is **kept**.
* `filter_by_distance(d)` now returns and records
  `diagnostics["distance_filter"] = {threshold, before, kept, dropped,
  without_distance_column, self_loops_dropped, ignored}`, where `ignored=True` means
  "threshold > 0 but not a single constraint was removed". Callers
  (`ranking/ranker.py`, `dry_run.py`) must report at **error level** that "this option
  was ignored for the current input" when **no constraint carries a distance column**
  (if the input does have distance columns and simply none of them falls below the
  threshold, the message is warning level, "not a single constraint was removed"),
  rather than claiming the opposite, that "all constraints will be dropped".
* Self-loops (donor == receptor, "donor = receptor lineage") are counted in
  `self_loop_transfers` and no longer vanish from the statistics; the text path is
  handled by the original `ranker`'s "to itself" branch.

## 5. Cache Keys Include a Topology Fingerprint

`_species.species_cache_identity(tree)` = `v{CACHE_KEY_VERSION}|` + a digest of the
label set + a **topology fingerprint** (`sha256` over the sorted
`parent_label->child_label` pairs). The file-level pickle cache keys of ALE /
RANGER-DTLx / ARTra / ecceTERA all use it, with the version prefix: entries from older versions that
contain only the label digest will **never be read**.

Consequently, when `--ale-cache-dir` reuses a directory across two trees with the same
label set but different parent-child relations, it will no longer hit donors from the
old topology (ground truth without cache: `E`; stale cache hit: `F`).

## 6. ALE Default Convention = `trf`

The `source` of `ALEAdapter` / `convert_from_ale` defaults to **`"trf"`**: the MaxTiC
integration in the official ALE repository itself
(the `maxtic/README.md:131-169` inside the ALE source package) specifies that the input
`constraints_from_transfers` is produced from **transfer events** (donor = the parent
of the donor). Measurements confirm that this path matches the official output
constraint-by-constraint (identical edge keys and identical weights), whereas
`"rec"` yields only 13,995 constraints with a total weight of 5482.98 vs 6438.6 —
the default changes the entire constraint set. `"rec"` remains fully available as an
option. When `convert()` finishes, it prints one line to stderr stating which
convention produced the output, along with the edge count and the total weight.

> The upper-layer wiring is already in sync: the `ale_source` parameter of `api.py`
> defaults to `"trf"`, and `--ale-source` in `cli.py` uses `default="trf"`
> (`choices=["rec", "trf"]`), so the CLI/API no longer override this package's
> default.
