# ALE Adapter Specification

> 🌐 中文版：[ale.md](ale.md)

- **Tool**: ALE_undated (AmalgamatedLikelihoodEstimation), an independent upstream project.
- **Adapter**: `maxtic_next.constraints.adapters.ale` (`ALEAdapter` / `convert_from_ale`).
- **Tool name (registry / CLI)**: `ale` (`--from ale`; `--from-ale` also accepted).
- **Provenance**: a faithful port of the control flow of the
  `constraints_from_reconciliations.py` shipped with the original MaxTiC (Python 2),
  including the fact that **rec constraints are produced only by nodes carrying a `T@`
  annotation** (in the original, the "2/ calculer la contrainte selon les reconciliations"
  block sits inside the `if e[:2]=="T@"` gate and executes once per T@ event).

## 1. Real Output Format

ALE_undated produces a sampled reconciliation file `<family>.ale.uml_rec` for each gene
family. Structure:

```
<N> reconciled G-s:
<blank line>
(gene_tree_1 with NHX-style event annotations);
(gene_tree_2 ...);
...
```

The first line is the **sample-count declaration line**: upstream ALE_undated actually
prints `<N> reconciled G-s:`, while a trimmed sample (e.g.
`tests/data/ale_example/rec.uml_rec`) only carries `<N> reconciled` — the adapter tests
**only whether the second token is exactly `reconciled`** (`ale.py:261-262`), the rest of
the declaration line playing no part, so both forms are recognized. Exactly one line
always separates the declaration line from the first reconciliation tree (not parsed;
blank in upstream output, the family name in this repository's samples). Internal node
event annotations carry a **leading dot** (ALE_undated always emits `.` + species
string); parsing splits events on `.`, e.g. `.59.T@61->62`, `.65`. Event types:

- `T@donor->receptor`: one transfer (donor = donor species, receptor = receptor species);
- `D@species`: duplication;
- Others: speciation (mapped species).

`N` is the number of sampled reconciliation trees and serves as the **support
denominator** (support = count / N).

## 2. Field Extraction and Semantic Mapping

- **Transfer events (trf, default)**: `T@donor->receptor` → constraint
  `(donor, receptor)`; the donor is taken as the **direct parent** of that node in the
  species tree (`parent_map.get(donnor)`, a `label -> parent_label` table precomputed from
  the species tree; `donnor_search`, which walks up the parent chain of the
  reconciliation tree, is used **only** by the rec branch below). This is
  exactly the direct input prescribed by the **official ALE MaxTiC integration**
  (`constraints_from_transfers` in `maxtic/constraints_from_reconciliations.py`), hence
  this package sets the default to `trf`.
- **Reconciliation events (rec)**: applied only to **nodes whose annotation contains a
  `T@` event** (matching the original gate); the receptor is the node's own mapped species
  (or, via child-direction matching plus recursive collection with `receptor_search`; a
  `D@` node recursively merges the receptors of both subtrees), and the donor is the
  result of walking up with `donnor_search`. Pure speciation/duplication nodes produce no
  rec constraints.
- **⚠️ `rec` and `trf` are not equivalent**: the constraint sets they produce differ in
  **both count and weights** (at the official example scale, `rec` yields only 13,999
  constraints with a total weight of 5483.31, vs 6438.6 for `trf`); switching
  `--ale-source` is equivalent to **replacing the entire input constraint set**, and all
  downstream numbers change accordingly. To reproduce older results that used `rec` as
  the default, pass `--ale-source rec` explicitly. The convention actually in effect is
  printed to stderr and recorded in `Result.run_metadata["adapter_diagnostics"]`.
- **Endpoint validity (permissive semantics, faithful to the original)**: the receptor
  must not be an extant species (leaf); the donor must resolve to a parent node in the
  species tree. Note that ALE's validation strength differs from the other four adapters
  (where both endpoints must be in the species tree label set) — this is the result of
  line-by-line fidelity to the original; when mixing outputs from multiple tools, make
  sure the upstream species naming is consistent with the user's species tree.
- **No endpoint hit-rate guard**: `convert_from_ale` has **no** `min_endpoint_hit_rate`
  parameter and finishes with
  `finalize_adapter_diagnostics(..., 0.0, raise_on_zero_hit=False)`, so an endpoint hit
  rate of **0** does **not** raise `LabelMismatchError` — it only prints an error-level
  notice to stderr; the other four adapters raise on a 0 hit rate. A
  `min_endpoint_hit_rate` passed through the registry / CLI is dropped by signature
  filtering and recorded in `diagnostics["adapter_params_not_applicable"]`.
- **Donor endpoint level (one species level away from the other four adapters)**: this
  adapter records `donor_endpoint_convention = "parent_of_donor"` (taking
  `parent(donor)`), whereas RANGER-DTLx / ARTra / ecceTERA / AleRax record
  `"donor_itself"` (taking the donor of `Mapping -->` itself / the
  `branchingOut speciesLocation` of the parent clade, etc.). The endpoints of the same
  transfer event **differ by one species level** under the two conventions; therefore
  when `--from-auto` encounters mixed input, each group is resolved with its own adapter
  and a warning "endpoint levels and weight conventions are not commensurable" is
  issued (the program emits this warning in Chinese; quoted here is its translation),
  instead of adding both directly onto the same edge.

## 3. Support / Weight / Distance

- Within a single family, counts are aggregated per edge key and **reciprocal
  subtraction** removes symmetric noise;
- **Three-level denominator degradation** (`_merge_payload`): (1) the `<N> reconciled`
  sample count declared in the file header
  (`weight_semantics = support_over_declared_samples`); (2) without a declared value, the
  number of reconciliation trees **actually parsed from this file**, `trees_seen`
  (`support_over_replicate_blocks`, with a warning); (3) when neither is available,
  `integer_count`: the weight is emitted as an **integer count**, in which case
  `min_support` **has no [0,1] support semantics** (the threshold degenerates to "appears
  at least once"), `diagnostics["min_support_is_fractional"] = False`, and stderr warns
  about it;
- `support = count / N`; only `support > min_support` is kept (default 0.05);
- `distance = distance_from(receptor, donor)`: species tree topological distance (branch
  lengths treated as 1); **self-loop exception**: when `donor == receptor`,
  `distance = None` is written (no distance information) and the constraint is not
  removed by `--min-transfer-distance`. This matches the behavior of RANGER-DTLx /
  ecceTERA / ARTra / AleRax — ALE was the **last** to be aligned; see the regression test
  `tests/test_ale_selfloop.py`, and [§3.2.3](../manual/en/03_cli_reference.md) for
  symptoms and measurements.

## 4. Error Handling and Validation

| Situation | Handling |
|-----------|----------|
| Family size `leaf count <= min_family_size` (default 5) | Skip the family |
| donor / receptor not in the species tree | The transfer is dropped |
| Endpoint hit rate of 0 | No exception raised (ALE has no hit-rate guard); an error-level notice on stderr + recorded in the diagnostics |
| Parse error on a single file (parallel mode) | Automatically falls back to sequential recomputation of that file |
| Cache hit (`cache_dir`) | Reused directly (resumable runs) |

## 5. Known Limitations

- ALE's donor/receptor are species labels and must match the species tree labels (they
  usually share the same origin, so they do).
- The `.uml_rec` event annotation format is sensitive to whitespace / separators; the
  port is byte-faithful to the original. Non-ALE_undated variant formats require separate
  adaptation.

## 6. Usage

```bash
maxtic-next species.tree fam1.ale.uml_rec fam2.ale.uml_rec --from ale \
    --ale-source rec --ale-min-support 0.05 --ale-min-family-size 5
```

```python
from maxtic_next.constraints.adapters import convert_from_ale
# The default of source is "trf" (DEFAULT_SOURCE, ALE's official MaxTiC integration
# convention); the explicit source="rec" below only demonstrates reproducing the old
# convention, it is not the default behavior.
cset = convert_from_ale(species_tree, ["fam1.ale.uml_rec"], source="rec",
                        cache_dir="cache/", parallel="process")
print(cset.diagnostics["donor_endpoint_convention"])  # 'parent_of_donor'
```

## 7. Samples

- `tests/data/ale_example/rec_real.uml_rec` + `species.tree`: the **real upstream output**
  with leading dots (header `1 reconciled G-s:`, event annotations `.59.T@61->62` / `.65`),
  matching the format described in §1;
- `tests/data/ale_example/rec.uml_rec` + `species.tree`: a trimmed sample (header
  `1 reconciled`, event annotations `59.T@61->62` / `X65`, no leading dots).

Both declaration-line forms are recognized by the adapter (`<N> reconciled G-s:` in
`rec_real.uml_rec` vs `<N> reconciled` in `rec.uml_rec`), and the two samples were
measured to produce constraint multisets over `(donor, receptor, weight)` that are
**identical edge by edge** on the same `species.tree`: one `65 -> 59, w = 1.0` under
`--ale-source rec` and one `65 -> 62, w = 1.0` under `--ale-source trf` in both cases (the
samples also differ in the leading-dot annotation style described in §1, which likewise
changes nothing). The only difference is `metadata["family"]`, which is taken from the file
name (`rec_real` / `rec`). The regression tests assert the expected constraint pair for
each sample **separately**:
`tests/test_ale_adapter.py::test_convert_endpoints_exist_in_species_tree` and
`test_convert_transfer_source_produces_65_62` use `rec.uml_rec`, while
`test_real_format_sample_end_to_end` uses `rec_real.uml_rec`; no test compares the two
samples against each other. `tests/test_ale_parallel.py` compares the sequential and
parallel results of the **same file**.
