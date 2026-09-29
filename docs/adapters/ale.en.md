# ALE Adapter Specification

> 🌐 中文版：[ale.md](ale.md)

- **Tool**: ALE_undated (AmalgamatedLikelihoodEstimation), an independent upstream project.
- **Adapter**: `maxtic_next.constraints.adapters.ale` (`ALEAdapter` / `convert_from_ale`).
- **Tool name (registry / CLI)**: `ale` (`--from ale`; `--from-ale` also accepted).
- **Provenance**: a faithful port of the control flow of the
  `constraints_from_reconciliations.py` shipped with the original MaxTiC (Python 2), including the fact
  that **rec constraints are produced only by nodes carrying a `T@` annotation** (in the
  original, the "2/ calculer la contrainte selon les reconciliations" block sits inside the
  `if e[:2]=="T@"` gate and executes once per T@ event).

## 1. Real Output Format

ALE_undated produces a sampled reconciliation file `<family>.ale.uml_rec` for each gene
family. Structure:

```
<N> reconciled G-s:
<blank line>
(gene_tree_1 with NHX-style event annotations);
(gene_tree_2 ...);
...(N sampled reconciliation trees in total)
```

Internal node event annotations carry a **leading dot** (ALE_undated always emits `.` +
species string); parsing splits events on `.`, e.g. `.59.T@61->62`, `.65`. Event types:

- `T@donor->receptor`: one transfer (donor = donor species, receptor = receptor species);
- `D@species`: duplication;
- Others: speciation (mapped species).

`N` is the number of sampled reconciliation trees and serves as the **support
denominator** (support = count / N).

## 2. Field Extraction and Semantic Mapping

- **Transfer events (trf, default)**: `T@donor->receptor` → constraint
  `(donor, receptor)`; the donor is taken as the parent of the node in the species tree
  (`donnor_search` walks up the parent chain to the first non-T@/D@ event). This is
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
  **both size and weight magnitude** (at the official example scale, edge counts and total
  weights differ by several times); switching `--ale-source` is equivalent to **replacing
  the entire input constraint set**, and all downstream numbers change accordingly. To
  reproduce older results that used `rec` as the default, pass `--ale-source rec`
  explicitly. The convention actually in effect is printed to stderr and recorded in
  `Result.run_metadata["adapter_diagnostics"]`.
- **Endpoint validity (permissive semantics, faithful to the original)**: the receptor
  must not be an extant species (leaf); the donor must resolve to a parent node in the
  species tree. Note that ALE's validation strength differs from the other four adapters
  (where both endpoints must be in the species tree label set) — this is the result of
  line-by-line fidelity to the original; when mixing outputs from multiple tools, make
  sure the upstream species naming is consistent with the user's species tree.
- **Donor endpoint level (one generation off from the other four adapters)**: this
  adapter records `donor_endpoint_convention = "parent_of_donor"` (taking
  `parent(donor)`), whereas RANGER-DTLx / ARTra / ecceTERA / AleRax record
  `"donor_itself"` (taking the donor of `Mapping -->` itself / the parent node of
  `branchingOut`, etc.). The endpoints of the same transfer event **differ by one species
  level** under the two conventions; therefore when `--from-auto` encounters mixed input,
  each group is resolved with its own adapter and a warning "endpoint levels and weight
  semantics are not commensurable" is issued, instead of adding both directly onto the
  same edge.

## 3. Support / Weight / Distance

- Within a single family, counts are aggregated per edge key and **reciprocal
  subtraction** removes symmetric noise;
- `support = count / N`; only `support > min_support` is kept (default 0.05);
- `distance = distance_from(receptor, donor)`: species tree topological distance (branch
  lengths treated as 1); **self-loop exception**: when `donor == receptor`,
  `distance = None` is written (no distance information) and the constraint is not
  removed by `--min-transfer-distance`. This matches the behavior of RANGER-DTLx /
  ecceTERA / ARTra / AleRax — ALE was the **last** to be aligned; see the regression test
  `tests/test_ale_selfloop.py`, and §3.2.3 of `docs/manual/en/03_cli_reference.md` for
  symptoms and measurements.

## 4. Error Handling and Validation

| Situation | Handling |
|-----------|----------|
| Family size `leaf count <= min_family_size` (default 5) | Skip the family |
| donor / receptor not in the species tree | The transfer is dropped |
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
# The default of source is "trf" (DEFAULT_SOURCE, the convention of ALE's official MaxTiC
# integration); source="rec" is written explicitly below to demonstrate how to reproduce
# the old convention — it is not the default behavior.
cset = convert_from_ale(species_tree, ["fam1.ale.uml_rec"], source="rec",
                        cache_dir="cache/", parallel="process")
print(cset.diagnostics["donor_endpoint_convention"])  # 'parent_of_donor'
```

## 7. Samples

`tests/data/ale_example/rec.uml_rec` + `species.tree`. See the regression tests
`tests/test_ale_adapter.py` / `tests/test_ale_parallel.py`.
