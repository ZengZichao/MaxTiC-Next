# ARTra Adapter Specification

> 🌐 中文版：[artra.md](artra.md)

- **Tool**: ARTra (Additive and Replacing Transfer classifier), from `ARTra-代码` /
  `ARTra-程序`. On top of RANGER-DTL-style reconciliations, it classifies each transfer as
  an **Additive** or **Replacing** transfer using rule-based heuristics + machine learning.
- **Adapter**: `maxtic_next.constraints.adapters.artra` (`ARTraAdapter` /
  `convert_from_artra`), parsing core `_recon_report` (shared with RANGER-DTLx).
- **Tool name (registry / CLI)**: `artra` (`--from artra`).

## 1. Real Output Format (verified from samples)

ARTra's final output (as evidenced by `ARTra-程序/output.txt`) is **identical in shape** to
the RANGER-DTLx reconciliation report, except that event names carry an extra
`Replacing` / `Additive` prefix:

```
 ------------ Reconciliation for Gene Tree 1 (rooted) -------------
Species Tree:
(...);
Gene Tree:
(...);
Reconciliation:
H117_0: Leaf Node
 = LCA[H117, H14]: Replacing Transfer, Mapping --> H15, Recipient --> H125
 = LCA[H159, H17]: Additive Transfer, Mapping --> H17, Recipient --> H159
```

## 2. Field Extraction and Semantic Mapping

- Consistent with RANGER: `Mapping --> D` = donor, `Recipient --> R` = receptor,
  producing constraint `(D, R)`;
- `Replacing Transfer` / `Additive Transfer` / plain `Transfer` are all recognized (the
  trigger condition is the presence of `Recipient -->`);
- Section headers such as `Species Tree:` / `Gene Tree:` and newick lines contain no
  `Recipient -->` and are ignored automatically.

## 3. Transfer-Class Filtering (ARTra-specific)

MaxTiC's temporal constraint semantics are indifferent to whether a transfer is additive
or replacing, so **all transfers are counted by default** (`transfer_kind="all"`). If a
study cares about only one class:

- `transfer_kind="replacing"`: count only replacing transfers;
- `transfer_kind="additive"`: count only additive transfers.

(Class filtering screens lines in the main process before parsing, without using
parallelism / cache, to avoid polluting the shared cache keys.)

## 4. Support / Weight / Distance / Error Handling

Same as RANGER-DTLx: per-block counting → reciprocal subtraction →
`support = count / declared sample count` → keep `> min_support`;
`distance = distance_from(receptor, donor)` (**self-loop exception**: when
`donor == receptor`, `distance = None` is written and the constraint is not removed by
`--min-transfer-distance`); invalid endpoints / families that are too small are dropped
and **counted and reported by reason**. ⚠️ ARTra prints
`Total number of optimal solutions: 24` yet outputs only **one** optimal solution, so
using the block count as the denominator necessarily degenerates for this tool; when the
denominator cannot be established the weight is an **integer count**, and `min_support`
has no [0,1] semantics. The donor endpoint likewise takes `donor_itself` (one generation
off from ALE).

## 5. Known Limitations

- Endpoints are species labels and must match the species tree labels;
- The machine-learning classification confidence is not folded into the weight (the
  weight remains the sampling support); if probability-weighted constraints are needed,
  generate weighted constraints in an upstream preprocessing step and take the
  text-constraint path.

## 6. Usage

```bash
maxtic-next species.tree FAM1.txt --from artra --artra-transfer-kind all
maxtic-next species.tree FAM1.txt --from artra --artra-transfer-kind replacing
```

```python
from maxtic_next.constraints.adapters import convert_from_artra
cset = convert_from_artra(species_tree, ["FAM1.txt"], transfer_kind="all", min_family_size=5)
```

## 7. Samples

`tests/data/artra_example/FAM1.txt` + `species.tree`. See the regression test
`tests/test_artra_adapter.py`.
