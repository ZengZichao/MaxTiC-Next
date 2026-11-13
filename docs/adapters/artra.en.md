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
m3 = LCA[H117,H14]: Replacing Transfer, Mapping --> H15, Recipient --> H125
m19 = LCA[H159,H17]: Replacing Transfer, Mapping --> H17, Recipient --> H159
m93 = LCA[H10,H90]: Additive Transfer, Mapping --> H90, Recipient --> H10
```

> **The block header is identical to RANGER-DTLx's** (it is not ARTra-specific, yet it is
> printed by ARTra's **own Java code**):
> `参考软件-ARTra-代码/JavaCode/dtl_algo/DTLReconciliation.java:619` (the stdout variant
> `printOutputOfDTRLConsole`) and `:1387` (the `PrintWriter` variant `printOutputOfDTRL`)
> print exactly
> ` ------------ Reconciliation for Gene Tree 1 (rooted) -------------`, **identical in
> shape** to what upstream `RANGER-DTLx-代码/Ranger-DTLx/DTL-algorithm.h:2170` outputs.
> Chunking matches only the `Reconciliation for Gene Tree` **prefix**, so both the full
> upstream header and the abbreviated `Reconciliation for Gene Tree 1:` used by the
> RANGER-DTLx sample are recognized (ARTra's own sample
> `examples/adapters/artra/FAM1.txt` carries the full `(rooted)` header). The event names
> are rewritten into `Replacing Transfer` / `Additive Transfer` by the transfer classifier
> selected with `-h`: `heuristic` defaults to `"ml"` at `DTLReconciliation.java:501`, so the
> labels in `ARTra-程序/output.txt` come from the machine-learning classifier
> `modifySolutionByML` (`:764` / `:767`), while the rule-based heuristics rewrite them at
> `:608` / `:610` (lost-gene), `:1465` / `:1467` (gene-frequency) and `:1509-1518`
> (mapping-count). The event lines are printed by `DTLReconciliation.java:637-641` (stdout)
> and `:1405-1409` (the `PrintWriter` path that writes `output.txt`) — ARTra **never writes**
> `Edge, Parent =` (that fragment
> only occurs in RANGER-DTLx's C++ printers). Inside `LCA[.., ..]` ARTra writes no space
> (`Node.java:466`) while RANGER-DTLx does; parsing only reads the labels after
> `Mapping -->` / `Recipient -->`, so both spellings work.

## 2. Field Extraction and Semantic Mapping

- Consistent with RANGER-DTLx: `Mapping --> D` = donor, `Recipient --> R` = receptor,
  producing constraint `(D, R)`;
- `Replacing Transfer` / `Additive Transfer` / plain `Transfer` are all recognized (the
  trigger condition is that the same line contains **both** `Mapping` and `Recipient`, see
  `parse_transfer_line`);
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
`support = count / denominator` (denominator = the number of samples, i.e. the tool's
declared value; when none is declared, the block count may serve as the denominator only
if there is **more than one block** and all blocks share an identical leaf set — a
single-block report never qualifies; ecceTERA does not check the leaf set, see
[eccetera.en.md](eccetera.en.md) §2)
→ keep `> min_support`;
`distance = distance_from(receptor, donor)` (**self-loop exception**: when
`donor == receptor`, `distance = None` is written and the constraint is not removed by
`--min-transfer-distance`); invalid endpoints / families that are too small are dropped
and **counted and reported by reason**. ⚠️ ARTra prints
`Total number of optimal solutions: 24` yet outputs only **one** optimal solution, so
using the block count as the denominator necessarily degenerates for this tool; when the
denominator cannot be established the weight is an **integer count**, and `min_support`
has no [0,1] semantics. The donor endpoint likewise takes `donor_itself` (one species
level away from ALE). Observability of drops matches the other adapters: every dropped
transfer is **counted by reason** (label misses in `dropped_by_label_miss`; numeric-ID
resolution runs only on the ecceTERA path, so `numeric_id_resolutions` is always 0 here),
a diagnostic summary is printed to stderr when the conversion finishes and the same dictionary
is attached to `ConstraintSet.diagnostics`; an endpoint hit rate of **0** raises
`LabelMismatchError`,
and one below `min_endpoint_hit_rate` (default `0.5`) triggers a warning
(`--quiet-adapters` only silences printing; counting and error raising are unaffected).

## 5. Known Limitations

- Endpoints are species labels and must match the species tree labels;
- The machine-learning classification confidence is not folded into the weight (the
  weight remains the sampling support); if probability-weighted constraints are needed,
  generate weighted constraints in an upstream preprocessing step and take the
  text-constraint path.
- This section is a summary: the **complete semantics** of the endpoint levels, weight
  conventions and the diagnostics counters are governed by the normative semantic notes
  [README_semantics.md](../../src/maxtic_next/constraints/adapters/README_semantics.md) /
  [README_semantics.en.md](../../src/maxtic_next/constraints/adapters/README_semantics.en.md).

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
