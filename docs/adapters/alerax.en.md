# AleRax Adapter Specification

> 🌐 中文版：[alerax.md](alerax.md)

- **Tool**: AleRax (Bayesian gene tree–species tree reconciliation and scoring), from
  `AleRax-代码`.
- **Adapter**: `maxtic_next.constraints.adapters.alerax` (`AleRaxAdapter` /
  `convert_from_alerax`).
- **Tool name (registry / CLI)**: `alerax` (`--from alerax`).

## 1. Real Output Format (verified against scripts)

Within its output directory, AleRax generates a transfer frequency file for each gene
family: older AleRax releases wrote `reconciliations/summaries/<family>_transfers.txt` (verified
via `AleRax-代码/scripts/extract_families_transfer.py`); in current versions, the
summaries-writing code in `AleOptimizer.cpp` generates `<family>_meanTransfers.txt`. Both
naming schemes are accepted by the adapter (case-insensitive). Each line is:

```
<donor_species> <receptor_species> <frequency>
```

Here `frequency` is the **expected frequency** of the transfer across posterior samples
(a float). Unlike ALE / RANGER's "count / number of samples", AleRax has already
**pre-aggregated within the family into frequencies**, so the adapter **uses the frequency
directly as the weight** (no longer dividing by a block count).
`weight_semantics = posterior_frequency`.

The donor endpoint is taken from the donor column of the table itself
(`donor_endpoint_convention = donor_itself`), which is **one species level away** from
ALE's `parent_of_donor`; when `--from-auto` mixes ALE and AleRax inputs, the groups are
resolved separately and a warning "endpoint levels and weight semantics are not
commensurable" is issued, instead of adding both directly onto the same edge.

## 2. Input Forms and Automatic Expansion

`inputs` accepts three forms, all uniformly expanded into a list of transfer frequency
files by `expand_alerax_inputs`:

1. The AleRax output **root directory** (`reconciliations/summaries` is appended
   automatically, then globbed);
2. The `reconciliations/summaries` **directory** (globbed directly);
3. Transfer frequency **files** given directly.

Auto-detection: a directory containing `reconciliations/summaries` or
`*_transfers.txt` → `alerax`.

## 3. Field Extraction and Semantic Mapping

- Each line `donor receptor frequency` → constraint `(donor, receptor)`, weight =
  frequency;
- Endpoints must be in the species tree label set, otherwise the transfer is dropped;
- Family name = file name minus the transfer-frequency file suffix;
- Across families, frequencies are **summed** per edge key (reciprocal subtraction is
  applied within each family first).

## 4. Support / Weight / Distance

- `weight = support = frequency` (no block-count denominator); only
  `frequency > min_support` is kept (default 0.05);
- `distance = distance_from(receptor, donor)`; **self-loop exception**: when
  `donor == receptor`, `distance = None` is written (no distance information) and the
  constraint is not removed by `--min-transfer-distance`.

## 5. Error Handling and Validation

| Situation | Handling |
|-----------|----------|
| Fewer than 3 fields on a line, or a non-numeric frequency | The line is skipped |
| donor / receptor not in the species tree | The transfer is dropped |
| Directory without summaries / no transfer frequency files | **Error** (never silently returns an empty constraint set) |

## 6. Known Limitations

- Transfer frequency files contain no gene tree leaf information, so **no
  `min_family_size` filtering is applied** (AleRax has already aggregated within
  families); the `min_family_size` parameter exists only as a placeholder for interface
  consistency.
- Transfer frequency files are structurally identical to MaxTiC's native
  whitespace-separated constraint format (`donor receptor weight`) — if not passed via
  `--from alerax` and lacking the conventional suffix, they will be parsed as generic
  text constraints (semantically equivalent).

## 7. Usage

```bash
maxtic-next species.tree alerax_run/ --from alerax          # root directory
maxtic-next species.tree run/reconciliations/summaries/ --from alerax
maxtic-next species.tree famA_transfers.txt famB_transfers.txt --from alerax
```

```python
from maxtic_next.constraints.adapters import convert_from_alerax
cset = convert_from_alerax(species_tree, ["alerax_run/"], min_support=0.05)
```

## 8. Samples

`tests/data/alerax_example/run/reconciliations/summaries/FAM{1,2}_transfers.txt` +
`species.tree`. See the regression test `tests/test_alerax_adapter.py`.
