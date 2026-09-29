# MaxTiC-Next Multi-Tool Compatibility Adapter Framework

> 🌐 中文版：[README.md](README.md)

> This directory contains the specification and API documentation of the **unified constraint
> adapter framework**. It converts the outputs of 5 upstream reconciliation / transfer-inference
> tools into MaxTiC-Next's internal constraint model `ConstraintSet`, and provides
> **plug-in style registration** and **auto-detection of output formats**. For the detailed
> adaptation specification of each tool, see the corresponding `docs/adapters/<tool>.md`.

## 1. Why an Adapter Framework Is Needed

The MaxTiC core consumes only one abstraction: a set of **weighted directed temporal
constraints** `donor → receptor` ("donor is not later than receptor"). Upstream tools,
however, emit heterogeneous output formats — from NHX-style `.uml_rec` files and
human-readable reconciliation reports, to standard recPhyloXML, and on to pre-aggregated
transfer frequency tables. The adapter framework's job is to map these heterogeneous
outputs, **without losing semantics**, onto the single `Constraint(donor, receptor, weight,
metadata)` abstraction.

Key design principles:

- **Semantic mapping, not regex guesswork**: every adapter is implemented on the basis of
  **verified source code / sample inspection of the tool's real output format**, not on
  guessed formats: the RANGER-DTLx / ecceTERA parsers are implemented against the
  upstream source code and real sample output, not against invented formats.
- **Endpoint validity checking**: all donors / receptors must exist in the species tree
  label set (leaf names + internal bootstrap labels); non-matching endpoints are dropped,
  and **incorrect constraints are never silently produced**.
- **Drops must be visible**: every dropped transfer is **counted by reason**; after the
  conversion finishes, a one-line diagnostic summary is printed to **stderr**, and the same
  dictionary is attached to `ConstraintSet.diagnostics` (as
  `result.run_metadata["adapter_diagnostics"]` on the Python API side). An endpoint hit
  rate of **0** raises `LabelMismatchError` directly; below `min_endpoint_hit_rate`
  (default 0.5) a warning is emitted. `--quiet-adapters` **only silences printing**;
  counting and error raising are unaffected.
- **Donor endpoint levels must not be mixed**: ALE takes `parent(donor)`
  (`donor_endpoint_convention = "parent_of_donor"`), while RANGER-DTLx / ARTra / ecceTERA /
  AleRax take the **donor itself** (`"donor_itself"`) — the two conventions **differ by one
  species level** for the same event. When `--from-auto` encounters mixed input, it resolves
  each group with its own adapter and warns that "endpoint levels and weight semantics are
  not commensurable", instead of adding both directly onto the same edge.
- **Weight semantics are not commensurable**: `weight_semantics` is one of four values —
  `support_over_declared_samples` (denominator = the number of samples **declared by the
  tool itself**, e.g. ALE's `N reconciled`, RANGER/ARTra's
  `Total number of optimal solutions: N`), `support_over_replicate_blocks` (valid only when
  all blocks share an identical leaf set, i.e. genuine repeated reconciliations of the same
  family), `integer_count` (denominator cannot be established → the weight is an **integer
  count**, in which case the `min_support` threshold **has no [0,1] support semantics**),
  and `posterior_frequency` (AleRax pre-aggregates upstream; this package no longer
  normalizes).
  ⚠️ **"Number of reconciliation blocks" is not the number of samples**: RANGER-DTLx /
  ARTra print a block header per input **gene tree**, and ARTra even prints
  `Total number of optimal solutions: 24` while outputting only 1 solution — the old
  "block count = support denominator" assumption was therefore wrong. The implementation
  now prefers the sample count declared by the tool, and honestly labels the weight as
  `integer_count` when it cannot be detected.
- **Unified metadata**: every constraint carries `metadata = {family, support, distance}`,
  where `distance` is the species tree topological distance used for
  `--min-transfer-distance` filtering. **Self-loop exception**: events with
  `donor == receptor` are written with `distance = None` (no distance information) in
  **all five** adapters (RANGER-DTLx / ecceTERA / ARTra / AleRax / ALE), so they pass the
  "constraints without a distance column are always kept" rule and the `to itself`
  statistic is no longer swallowed by `--d`; the text-input path has always behaved this
  way. ALE was the last to be aligned; see the regression test
  `tests/test_ale_selfloop.py`, and [ale.md](ale.md) plus §3.2.3 of
  `docs/manual/en/03_cli_reference.md` for symptoms and measurements.

## 2. Supported Tools at a Glance

| Tool | Input | Real format (verified source) | Transfer semantics |
|------|-------|-------------------------------|--------------------|
| **ale** | `*.uml_rec` | Sampled reconciliations, NHX `T@donor->receptor` / `D@` (`constraints_from_reconciliations.py`) | `T@` events donor→receptor |
| **ranger** | Reconciliation report text | `m# = LCA[..]: Transfer, Mapping --> D, ..., Recipient --> R` (`DTL-algorithm.h:1840`) | Mapping=donor, Recipient=receptor |
| **eccetera** | recPhyloXML | `branchingOut speciesLocation=D` + child `transferBack destinationSpecies=R` (`DTLGraph.cpp`) | Parent branchingOut=donor, child transferBack=receptor |
| **artra** | `output.txt` | Same shape as RANGER, with a `Replacing/Additive Transfer` prefix (`ARTra-程序/output.txt`) | Mapping=donor, Recipient=receptor |
| **alerax** | `reconciliations/summaries/*_transfers.txt` or the output directory | `donor receptor frequency` (`extract_families_transfer.py`) | Each line donor→receptor, frequency as weight |

> RANGER-DTLx and ARTra outputs are **identical in shape** and share the parsing core
> `_recon_report`; therefore `--from ranger` and `--from artra` produce consistent results
> on the same report (ARTra additionally supports filtering by transfer class).

## 3. Unified API

### 3.1 Registry dispatch

```python
from maxtic_next.tree.tree import Tree
from maxtic_next.constraints.adapters import registry

species = Tree(); species.read_newick(open("species.tree").readline())

# Explicitly specify the tool
cset = registry.convert("ranger", species, ["FAM1.dtl"], min_family_size=0)

# Auto-detection (by content / extension)
cset = registry.convert_auto(species, ["out.recphyloxml"])

print(registry.REGISTRY.names())      # ['ale', 'alerax', 'artra', 'eccetera', 'ranger']
print(registry.detect_format("run/")) # 'alerax'
```

### 3.2 Convenience functions (per tool)

```python
from maxtic_next.constraints.adapters import (
    convert_from_ale, convert_from_ranger_dtl, convert_from_eccetera,
    convert_from_artra, convert_from_alerax,
)
cset = convert_from_eccetera(species, ["FAM1.recphyloxml"], min_family_size=5)
```

### 3.3 Top-level `api.rank` (end-to-end)

```python
from maxtic_next import api
# Rank directly from upstream tool output (one step)
api.rank("species.tree", ["FAM1.dtl"], from_tool="ranger")
# Two-stage: only generate the unified constraint file (for reuse by Snakemake/Nextflow), no ranking
api.rank("species.tree", ["run/"], from_tool="alerax", constraints_out="constraints.tsv")
```

### 3.4 Common parameters

| Parameter | Meaning | Default | Applies to |
|-----------|---------|---------|------------|
| `min_support` | Minimum support / frequency threshold within a single family; keeps `> min_support` | 0.05 | all |
| `min_family_size` | Minimum gene family size (leaf count); families with `<=` the threshold are skipped | 5 | ale/ranger/eccetera/artra (not applicable to alerax) |
| `cache_dir` | File-level pickle cache directory (resumable runs) | None | ale/ranger/eccetera/artra |
| `parallel` | `"process"` (process pool, bypasses the GIL) / `"thread"` | `"process"` | all |
| `max_workers` | Number of parallel workers | auto | all |
| `output_path` | Write out unified constraints (`family,donor,receptor,weight,distance`) | None | all |
| `source` | ALE constraint source `"rec"`/`"trf"` | **`"trf"`** (`DEFAULT_SOURCE`; the convention of ALE's official MaxTiC integration.) | ale |
| `transfer_kind` | ARTra transfer class `"all"`/`"replacing"`/`"additive"` | `"all"` | artra |
| `min_endpoint_hit_rate` | Endpoint hit-rate warning threshold; a hit rate of 0 raises an error directly (`<= 0` disables warnings) | `0.5` | all |
| `quiet` | Do not print the stderr diagnostic summary (counting and error raising are unaffected) | `False` | all |

## 4. Command Line

```bash
# Explicit tool
maxtic-next species.tree FAM1.dtl        --from ranger
maxtic-next species.tree out.recphyloxml --from eccetera
maxtic-next species.tree out.txt         --from artra --artra-transfer-kind replacing
maxtic-next species.tree alerax_run/     --from alerax
maxtic-next species.tree rec1.uml_rec rec2.uml_rec --from ale   # equivalent to --from-ale

# Auto-detection
maxtic-next species.tree some_output      --from-auto

# Two-stage: only produce constraints (for pipelines)
maxtic-next species.tree FAM1.dtl --from ranger -o constraints.tsv
```

`--ale-min-support` / `--ale-min-family-size` / `--ale-cache-dir` / `--ale-parallel` apply
to all reconciliation adapters (the historical `--ale-` prefix is kept for backward
compatibility).

## 5. Auto-Detection Rules

`detect_format(path)` order (from most specific to most permissive):

1. **Directory**: contains `reconciliations/summaries` or `*_transfers.txt` → `alerax`;
2. **Extension**: `.uml_rec` → ale; `_transfers.txt` → alerax; `.recphyloxml` → eccetera;
3. **Content sniffing** (first 64 KB): contains `<recPhylo`/`<recGeneTree` → eccetera;
   contains `T@` and `reconciled`/`D@` → ale; contains `Recipient -->` and
   `Replacing/Additive Transfer` → artra; contains only `Recipient -->` → ranger;
4. No match → `None` (treated as generic text constraints and parsed by the native path).

`detect_formats(paths)`: multiple paths must resolve to the **same** tool, otherwise a
`ValueError` is raised (mixed input requires an explicit `--from <tool>`).

## 6. Extension: Registering a Custom Adapter

```python
from maxtic_next.constraints.adapters.registry import REGISTRY, AdapterEntry

def convert_mytool(species_tree, inputs, **kwargs):
    ...  # returns a ConstraintSet
    return cset

REGISTRY.register(AdapterEntry(
    "mytool", convert_mytool, "My reconciliation tool", "*.myrec files"))
```

Note: the content-sniffing rules for auto-detection are built into `detect_format`; if a
custom tool needs auto-detection, decide in the application layer first and then dispatch
explicitly with `registry.convert("mytool", ...)`.

## 7. Known Limitations (General)

- **Internal node label alignment**: constraint endpoints must match the species tree
  labels that the user passes into MaxTiC. ecceTERA labels internal species nodes with its
  own numeric IDs (see `docs/adapters/eccetera.md`); if they do not match the species tree
  bootstrap labels they are filtered out by `valid_labels` — this is **conservatively
  correct** (no false constraints are produced), but may reduce recall.
- **Family size estimation**: the family size in reconciliation reports / recPhyloXML is
  estimated from leaf event counts, not the exact gene tree leaf count.
- **Cache trust**: the file-level cache is a local pickle (consistent with ALE) and only
  reads trusted files created by the current run inside `cache_dir`; do not point it at
  untrusted sources.
