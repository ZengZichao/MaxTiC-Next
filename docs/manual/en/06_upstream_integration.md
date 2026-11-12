# 06 · Upstream integration recipes

This chapter is the end-to-end integration manual between MaxTiC-Next and the **five upstream
reconciliation / transfer-inference tools**. The core idea:

> Upstream tool (infers HGT transfer events) → an **adapter** unifies them into
> `Constraint(donor, receptor, weight, metadata)` → MaxTiC-Next ranks. You can either do it in
> one shot (`--from <tool>` ranks directly) or in two stages (`-o` writes a unified constraint
> file first, then rank).

## 6.0 One mental model

Every adapter emits the same isomorphic `Constraint(donor, receptor, weight, metadata={family, support, distance})` shape:

- **donor → receptor**: always means "donor no later than receptor";
- **Endpoint validity**: donor/receptor must exist in **the species tree you pass to MaxTiC**
  (internal bootstrap labels + leaf names), otherwise the constraint is dropped (**never silently invents a wrong constraint**);
- **Weight/support**: typically "sample count / #samples" or a posterior frequency;
- **Distance**: the species-tree topological `distance_from(receptor, donor)`, used by `--d`.

Common parameters (the `--ale-` prefix is historical, but they do not take effect for every
adapter — two exceptions, see below):

| Option | Meaning | Default |
|--------|---------|---------|
| `--ale-min-support` | keep `support > threshold` | `0.05` |
| `--ale-min-family-size` | skip families with leaves `<=` threshold (n/a for alerax) | `5` |
| `--ale-cache-dir` | file-level pickle cache (resume) | none |
| `--ale-parallel` | `process` (default) / `thread` | `process` |

> **Two exceptions** (same account as 03 §3.2): the `--min-endpoint-hit-rate` threshold guard
> applies to `ranger/eccetera/artra/alerax` only — it has **no effect for `ale`**, which keeps the
> original tool's endpoint-check semantics and never warns however low the hit rate is
> (`ale.py:842-851` hard-codes `threshold=0.0, raise_on_zero_hit=False`); and
> `--ale-min-family-size` / `--ale-cache-dir` are mere interface placeholders for `alerax`
> (`alerax.py:222-223` assigns them and nothing ever reads them), so passing them is accepted but
> **has no effect**.

> **Small-sample tip**: tiny demo families are easily filtered by `--ale-min-family-size 5`; set
> `--ale-min-family-size 0` temporarily. Use the default (or a tuned value) for real data.

Auto-detection (`--from-auto`) runs three passes, from the most specific rule to the loosest
(`registry.py:288-320`): (1) **directory** — contains a `reconciliations/summaries` subdirectory,
or any `*_transfers.txt` / `*_meanTransfers.txt` inside → alerax; (2) **extension** (compression
suffix stripped first) — `.uml_rec` → ale, `*_transfers.txt` / `*_meanTransfers.txt` → alerax,
`.recphyloxml` → eccetera; (3) **content sniffing** — `<recPhylo` → eccetera, `T@` together with
`reconciled`/`D@` → ale, `Recipient -->` with `Replacing/Additive Transfer` → artra, only
`Recipient -->` → ranger. When the inputs resolve to **several** tools, `--from-auto` groups
them, converts each group and merges them with a warning; it is an explicit `--from` that
requires a single tool.

---

## 6.1 ALE (`ALEml_undated` / `ALEmcmc_undated`)

**Role**: Amalgamated Likelihood Estimation; sampled reconciliation per gene family, output `.uml_rec`.

**Real format**: a `<N> reconciled …` declaration line (ALE_undated actually prints
`<N> reconciled G-s:`, and the repository sample `examples/adapters/ale/rec.uml_rec` starts with
`1 reconciled`), followed by N NHX-style sampled reconciled trees; events annotate node labels,
`T@donor->receptor` (transfer), `D@species` (duplication). The adapter only tests whether the
**second token is exactly `reconciled`** (`ale.py:261-262`) — the rest of the declaration line plays
no part. The first token is `N`, used as the support denominator; when it is missing the
denominator falls back to the number of trees actually parsed, and if neither is available the
weight is an **integer count**, declared in the diagnostics and on stderr (`ale.py:615-661`).

**Upstream commands (illustrative)**:

```bash
# 1. observe & reconcile (per gene family)
ALEobserve gene_family_XXX.treelist
ALEml_undated species_tree.nwk gene_family_XXX.treelist.ale
# yields gene_family_XXX.treelist.ale.uml_rec

# 2. collect all .uml_rec
ls *.uml_rec > all_rec_files    # or just list them as arguments
```

**MaxTiC-Next recipe (one shot)**:

```bash
maxtic-next species_tree.nwk fam1.ale.uml_rec fam2.ale.uml_rec ... \
    --from ale --ale-source rec \
    --ale-min-support 0.05 --ale-min-family-size 5 \
    --ale-cache-dir ale_cache/
```

**Two stages (recommended for large batches / pipelines)**:

```bash
# Stage 1: generate the unified constraint file only (reusable, cacheable/resumable)
maxtic-next species_tree.nwk *.uml_rec --from ale \
    --ale-cache-dir ale_cache/ -o constraints.tsv
# Stage 2: rank the constraints
maxtic-next species_tree.nwk constraints.tsv --ls 180
```

**Key points**:
- `--ale-source` defaults to **`trf`**: each `T@donor->receptor` event
  yields `parent(donor) -> receptor`, which is the convention of **ALE's own MaxTiC
  integration** (`constraints_from_transfers`). Its full statement is "**the parent of the donor
  branch** is older than **the child of the receptor branch**" — the receptor endpoint stands for
  the child of the receptor branch, not the receptor branch itself (same wording as the
  `--ale-source` help text and 03 §3.2.1). `rec` draws on the reconciliation-record block
  of `.uml_rec` (in the upstream script that block sits inside the `T@` gate, so pure
  speciation/duplication nodes contribute nothing). **These two conventions are not
  interchangeable**: constraint **counts and total weights both differ** (on the official example
  `rec` yields only 13,995 constraints with a total weight of 5482.98, vs 6438.6 for `trf` —
  about 17% apart; see [../../adapters/ale.en.md](../../adapters/ale.en.md) §2), so switching
  replaces the entire input. Pass `--ale-source rec` explicitly to reproduce results computed
  from reconciliation records. The effective convention is printed to stderr.
- The donor endpoint convention is **`parent_of_donor`** (the donor's parent) — it differs by
  **one species level** from `donor_itself` used by RANGER-DTLx / ARTra / ecceTERA / AleRax.
  With mixed `--from-auto` input the files are grouped per adapter and a warning is emitted,
  instead of being silently summed. The merge runs **two** consistency checks
  (`registry.py:464-483`): different endpoint levels are recorded as `convention_conflict`, and a
  different **weight convention** (`weight_semantics` — `[0,1]` support, integer counts and
  posterior frequencies are incommensurable) as `semantics_conflict`; either one warns on stderr
  and is written into the diagnostics.
- The receptor must not be a leaf (extant species); donor/receptor must be in the species-tree
  label set.
- In-family **reciprocal subtraction** removes symmetric noise; `support = count / N`, where N
  is the `N reconciled` figure the report itself declares. The keep test is
  `support > threshold`, so constraints **at or below** `--ale-min-support` (default 0.05) are
  dropped and counted.
- The `--ale-cache-dir` cache key includes a **species-tree topology fingerprint**:
  two trees sharing a label set but differing in topology never hit the same
  cache entry.
- Details: [../../adapters/ale.en.md](../../adapters/ale.en.md).

---

## 6.2 RANGER-DTLx

**Role**: DTL reconciliation + AggregateRanger; a **human-readable** reconciliation report.

**Real format** (source `DTL-algorithm.h:1840`): per gene node; transfer lines look like

```
 = LCA[g1, g4]: Transfer, Mapping --> 61, Edge, Parent = 65, Recipient --> 62
```

`Mapping --> D` = donor, `Recipient --> R` = receptor → constraint `(D, R)`. One file may contain
multiple `Reconciliation for Gene Tree N` blocks, summed per edge key across blocks. **The block
count is not automatically the support denominator**: the denominator is settled in three tiers
(`_recon_report.py:243-279`) — the sample count the report declares itself (a line such as
`Total number of optimal solutions: N`, the maximum when several are declared) wins; with no
declaration, blocks are used only when there is **more than one block** and all blocks share the
**same leaf set**; otherwise there is no denominator, the weight degrades to an **integer count**,
and a warning is emitted to stderr.

**Upstream commands (illustrative)**:

```bash
Ranger-DTL -i gene_family_XXX.nwk -s species_tree.nwk -o FAM_XXX.dtl
# multiple optima / samples may be merged as multiple blocks in one file,
# or split across files and summed across files by the adapter
```

**MaxTiC-Next recipe**:

```bash
# one shot
maxtic-next species_tree.nwk FAM1.dtl FAM2.dtl --from ranger \
    --ale-min-family-size 5

# two stages
maxtic-next species_tree.nwk *.dtl --from ranger -o constraints.tsv
maxtic-next species_tree.nwk constraints.tsv --ls 180
```

**Key point**: internal nodes are often numeric labels (e.g. `61`) and must match the species-tree
bootstrap labels. Details: [../../adapters/ranger_dtlx.en.md](../../adapters/ranger_dtlx.en.md).

---

## 6.3 ecceTERA

**Role**: DTL reconciliation; standard **recPhyloXML** output (Duchemin et al. 2018).

**Real format** (source `DTLGraph.cpp`): donor node `<branchingOut speciesLocation="D">`, and the
receptor child node `<transferBack destinationSpecies="R">` → constraint `(D, R)`.

**Upstream commands (illustrative)**:

```bash
ecceTERA species.tree=species_tree.nwk gene.trees=gene_family_XXX.nwk \
    output.dir=out/ recPhyloXML.reconciliation=true
# yields out/..._recPhyloXML.xml (recPhyloXML)
```

**MaxTiC-Next recipe**:

```bash
maxtic-next species_tree.nwk fam1.recphyloxml fam2.recphyloxml --from eccetera
# or auto-detect (recPhyloXML is recognized as eccetera)
maxtic-next species_tree.nwk out.recphyloxml --from-auto
```

**Key limitation (important)**: ecceTERA labels **internal species nodes** with its own **numeric
IDs** (leaves use names); these IDs may differ from your species-tree bootstrap labels — mismatched
endpoints are conservatively dropped (no spurious constraints, but possibly lower recall). Align
your species-tree internal labels with ecceTERA's naming. XML containing `<!DOCTYPE>`/`<!ENTITY>`
is **rejected** (XXE guard). Details: [../../adapters/eccetera.en.md](../../adapters/eccetera.en.md).

---

## 6.4 ARTra (Additive / Replacing transfer classification)

**Role**: classifies each transfer on a RANGER-style reconciliation as **Additive** or **Replacing**.

**Real format**: **isomorphic** to a RANGER report, only event names gain a `Replacing`/`Additive` prefix:

```
 = LCA[H117, H14]: Replacing Transfer, Mapping --> H15, Recipient --> H125
 = LCA[H159, H17]: Additive Transfer, Mapping --> H17, Recipient --> H159
```

**MaxTiC-Next recipe**:

```bash
# count all transfers (default)
maxtic-next species_tree.nwk FAM1.txt --from artra --artra-transfer-kind all
# only replacing transfers
maxtic-next species_tree.nwk FAM1.txt --from artra --artra-transfer-kind replacing
# only additive transfers
maxtic-next species_tree.nwk FAM1.txt --from artra --artra-transfer-kind additive
```

**Key point**: MaxTiC's time-constraint semantics are independent of additive vs. replacing, so all
transfers are counted by default; the category filter only lets you focus on one class. It shares the
parsing core with RANGER, so results agree on the same report.
Details: [../../adapters/artra.en.md](../../adapters/artra.en.md).

---

## 6.5 AleRax

**Role**: Bayesian gene-tree/species-tree reconciliation & scoring; transfers are **pre-aggregated**
into expected frequencies within each family.

**Real format** (script `extract_families_transfer.py`): one
`reconciliations/summaries/<family>_transfers.txt` per family in the output directory; current
AleRax (`AleOptimizer.cpp`) writes `<family>_meanTransfers.txt` instead, and **both suffixes are
accepted** (`alerax.py:67`, `registry.py:265`). Each line is `donor receptor frequency`.
Frequencies are already expectations, used **directly as weights** (not divided by block count).

**MaxTiC-Next recipe** (three input forms, auto-expanded):

```bash
# 1. AleRax output root dir (auto-appends reconciliations/summaries)
maxtic-next species_tree.nwk alerax_run/ --from alerax
# 2. summaries dir
maxtic-next species_tree.nwk alerax_run/reconciliations/summaries/ --from alerax
# 3. explicit *_transfers.txt / *_meanTransfers.txt files
maxtic-next species_tree.nwk famA_transfers.txt famB_transfers.txt --from alerax
```

**Key point**: neither `*_transfers.txt` nor `*_meanTransfers.txt` carries leaf info, so **no
`min_family_size` filtering**; only keeps `frequency > --ale-min-support` (default 0.05);
frequencies are summed across families per edge key.
Details: [../../adapters/alerax.en.md](../../adapters/alerax.en.md).

---

## 6.6 Five-tool cheat sheet

| Tool | Input | Donor/receptor source | Weight | Family filter | CLI |
|------|-------|-----------------------|--------|---------------|-----|
| ALE | `*.uml_rec` | `T@D->R` / reconciliation-event sets | count ÷ declared N (fallbacks: 6.1) | yes | `--from ale` |
| RANGER-DTLx | reconciliation report | `Mapping-->D` / `Recipient-->R` | count ÷ three-tier denominator (6.2) | yes | `--from ranger` |
| ecceTERA | recPhyloXML | `branchingOut` / `transferBack` | count ÷ declared samples, else blocks | yes | `--from eccetera` |
| ARTra | `output.txt` | same as RANGER | same as RANGER (three-tier) | yes | `--from artra` |
| AleRax | `*_transfers.txt` / `*_meanTransfers.txt` | each line `donor receptor` | frequency | no | `--from alerax` |

## 6.7 End-to-end best practices

1. **Dry-run first**: `--dry-run` checks tree/label uniqueness/endpoint validity (especially ecceTERA ID alignment).
2. **Two stages**: for large family batches, first `-o constraints.tsv` (with `--ale-cache-dir` resume),
   then re-tune ranking (`--ls`/`--ts`/`--d`) on the constraints without re-parsing.
3. **Auto-detect**: use `--from-auto` when unsure; **mixed-tool inputs** are grouped, converted per
   tool and merged with a warning, but the merged weights are not comparable — to keep them clean,
   run one tool at a time with an explicit `--from`.
4. **Distance filtering**: use `--d` only when constraints actually carry a distance column (ALE 5-col / space 4-col).
5. **Reproducibility**: fix `--seed`, record version and parameters (both the HTML report and stdout header carry run metadata).

## 6.8 Custom adapter (extension)

```python
from maxtic_next.constraints.adapters.registry import REGISTRY, AdapterEntry

def convert_mytool(species_tree, inputs, **kwargs):
    ...  # parse your tool's output, return a ConstraintSet
    return cset

REGISTRY.register(AdapterEntry("mytool", convert_mytool, "我的调和工具", "*.myrec"))
```

The third positional argument is `AdapterEntry.description` (`registry.py:95-99`): it is stored
verbatim and echoed by `AdapterRegistry.describe()`, the human-readable tool list meant for CLI
help / documentation (`registry.py:122-131`), and all five built-in registrations keep it in
Chinese (`registry.py:235-267`; e.g. `"AleRax 贝叶斯调和"` = "AleRax Bayesian
reconciliation", `"ecceTERA DTL 调和"` = "ecceTERA DTL reconciliation"). The literal
above therefore stays Chinese — an English string in your snippet would never match what
the tool list prints (*"my reconciliation tool"*).

> Prev: [05 · Input / output formats](05_io_formats.md) | Next: [07 · Workflows & deployment](07_workflows_deployment.md)
