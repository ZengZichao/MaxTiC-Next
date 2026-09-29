# Multi-Tool Compatibility Adapters — Sample Inputs for Each Upstream Tool

> 🌐 中文版：[README.md](README.md)

This directory provides **real-format** sample inputs for 5 upstream reconciliation /
transfer-inference tools, demonstrating how to feed their output directly into
MaxTiC-Next ranking via `--from <tool>` / `--from-auto`.
The species tree is uniformly `../minitree.tree` (internal node labels are bootstrap
numbers such as 61/62/67).

## Contents

| Subdirectory / File | Tool | Format |
|---------------------|------|--------|
| `ranger/FAM1.dtl` | RANGER-DTLx | `Mapping --> / Recipient -->` reconciliation report |
| `eccetera/FAM1.recphyloxml` | ecceTERA | recPhyloXML (`branchingOut` / `transferBack`) |
| `artra/FAM1.txt` | ARTra | `Replacing/Additive Transfer` report |
| `alerax/run/reconciliations/summaries/*_transfers.txt` | AleRax | `donor receptor freq` |
| `run_adapters.sh` | — | one-command demo of all tools + auto-detection + two-stage |

## One-Command Run

```bash
micromamba run -n python-3.11 bash examples/adapters/run_adapters.sh
```

## Individual Run Examples

```bash
# Explicit tool
maxtic-next examples/minitree.tree examples/adapters/ranger/FAM1.dtl \
    --from ranger --ale-min-family-size 0

# recPhyloXML auto-detection
maxtic-next examples/minitree.tree examples/adapters/eccetera/FAM1.recphyloxml \
    --from-auto --ale-min-family-size 0

# ARTra: count replacing transfers only
maxtic-next examples/minitree.tree examples/adapters/artra/FAM1.txt \
    --from artra --artra-transfer-kind replacing --ale-min-family-size 0

# AleRax directory (automatically expands *_transfers.txt)
maxtic-next examples/minitree.tree examples/adapters/alerax/run \
    --from alerax

# Two-stage: produce only the unified constraint file (for reuse by Snakemake/Nextflow)
maxtic-next examples/minitree.tree examples/adapters/ranger/FAM1.dtl \
    --from ranger -o constraints.tsv --ale-min-family-size 0
```

> Note: these are **small demonstration samples**, hence `--ale-min-family-size 0` is
> added so that the family-size filter does not drop them; for real data, use the
> default `5` or adjust as needed. See `docs/adapters/<tool>.md` for each tool's
> adapter specification and known limitations.
