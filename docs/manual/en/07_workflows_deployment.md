# 07 · Workflows & deployment

MaxTiC-Next ships container images and wrappers for the two mainstream workflow engines.

> The commands below are aligned word-for-word with the repository's `Dockerfile`,
> `Singularity.def` and `workflows/`. **Note**: this chapter was written in an environment with
> **no** `nextflow`, `snakemake` or `docker` installed, so the engine-side files are checked
> statically (syntax and path resolution) but have **not been executed end to end**. On a real
> cluster, please work through the self-check list in
> 7.8 first.

## 7.1 The two-stage workflow model

Every engine follows the same two-stage model:

```
[Stage 1 · optional] upstream tool output (.uml_rec / .dtl / recphyloxml / ...)
                 │  maxtic-next ... --from <tool> -o constraints.tsv
                 ▼
              unified constraint file constraints.tsv (5-column ALE comma format)
                 │  maxtic-next species.tree constraints.tsv --ls ...
                 ▼
[Stage 2]     three output files (+ a 4th when --random-trees>0) + HTML report
```

- With text constraints already in hand, start directly at Stage 2.
- Stage 1 products are reusable and cacheable (`--ale-cache-dir` for resumption; the cache key
  includes a species-tree topology fingerprint).
- `-o` is **only meaningful in adapter mode**: passing `-o` with text constraints is a hard
  error instead of a silent no-op.
- Stage 1 and Stage 2 are **two separate process invocations**: `-o` writes the constraints and
  returns without ranking.

## 7.2 Docker

The image is based on `python:3.11-slim`; the entrypoint is `maxtic-next` (the CLI exit code
becomes the container exit code).

```bash
# build (at the repository root; name/tag identical to nextflow.config, Singularity and ch.01)
docker build -t maxtic-next:0.1.1 .

# basic run (mount a data directory)
docker run --rm -v "$PWD/examples":/data maxtic-next:0.1.1 \
    /data/minitree.tree /data/Cyano_CUTConstraints.tsv --seed 42

# enable MCMC (preliminary; convergence diagnostics not validated — samples are not posteriors)
docker run --rm -v "$PWD/examples":/data maxtic-next:0.1.1 \
    /data/minitree.tree /data/Cyano_CUTConstraints.tsv --seed 42 \
    --mcmc --mcmc-iters 2000
```

Notes:

- **Image names must be lowercase**: `docker build -t MaxTiC-Next .` fails outright with
  `invalid reference format`. Use
  `maxtic-next:0.1.1` consistently.
- `jinja2` + `plotly` are installed in the image, so the HTML report works out of the box.
- **`.dockerignore`**: the `Dockerfile` uses `COPY . /app`, so `.dockerignore`
  excludes `.git`, `tests`, `__pycache__`, `.pytest_cache`, `.DS_Store`,
  `benchmarks/datasets` and `docs`. Everything `pip install .` and the in-image smoke test
  actually need is kept: `pyproject.toml`, `src/`, `examples/`, `requirements*.txt`, `LICENSE`,
  `README.md`. No build or run command needs adjusting — the image is simply
  smaller and carries no repository metadata.
- **Writing into a mounted directory**: if a previous run left products there, the CLI refuses
  to overwrite and exits with code **3**. Add `-f/--force`, or change the prefix
  with `-p /data/run2`.
- Smoke test: `docker run --rm maxtic-next:0.1.1 --version` → `MaxTiC-Next 0.1.1`.

## 7.3 Singularity / Apptainer

For HPC clusters (rootless). The build context is the directory containing
`Singularity.def`, i.e. **the repository root**:

```bash
# build at the repository root
sudo singularity build maxtic-next_0.1.1.sif Singularity.def
# or: apptainer build maxtic-next_0.1.1.sif Singularity.def

# the image ships examples/, so this works as-is
singularity run maxtic-next_0.1.1.sif \
    examples/minitree.tree examples/Cyano_CUTConstraints.tsv --seed 42 --no-html

# real data: bind-mount a host directory
singularity run -B "$PWD/data":/data maxtic-next_0.1.1.sif \
    /data/minitree.tree /data/Cyano_CUTConstraints.tsv --seed 42
```

`%files` is an explicit whitelist (unrelated to Docker's `.dockerignore`) and must include
`pyproject.toml`, `README.md` (which the `readme` field in `pyproject.toml` reads), `src/` and
`examples/` — otherwise `pip install .` or the smoke test fails.

## 7.4 Snakemake

Wrapper files: `workflows/Snakefile` + `workflows/config.yaml`
(conda environment definition: `workflows/envs/maxtic_next.yaml`).

Key `config.yaml` fields (named after the CLI options):

```yaml
species_tree: "examples/minitree.tree"
constraints: "examples/Cyano_CUTConstraints.tsv"
seed: 42
from_ale: false          # enables Stage 1 (ALE .uml_rec -> constraints) when true
local_search: 0
min_transfer_distance: 0
threshold_constraints: 0
random_trees: 0          # >0 produces the 4th file *.mt.random_dist.tsv
output_style: "short"    # the CLI default; legacy = the original long names
use_conda: false         # when true, attaches envs/maxtic_next.yaml to the rules' conda:
outdir: "maxtic_next_out"
ale:
  rec_files:
    - "examples/adapters/ale/rec.uml_rec"
  ale_cache_dir: "ale_cache"
  ale_source: "trf"      # trf = default (ALE's own MaxTiC integration) / rec = reconciliations
```

Running (**both working directories work**):

```bash
# from the repository root
snakemake -s workflows/Snakefile --configfile workflows/config.yaml --cores 4
# or from inside workflows/
cd workflows && snakemake --configfile config.yaml --cores 4
```

- **Paths resolve against the repository root**: path resolution is anchored to the repository
  root rather than to the current working directory, so "run from
  inside `workflows/`" still finds `examples/...` instead of looking for `workflows/examples/...`.
  Relative paths in `config.yaml` are evaluated against the **repository root**; absolute
  paths pass through.
- `from_ale: false`: runs Stage 2 only (rank text constraints).
- `from_ale: true`: runs `generate_constraints` (ALE → constraints) first, then ranks.
- Products (three files + HTML, plus one more when `random_trees > 0`) are written to `outdir`;
  the `rule all` target set stays in sync with `--output-style`.
- Both rules pass `--force`: `outdir` belongs to the workflow manager, so a rerun must be able
  to overwrite previous products (otherwise the CLI's overwrite protection fails the task with
  exit code 3).
- With `use_conda: true`, Snakemake creates `envs/maxtic_next.yaml` per rule. That environment
  **includes `jinja2` / `plotly`**, not just an interpreter, so the
  HTML report is produced inside the pipeline instead of silently degrading into
  "the pipeline succeeded but there is no report".

> Note: the snakemake engine is not shipped with the repository; `pip install snakemake` it
> yourself.

## 7.5 Nextflow

Wrapper files: `workflows/main.nf` + `workflows/nextflow.config`.

```bash
# from the repository root
nextflow -C workflows/nextflow.config run workflows/main.nf \
    --species_tree examples/minitree.tree \
    --constraints examples/Cyano_CUTConstraints.tsv \
    --seed 42 --outdir results/
```

ALE two-stage mode (the species tree must match ALE's species naming):

```bash
nextflow -C workflows/nextflow.config run workflows/main.nf \
    --species_tree examples/adapters/ale/species.tree --from_ale \
    --ale_rec_files 'examples/adapters/ale/rec.uml_rec'
```

- The species tree and constraint files are staged into the task directory via channel inputs;
  the parameter is `--species_tree` (mapping to `params.species_tree` in `main.nf`), passed
  through to `maxtic-next` one-for-one. Suited to cloud / cluster batch runs.
- **`--output-style` defaults to `short` (matching the CLI)**, and `RANK`'s output declarations
  are **derived from** `params.output_style` instead of being hard-coded long names; `main.nf`
  and `workflows/config.yaml` carry the same `short` default, so the declared names and the
  names actually written cannot drift apart — a drift would fail the pipeline with
  "output file not found". Pass `--output_style legacy` for the original long names.
- `RANK`'s declarations include the **optional 4th file** (present only with
  `--random-trees > 0`), expressed with `optional: true`, so `--random_trees 0` does not fail
  the pipeline.
- Also passes `--force` (same reason as above).
- To run via Docker: `docker build -t maxtic-next:0.1.1 .` first, then set `docker.enabled` to
  true in `nextflow.config` and uncomment `image = 'maxtic-next:0.1.1'`.

DSL2 structure:

- Both `process` blocks are declared at **top level**; `workflow` only uses
  `if (params.from_ale)` to decide whether to *call* Stage 1. A `process` declared inside an
  `if` block is something DSL2 refuses to parse.
- `GENERATE_CONSTRAINTS` takes `tuple path(species), path(rec)` and is invoked as
  `GENERATE_CONSTRAINTS(species_ch, rec_ch)`, **in the same order** — passing
  `(rec_ch, species_ch)` swaps the arguments and the task fails immediately.
- The species tree is shared by both stages, so it is carried by `Channel.value(...)`
  (re-subscribable, broadcasts automatically) rather than the single-subscription
  `Channel.fromPath` hot channel.

## 7.6 Environment definition (conda)

`workflows/envs/maxtic_next.yaml` defines the conda environment (named `MaxTiC-Next`), for
Snakemake's `conda:` directive (`use_conda: true`) or manual creation. The package is not yet on
PyPI, so it must still be installed after creating the environment (the core needs no
third-party dependency):

```bash
micromamba env create -f workflows/envs/maxtic_next.yaml
micromamba run -n MaxTiC-Next pip install -e <repository root>
```

The environment already declares `jinja2 >=3.0` and `plotly >=5.0`
(equivalent to `pip install "MaxTiC-Next[report]"`).

## 7.7 Bulk-deployment recommendations

| Scenario | Recommendation |
|----------|----------------|
| Thousands of gene families → constraints | In Stage 1 use `--ale-cache-dir` (resumable) and `--ale-parallel process` (multi-core) |
| Retuning ranking parameters | Freeze the Stage 1 product `constraints.tsv` and rerun only Stage 2 |
| Long local search | incremental scoring is on by default (3.1–8.2× faster) + `--checkpoint` against interruption (see chapter 08) |
| Search results that must match across machines | `--local-search-max-iters N` (a wall-clock budget is not portable) |
| Cluster reproducibility | Fix `--seed`, pin the Python/dependency versions, containerize and **pin the image tag** (`:0.1.1`, not `:latest`) |
| Suspicious upstream output | Run `--dry-run` first (in adapter mode it genuinely parses the output, so zero-constraint, hit-rate **and parse-failure** cases all come out as errors and exit 1), then read the `[tool] …` diagnostic line on stderr |
| Upstream examples arrive as archives | `.gz` / gzip streams / single-member archives decompress **transparently**; an official `.tgz` bundle (e.g. ALE's `reconciliations.tgz`) can be fed whole in the constraint position — it is expanded into its members and recorded in `run_metadata["archives_expanded"]` (see 05 §5.7) |
| Wider near-optimal neighbourhood needed | `--near-optimal-top-k K` (default 50) enlarges the support set of the robustness/sensitivity summary; memory and summary cost grow like K·n² (see 03 §3.6b) |
| Exit codes | `0` success / `1` preflight error (including an adapter that cannot parse in `--from` mode) / `2` invalid argument / `3` overwrite refused (`--force` missing) / `4` precondition unmet. Non-zero is captured correctly for Snakemake/Nextflow retries |

> The desktop app `maxtic-studio` (separate repository `MaxTiC-Next-Studio`) reuses the same
> exit-code table, where **`3`** has one extra meaning: with PySide6 missing it prints an install
> hint and exits 3 rather than raising `ModuleNotFoundError` (see 03 §3.10). Workflow engines only
> invoke `maxtic-next`, so they never see that entry point.

## 7.8 Self-check list after these changes

After editing the wrappers on a machine without the engines, verify as follows:

```bash
# 1) image builds, entrypoint works, exit code propagates
docker build -t maxtic-next:0.1.1 .
docker run --rm maxtic-next:0.1.1 --version
docker run --rm maxtic-next:0.1.1 examples/minitree.tree \
    examples/Cyano_CUTConstraints.tsv --seed 42 --no-html; echo "exit=$?"

# 2) the image really has no .git / tests / docs
docker run --rm --entrypoint sh maxtic-next:0.1.1 -c 'ls /app'

# 3) Nextflow syntax (lint even if you do not run it)
nextflow lint -C workflows/nextflow.config workflows/main.nf
nextflow -C workflows/nextflow.config run workflows/main.nf --outdir /tmp/nf_check

# 4) Snakemake dry-run (must work from both working directories)
snakemake -s workflows/Snakefile --configfile workflows/config.yaml -n -p
cd workflows && snakemake --configfile config.yaml -n -p
```

> Prev: [06 · Upstream integration](06_upstream_integration.md) ｜ Next: [08 · Advanced features](08_advanced_features.md)
