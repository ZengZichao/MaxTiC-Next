# 01 · Installation

## 1.1 Requirements

| Item | Requirement |
|------|-------------|
| Python | ≥ 3.9 (**3.11 recommended**) |
| Core algorithm deps | standard library only (`argparse` / `random` / `math` / `functools`), no network |
| HTML report (optional) | `jinja2` ≥ 3.0, `plotly` ≥ 5.0 |
| Pruning / graph operations (optional) | no third-party dependency: `prune.py` uses the standard library only, deliberately no ete3 / networkx |
| Dev/test (optional) | `pytest` ≥ 8.0 |

> **Important**: the core ranking functionality needs **no third-party libraries**. Only the
> interactive HTML report requires `jinja2` + `plotly`.

## 1.2 Option A: pip install (recommended)

```bash
# no clone needed up front: install straight from the repository
pip install "git+https://github.com/ZengZichao/MaxTiC-Next.git"

# for the interactive HTML report (jinja2 + plotly), clone and install with extras
git clone https://github.com/ZengZichao/MaxTiC-Next.git
cd MaxTiC-Next
pip install -e ".[report]"
```

Verify:

```bash
maxtic-next --help
python -c "import maxtic_next; print(maxtic_next.__version__)"   # 0.1.0
```

## 1.3 Option B: no install (PYTHONPATH)

Offline / without installing, run directly from source:

```bash
git clone https://github.com/ZengZichao/MaxTiC-Next.git
cd MaxTiC-Next
export PYTHONPATH=src
python -m maxtic_next --help
```

## 1.4 Option C: micromamba / conda (recommended for this project)

All scripts and tests in this project run and are validated in the micromamba `python-3.11`
environment.

```bash
# if the env does not exist yet
micromamba create -n python-3.11 python=3.11 -y
micromamba run -n python-3.11 pip install -e ".[report]"

# then wrap every command with micromamba run:
micromamba run -n python-3.11 maxtic-next examples/minitree.tree \
    examples/Cyano_CUTConstraints.tsv --seed 42
```

## 1.5 Option D: containers (Docker / Singularity)

See [07 · Workflows & deployment](07_workflows_deployment.md). Base image `python:3.11-slim`:

```bash
docker build -t maxtic-next:0.1.0 .
docker run --rm -v "$PWD/examples":/data maxtic-next:0.1.0 \
    /data/minitree.tree /data/Cyano_CUTConstraints.tsv --seed 42
```

> Image names must be lowercase (`MaxTiC-Next` fails with `invalid reference format`); the tag
> matches `Dockerfile`, `Singularity.def` and `workflows/nextflow.config`. The build context is
> narrowed by `.dockerignore` (`.git`, `tests`, `docs`, … stay out of the image).

## 1.6 Run the test suite (optional)

```bash
python -m pytest tests/ -q
# Expected: all pass (currently 486 passed, 1 skipped; covers adapters, the
#           reference-implementation differential, boundary cases and CLI switches.
#           The single skip needs python2 plus the original MaxTiC.py to run)
```

## 1.7 Uninstall

```bash
pip uninstall MaxTiC-Next
```

## 1.8 Common install issues

| Symptom | Cause / fix |
|---------|-------------|
| `maxtic-next: command not found` | not pip-installed; use `PYTHONPATH=src python -m maxtic_next` |
| `ModuleNotFoundError: jinja2` | report deps missing; `pip install -e ".[report]"` or add `--no-html` |
| Report fails but ranking works | the report is optional; failure prints a warning and does not block |
| Path with spaces / non-ASCII fails | quote the path; prefer running in a space-free directory |

> Prev: [00 · Overview](00_index.md) ｜ Next: [02 · Quickstart](02_quickstart.md)
