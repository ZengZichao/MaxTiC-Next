# Security Policy

MaxTiC-Next is research software distributed under the CeCILL-2.1 license
without warranty; see [LICENSE](LICENSE).

## Supported versions

| Version | Supported |
| --- | --- |
| 0.2.x | Yes |
| < 0.2.0 | No |

## Reporting a vulnerability

Please report suspected security issues privately via GitHub's
[security advisory form](https://github.com/ZengZichao/MaxTiC-Next/security/advisories/new)
rather than opening a public issue. If that form is unavailable, contact the
corresponding author listed in [CITATION.cff](CITATION.cff).

Include the version or commit, a minimal reproduction, and the expected versus
actual behaviour. Do not attach real data — use a synthetic or minimal example.

## Scope

In scope: the core Python package (`src/maxtic_next/`), the bundled Snakemake
and Nextflow wrappers (`workflows/`), and the distribution artefacts (PyPI,
conda recipe, Docker and Singularity images).

Out of scope: third-party dependencies, which should be reported to their
upstream projects, and the separate desktop front-end repository.
