# Contributing to MaxTiC-Next

Thanks for your interest in contributing! This document covers the development
setup, quality gates, and release conventions. 中文速览见文末。

## Development setup

MaxTiC-Next targets Python 3.9–3.14 and keeps **zero runtime dependencies**
(standard library only). Optional extras are declared in `pyproject.toml`:

```bash
git clone https://github.com/ZengZichao/MaxTiC-Next
cd MaxTiC-Next
python -m pip install -e ".[dev,report]"   # dev toolchain + optional report deps
```

## Running the test suite

```bash
pytest tests/ -v                               # full suite
pytest tests/test_architecture_gates.py -v     # architecture & documentation gates
python tools/coverage_gate.py core             # coverage floor check
```

CI runs the core matrix on Python 3.9–3.14, once with and once without
`numpy>=1.26`. The numpy and pure-stdlib implementations of
`ranking/reachability.py` and `tree/cache.py` must stay output-equivalent —
this is guarded by dedicated equivalence tests.

## Style gates (all blocking in CI)

```bash
ruff check src tests tools benchmarks
ruff format --check src tests tools benchmarks
mypy
```

Rule sets live in `pyproject.toml` (`[tool.ruff.lint]` and mypy overrides).

## Commit & pull request conventions

- Follow Conventional Commits (`feat:`, `fix:`, `docs:`, `chore:`, …), as used
  in the existing history.
- Do **not** add runtime dependencies under `src/` — stdlib-only is a hard
  design goal of the rewrite.
- Open pull requests against `main`. All CI checks (core matrix, lint & types,
  sdist/wheel install) are required before merge; branch protection enforces
  this for everyone, including maintainers.
- Use the issue templates (`bug report` / `feature request`) when opening a
  new issue.

## Release process

1. Bump the version in **both** places: `pyproject.toml` (`version = …`) and
   `src/maxtic_next/__init__.py` (`__version__`).
2. Update `CITATION.cff` (`version`, `date-released`) and mint/update the
   Zenodo DOI if applicable.
3. Tag with an **annotated** tag: `git tag -a vX.Y.Z -m "Release vX.Y.Z"`.
   Annotated tags carry the tagger and message metadata and keep the history
   consistent (v0.1.1 was lightweight by accident — please do not repeat that).
4. Push the tag and publish a GitHub release with notes; CI must be green on
   the release commit.

## License

By contributing, you agree that your contributions are licensed under the
repository license (CeCILL-2.1).

---

## 中文速览

- 开发安装：`pip install -e ".[dev,report]"`。核心算法零运行时依赖，请勿在
  `src/` 引入第三方运行库。
- 提交前自检：`pytest tests/`、`ruff check`、`ruff format --check`、`mypy`、
  `python tools/coverage_gate.py core`。
- 提交信息使用 Conventional Commits；向 `main` 发起 PR，全部 CI 检查通过后
  才能合并（分支保护对所有人生效，含管理员）。
- 发版：同时更新 `pyproject.toml` 与 `__init__.py` 中的版本号、`CITATION.cff`，
  使用 annotated tag（`git tag -a vX.Y.Z`），CI 绿灯后再发 GitHub Release。
- 详细规则以英文正文为准。
