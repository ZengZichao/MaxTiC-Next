# Legacy Benchmark

> 🌐 中文版：[README.md](README.md)

`run_benchmark.py` is an ad-hoc benchmark script from 2026-09, retained for
historical comparison only. It uses relative `examples/` paths, so it must be
run from the repository root; the temporary directory is cleaned up via `atexit`.

The maintained benchmark entry point is:

```bash
python -m maxtic_next.benchmarks.suite --max-n 1000 --max-e 6000
```

Please prefer `maxtic_next.benchmarks.suite` (supporting `--real`, `--repeats`,
`--warmup`).
