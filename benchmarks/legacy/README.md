# Legacy benchmark

`run_benchmark.py` 是 2026-09 的临时基准脚本，仅保留作历史对照。它使用相对
`examples/` 路径，必须从仓库根运行；临时目录通过 `atexit` 清理。

现行基准入口是：

```bash
python -m maxtic_next.benchmarks.suite --max-n 1000 --max-e 6000
```

请优先使用 `maxtic_next.benchmarks.suite`（支持 `--real`、`--repeats`、`--warmup`）。
