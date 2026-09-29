#!/usr/bin/env python3
"""覆盖率地板门禁。

为什么不用 `coverage --fail-under` 一把梭：核心算法是零依赖的，任何作业都必须
能单独断言它的地板；把阈值写死在配置里，会被可选依赖作业的执行情况稀释成
一个既不严格也不可信的总数。

图形前端位于独立仓库 MaxTiC-Next-Studio，其覆盖率由该仓库自己的 CI 度量，
本仓库只有一个作用域：

    python tools/coverage_gate.py           # 全仓 ≥ 85%
    python tools/coverage_gate.py core      # 同上

前置：同目录已存在 `.coverage` 数据文件（由 `pytest --cov=maxtic_next` 产生）。
本脚本只用标准库，读 `coverage json` 的输出。
"""

from __future__ import annotations

import json
import os
import subprocess
import sys
import tempfile

CORE_FLOOR = 85.0  # 改前实测 85.7%（ data/coverage-summary.md）

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)


def _load_coverage_json() -> dict:
    """跑 `coverage json` 导出明细（不依赖 pyproject 的 omit）。

    中间 json 写到系统临时目录，避免在仓库里留下未跟踪垃圾文件。
    """
    out = os.path.join(tempfile.gettempdir(), "maxtic-coverage-scope.json")
    cmd = [sys.executable, "-m", "coverage", "json", "-o", out, "-q"]
    proc = subprocess.run(cmd, capture_output=True, text=True, cwd=ROOT)
    if proc.returncode != 0 or not os.path.isfile(out):
        sys.exit(
            "找不到覆盖率数据：请先运行 `pytest --cov=maxtic_next`（"
            f"coverage json 返回：{proc.stdout.strip()} {proc.stderr.strip()}）"
        )
    with open(out, encoding="utf-8") as fh:
        return json.load(fh)


def _ratio(data: dict):
    stmts = missing = 0
    files = []
    for path, info in data["files"].items():
        norm = path.replace("\\", "/")
        summary = info["summary"]
        stmts += summary["num_statements"]
        missing += summary["missing_lines"]
        files.append(
            (
                summary["percent_covered"],
                norm.rsplit("maxtic_next/", 1)[-1],
                summary["num_statements"],
            )
        )
    if not stmts:
        return None, None, []
    return (stmts - missing) * 100.0 / stmts, stmts, sorted(files)


def main(argv: list) -> int:
    scope = (argv[1] if len(argv) > 1 else "core").lower()
    if scope not in ("core", "all"):
        sys.exit(f"未知作用域 {scope!r}（本仓库仅有：core）")
    ratio, stmts, files = _ratio(_load_coverage_json())
    if ratio is None:
        print("覆盖率数据为空")
        return 1
    mark = "OK  " if ratio >= CORE_FLOOR else "FAIL"
    print(f"  {mark} {'core':12s} {ratio:5.1f}% / 地板 {CORE_FLOOR:.0f}%（{int(stmts)} 条语句）")
    for pct, name, n in [f for f in files if f[0] < CORE_FLOOR][:8]:
        print(f"       {pct:5.1f}%  {n:4d} stmt  {name}")
    if ratio < CORE_FLOOR:
        print(f"\n覆盖率门禁未通过：core 覆盖率 {ratio:.1f}% 低于地板 {CORE_FLOOR:.0f}%")
        return 1
    print("覆盖率门禁通过")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
