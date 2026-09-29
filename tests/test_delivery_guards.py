"""交付面守卫的回归测试。

覆盖的不变量：

* 墙钟停止条件的局部搜索会给出**可执行的复现配方**（迭代数 + stderr 提示 + 元数据标志）。
* HTML 报告与三件套共用 `--force` 覆盖守卫。
* 报告生成失败除 stdout 警告外，还进 ``Result.warnings`` 与 ``run_metadata``。
* 适配器签名过滤丢弃公共参数时写进诊断，不静默。

图形前端（文件对话框、i18n ``_`` 解包、offscreen 下对话框可调用，以及
``--incremental`` 默认值与图形参数的一致性校验）位于独立仓库 MaxTiC-Next-Studio，
其回归测试是 Studio 的 ``tests/test_gui_dialogs.py`` 与 ``tests/test_params_validation.py``。
"""

from __future__ import annotations

import io
import os
import re
import subprocess
import sys
from contextlib import redirect_stderr
from pathlib import Path

import pytest

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
SRC = ROOT / "src" / "maxtic_next"
if str(SRC.parent) not in sys.path:
    sys.path.insert(0, str(SRC.parent))

from maxtic_next import api  # noqa: E402

TREE = str(ROOT / "examples" / "minitree.tree")
CONS = str(ROOT / "examples" / "Cyano_CUTConstraints.tsv")


# --------------------------------------------------------------------------- #
# 墙钟停止条件的可复现配方
# --------------------------------------------------------------------------- #
def _rank(tmp_path, **kw):
    kw.setdefault("html_report", False)
    prefix = str(
        tmp_path / ("p" + re.sub(r"\W", "", str(sorted(kw.items(), key=lambda kv: str(kv)))))
    )
    return api.rank(TREE, [CONS], output_prefix=prefix, print_summary=False, **kw)


def test_wall_clock_local_search_reports_reproduction_recipe(tmp_path) -> None:
    err = io.StringIO()
    with redirect_stderr(err):
        r = _rank(tmp_path, seed=42, local_search=0.2, temperature=0.001)
    meta = r.run_metadata
    joined = " ".join(r.warnings or [])
    assert meta["local_search_deterministic"] is False
    assert meta.get("local_search_iterations") is not None
    assert "--local-search-max-iters" in joined, joined
    assert "--local-search-max-iters" in err.getvalue()


def test_iteration_bounded_local_search_is_marked_deterministic(tmp_path) -> None:
    r = _rank(
        tmp_path, seed=42, local_search=30.0, temperature=0.001, local_search_max_iterations=500
    )
    assert r.run_metadata["local_search_deterministic"] is True
    assert r.run_metadata["local_search_iterations"] == 500
    assert not [w for w in (r.warnings or []) if "--local-search-max-iters" in w and "WARNING" in w]


# --------------------------------------------------------------------------- #
# / HTML 报告的覆盖守卫与失败可见性
# --------------------------------------------------------------------------- #
def test_write_html_report_honours_force(tmp_path) -> None:
    """单元级：既有报告需要 ``force=True`` 才覆盖，且写入是原子的。"""
    from maxtic_next.report.html import write_html_report

    r = api.rank(
        TREE,
        [CONS],
        seed=42,
        print_summary=False,
        html_report=False,
        output_prefix=str(tmp_path / "u"),
        force=True,
    )
    marker = "<!-- 用户自己的文件，不许被静默覆盖 -->"
    (tmp_path / "u.html").write_text("占位" + marker, encoding="utf-8")
    with pytest.raises(FileExistsError):
        write_html_report(str(tmp_path / "u"), r)
    assert marker in (tmp_path / "u.html").read_text(encoding="utf-8")
    path = write_html_report(str(tmp_path / "u"), r, force=True)
    assert marker not in Path(path).read_text(encoding="utf-8")


def test_only_report_present_degrades_without_failing_run(tmp_path) -> None:
    """集成级：只剩报告是"别人的文件"时，三件套照常交付、报告降级并留痕。

    注意"三件套已存在 + 无 force"整条命令会以退出码 3 中止（仓库既有语义，由
    tests/test_correctness_guards.py 的  用例看守），那种情况下根本走不到报告；
    所以这里刻意只保留既有报告、删掉三件套。子进程运行是为了绕开
    ``_WRITTEN_BY_THIS_PROCESS`` 的"同进程自我覆盖允许"设计。
    """
    prefix = str(tmp_path / "half")
    api.rank(
        TREE,
        [CONS],
        seed=42,
        print_summary=False,
        html_report=True,
        output_prefix=prefix,
        force=True,
    )
    marker = "<!-- 别人的报告 -->"
    with open(prefix + ".html", "a", encoding="utf-8") as fh:
        fh.write(marker)
    for suffix in (".mt.informative.tsv", ".mt.conflicts.tsv", ".mt.partial_order.tsv"):
        os.remove(prefix + suffix)

    child = """
import sys, os
sys.path.insert(0, {src!r})
from maxtic_next import api
r = api.rank({tree!r}, [{cons!r}], seed=42, print_summary=False, html_report=True,
             output_prefix={prefix!r}, force=False)
print("REPORT_WRITTEN", bool(r.html_report_file), flush=True)
print("WARN", any("FileExistsError" in w for w in (r.warnings or [])), flush=True)
print("META", "FileExistsError" in str(r.run_metadata.get("html_report_error")), flush=True)
print("TSVOK", os.path.isfile(r.informative_file), flush=True)
"""
    out = subprocess.run(
        [
            sys.executable,
            "-c",
            child.format(src=str(SRC.parent), tree=TREE, cons=CONS, prefix=prefix),
        ],
        capture_output=True,
        text=True,
    )
    log = out.stdout + out.stderr
    assert out.returncode == 0, log
    assert "REPORT_WRITTEN False" in log, log
    assert "WARN True" in log, log
    assert "META True" in log, log
    assert "TSVOK True" in log, log
    assert marker in Path(prefix + ".html").read_text(encoding="utf-8")


def test_html_report_failure_is_recorded_not_only_printed(tmp_path, monkeypatch) -> None:
    """报告失败仍不阻断主流程（设计意图），但必须留下机器可读的痕迹。"""
    import maxtic_next.report.html as html_mod

    def boom(*a, **k):
        raise RuntimeError("模拟模板渲染崩溃")

    monkeypatch.setattr(html_mod, "render_html_report", boom)
    out = io.StringIO()
    with redirect_stderr(out):
        r = api.rank(
            TREE,
            [CONS],
            seed=42,
            print_summary=True,
            html_report=True,
            output_prefix=str(tmp_path / "fail"),
            force=True,
        )
    assert r.html_report_file in (None, "")
    assert any("HTML 报告生成失败" in w for w in (r.warnings or [])), r.warnings
    assert "RuntimeError" in str(r.run_metadata.get("html_report_error"))
    assert os.path.isfile(r.informative_file)  # 主产物不受影响


# --------------------------------------------------------------------------- #
# 适配器签名过滤不再静默
# --------------------------------------------------------------------------- #
def test_public_kwargs_reports_dropped_common_params() -> None:
    from maxtic_next.constraints.adapters.registry import _public_kwargs

    def adapter(min_support=0.0, quiet=False):
        return min_support, quiet

    sink: list[str] = []
    got = _public_kwargs(
        adapter,
        {"min_support": 0.5, "min_endpoint_hit_rate": 0.9, "quiet": True, "unrelated": 1},
        sink,
    )
    assert got == {"min_support": 0.5, "quiet": True}
    assert sink == ["min_endpoint_hit_rate"], sink


def test_note_dropped_writes_diagnostics() -> None:
    from maxtic_next.constraints.constraint import ConstraintSet
    from maxtic_next.constraints.adapters.registry import _note_dropped

    cset = ConstraintSet([])
    _note_dropped(cset, "ale", ["min_endpoint_hit_rate"])
    note = cset.diagnostics.get("adapter_params_not_applicable")
    assert note and "ale" in note and "min_endpoint_hit_rate" in note
    # 没有丢弃时不得凭空写键
    clean = ConstraintSet([])
    _note_dropped(clean, "ale", [])
    assert "adapter_params_not_applicable" not in clean.diagnostics


# --------------------------------------------------------------------------- #
# `--incremental` 默认开启（逃生门 `--no-incremental` 必须仍然可用）
# --------------------------------------------------------------------------- #
def _parse(argv):
    from maxtic_next.cli import build_parser

    return build_parser().parse_args(argv)


def test_incremental_is_now_the_default() -> None:
    base = [TREE, CONS]
    assert _parse(base).incremental is True
    assert _parse(base + ["--local-search", "1"]).incremental is True
    assert _parse(base + ["--incremental"]).incremental is True  # 兼容旧命令行
    assert _parse(base + ["--no-incremental"]).incremental is False  # 逃生门
    assert _parse(base + ["--incremental", "--no-incremental"]).incremental is False


def test_local_search_reports_speed_ratio_metadata(tmp_path) -> None:
    """同一输入下增量与全量的目标值必须一致。"""
    kw = dict(
        seed=42,
        local_search=30.0,
        temperature=0.001,
        local_search_max_iterations=4000,
        print_summary=False,
        html_report=False,
    )
    a = api.rank(TREE, [CONS], output_prefix=str(tmp_path / "full"), incremental=False, **kw)
    b = api.rank(TREE, [CONS], output_prefix=str(tmp_path / "inc"), incremental=True, **kw)
    assert a.best_order == b.best_order
    assert a.values == b.values
    assert Path(a.informative_file).read_text(encoding="utf-8") == Path(
        b.informative_file
    ).read_text(encoding="utf-8")
