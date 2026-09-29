"""Jinja 渲染回归护栏。

目标：在 jinja2 完备环境下，HTML 报告必须真正走 Jinja2 渲染路径，绝不触发
``warnings.warn`` 降级（例如 ``robustness.j2`` 那种只在 jinja2 完备时才暴露的模板
bug），并且产出文件必须包含 Jinja 模板专属标记（而非降级字符串模板）。

判别依据：
* ``render_html_report`` 在缺依赖/渲染失败时会 ``warnings.warn(UserWarning)``，
  正常 Jinja 路径不产生任何警告 —— 故用 ``warnings.simplefilter("error",
  UserWarning)`` 把任何降级告警变成测试失败。
* Jinja 模板 ``report/templates/partials/summary.j2`` 渲染出 ``<div class="meta">``
  （**双引号**）；而降级字符串模板 ``report/html.py::_fallback_render`` 只用单引号
  ``<div class='meta'>``。``class="meta"``（双引号）因此成为「真 Jinja 渲染」的专属
  标记，降级模板永远不可能产出它。
"""

import os
import warnings

from maxtic_next.api import rank

DATA_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")
TREE = os.path.join(DATA_DIR, "minitree.tree")
CONS = os.path.join(DATA_DIR, "Cyano_CUTConstraints.tsv")

# Jinja 模板专属标记：summary.j2 渲染的 <div class="meta">（双引号）。
_JINJA_META_MARKER = 'class="meta"'

# MCMC-only 用语（局部搜索 / 报告路径禁止出现）
_FORBIDDEN = ["后验", "置信区间", "置信", "posterior", "confidence"]


def _assert_no_forbidden(text: str, where: str) -> None:
    low = text.lower()
    for term in _FORBIDDEN:
        assert term.lower() not in low, f"{where} 出现禁用 MCMC 用语：{term!r}"


def test_html_report_uses_jinja_when_available(tmp_path):
    """jinja2 完备时，HTML 报告必须走 Jinja 渲染、无 UserWarning、含 Jinja 专属标记。

    同时保留「无 MCMC 禁用词」断言（术语铁律）。
    """
    prefix = str(tmp_path / "r")
    # 把任何 UserWarning（含 render_html_report 的降级告警）升级为异常：
    # 一旦走了 fallback，测试立即失败。
    with warnings.catch_warnings():
        warnings.simplefilter("error", UserWarning)
        r = rank(
            TREE,
            CONS,
            seed=42,
            local_search=0.02,
            output_prefix=prefix,
            html_report=True,
            print_summary=False,
        )

    assert r.html_report_file, "未生成 HTML 报告文件"
    with open(r.html_report_file, "r", encoding="utf-8") as fh:
        html = fh.read()

    # 1) Jinja 专属标记必须存在（证明走的是 Jinja 而非降级字符串模板）
    assert _JINJA_META_MARKER in html, (
        f"HTML 报告未包含 Jinja 专属标记 {_JINJA_META_MARKER!r}，"
        f"可能走了降级模板：{r.html_report_file}"
    )
    # 2) 降级模板使用单引号 class='meta'，确认本次确非降级
    assert "class='meta'" not in html, (
        f"HTML 报告含降级模板标记 class='meta'，疑似未走 Jinja 渲染：{r.html_report_file}"
    )
    # 3) 术语铁律：MCMC 关闭（默认）时不得出现 MCMC-only 用语
    _assert_no_forbidden(html, "html_report")
