"""术语铁律守卫测试：MCMC 关闭（默认）时，排序报告与稳健性/敏感性摘要
不得出现「后验 / 置信区间 / posterior / confidence / 置信」等 MCMC-only 用语。

这些用语仅允许出现在 ``robustness/mcmc.py``。
局部搜索产出的摘要严格称为「基于局部搜索访问解的稳健性/敏感性摘要」，不涉及概率推断。
"""

import os

from maxtic_next.api import rank
from maxtic_next.io.output import format_summary

DATA_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")
TREE = os.path.join(DATA_DIR, "minitree.tree")
CONS = os.path.join(DATA_DIR, "Cyano_CUTConstraints.tsv")

# MCMC-only 用语（局部搜索 / 报告路径禁止出现）
_FORBIDDEN = ["后验", "置信区间", "置信", "posterior", "confidence"]


def _assert_no_forbidden(text: str, where: str) -> None:
    low = text.lower()
    for term in _FORBIDDEN:
        assert term.lower() not in low, f"{where} 出现禁用 MCMC 用语：{term!r}"


def test_summary_no_mcmc_terms():
    """stdout 摘要（MCMC 关）不含禁用 MCMC 用语。"""
    r = rank(TREE, CONS, seed=42, print_summary=False)
    out = format_summary(r)
    _assert_no_forbidden(out, "format_summary")


def test_sensitivity_summary_no_mcmc_terms(tmp_path):
    """局部搜索产生的稳健性/敏感性摘要（MCMC 关）不含禁用 MCMC 用语。"""
    prefix = str(tmp_path / "r")
    r = rank(TREE, CONS, seed=42, local_search=0.02, output_prefix=prefix, print_summary=False)
    assert r.sensitivity_summary is not None
    # 摘要标题应仅为「稳健性/敏感性摘要」（术语铁律）
    assert "稳健性/敏感性" in r.sensitivity_summary["title"]
    _assert_no_forbidden(str(r.sensitivity_summary), "sensitivity_summary")


def test_html_report_no_mcmc_terms(tmp_path):
    """交互式 HTML 报告（MCMC 关）不含禁用 MCMC 用语。"""
    prefix = str(tmp_path / "r")
    r = rank(
        TREE,
        CONS,
        seed=42,
        local_search=0.02,
        output_prefix=prefix,
        html_report=True,
        print_summary=False,
    )
    assert r.html_report_file
    with open(r.html_report_file, "r", encoding="utf-8") as fh:
        html = fh.read()
    _assert_no_forbidden(html, "html_report")
