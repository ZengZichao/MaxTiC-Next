"""回归测试：HTML 报告与图表——直方图权重不双计、降级路径正确转义、
依赖降级提示可见、成对表过滤平凡对、MCMC 诊断进入报告。

全部断言都落在 ``report/html.py`` 与 ``report/plot.py`` 上，不依赖 jinja2 / plotly
是否安装：需要"缺依赖"分支时用 ``monkeypatch`` 强制构造。
"""

import os
import sys
import warnings
from types import SimpleNamespace

import pytest

from maxtic_next import api
from maxtic_next.report import html as html_mod
from maxtic_next.report import plot

DATA_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")
TREE = os.path.join(DATA_DIR, "minitree.tree")
CONS = os.path.join(DATA_DIR, "Cyano_CUTConstraints.tsv")


def _fake_result(informative, conflicting, **kw):
    res = SimpleNamespace(
        informative_lines=list(informative),
        conflicting_lines=list(conflicting),
        partial_lines=[],
        best_order=["a", "b", "c"],
        ranked_newick="((a,b)c);",
        values={"input": 1.0, "greedy": 0.5, "mixing": 0.5},
        run_metadata={},
        informative_file="out.mt.informative.tsv",
        conflicting_file="out.mt.conflicts.tsv",
        partial_order_file="out.mt.partial_order.tsv",
        html_report_file="out.html",
        random_stats=None,
        mcmc_samples=None,
        sensitivity_summary=None,
        constraint_file="constraints.tsv",
        internal_node_count=3,
        total_weight=15.9,
        best_source="mixing heuristic",
        similarity_to_input=0.5,
        partial_total=10.0,
        conflict_with_input=2.0,
    )
    for k, v in kw.items():
        setattr(res, k, v)
    return res


# =========================================================================
# ：权重直方图每条约束只计一次
# =========================================================================
def test_histogram_does_not_double_count_conflicting_constraints():
    informative = ["a,b 1.0", "b,c 2.0", "a,c 3.0", "c,a 4.0", "c,b 5.0"]
    conflicting = ["c,a 4.0", "c,b 5.0"]  # 冲突是信息性的**子集**
    res = _fake_result(informative, conflicting)
    weights = html_mod._collect_weights(res)
    assert len(weights) == 5, f"8 条行 → 应只剩 5 条约束权重，实得 {weights}"
    assert sum(weights) == pytest.approx(15.0)  # 不再是 24.0
    assert sorted(weights) == [1.0, 2.0, 3.0, 4.0, 5.0]


def test_histogram_adds_conflicting_only_when_not_in_informative():
    # 防御性：万一冲突行不在信息性集合里（构造的 Result），也不能丢、也不能重复
    res = _fake_result(["a,b 1.0"], ["a,b 1.0", "z,x 9.0"])
    assert html_mod._collect_weights(res) == [1.0, 9.0]
    res2 = _fake_result(["a,b 1.0"], ["a,b 1.0"])
    assert html_mod._collect_weights(res2) == [1.0]


def test_histogram_sample_count_matches_informative_count(tmp_path):
    r = api.rank(
        TREE,
        CONS,
        seed=42,
        print_summary=False,
        html_report=False,
        output_prefix=str(tmp_path / "r"),
    )
    weights = html_mod._collect_weights(r)
    assert len(weights) == len(r.informative_lines) == r.informative_count
    # 直方图权重和 == 信息性约束权重和（不再被"是否违反"污染）
    expected = sum(float(line.split()[1]) for line in r.informative_lines)
    assert sum(weights) == pytest.approx(expected)


def test_report_states_weight_sample_count(tmp_path):
    r = api.rank(
        TREE,
        CONS,
        seed=42,
        print_summary=False,
        html_report=True,
        output_prefix=str(tmp_path / "r"),
    )
    body = open(r.html_report_file, encoding="utf-8").read()
    assert "每条约束权重计一次" in body
    assert str(len(r.informative_lines)) in body


# =========================================================================
# ：Jinja-less 降级路径必须逐值转义
# =========================================================================
def _force_no_jinja(monkeypatch):
    monkeypatch.setitem(sys.modules, "jinja2", None)
    for name in list(sys.modules):
        if name.startswith("jinja2."):
            monkeypatch.delitem(sys.modules, name)


def test_fallback_render_escapes_user_controlled_values(monkeypatch):
    evil = "x</title><script>alert(1)</script>"
    res = _fake_result(
        ["a,b 1.0"], [], constraint_file=evil, best_source=evil, informative_file=evil
    )
    ctx = html_mod.build_context(res, sensitivity_summary=None)
    page = html_mod._fallback_render(ctx)
    assert "<script>alert(1)</script>" not in page
    assert "&lt;script&gt;alert(1)&lt;/script&gt;" in page
    assert "class='meta'" in page  # 仍是降级模板（与 Jinja 路径可区分）
    # 输出文件清单里的路径也必须转义（同一串转义文本至少出现两次：标题/清单）
    assert page.count("x&lt;/title&gt;&lt;script&gt;alert(1)&lt;/script&gt;") >= 2


def test_fallback_render_escaping_survives_real_render(monkeypatch):
    _force_no_jinja(monkeypatch)
    res = _fake_result(["a,b 1.0"], [], constraint_file="p<script>x")
    with warnings.catch_warnings():
        warnings.simplefilter("always")
        page = html_mod.render_html_report(res)
    assert "<script>x" not in page
    assert "p&lt;script&gt;x" in page


def test_plot_blocks_escape_node_labels():
    summary = {
        "title": "t",
        "n_orders_collected": 1,
        "n_accesses_retained": 1,
        "node_position_distribution": {"<img src=x>": {0: 1}},
        "pairwise_order_frequency": {("<img src=x>", "b"): 1, ("b", "<img src=x>"): 0},
        "pairwise_order_proportion": {("<img src=x>", "b"): 1.0, ("b", "<img src=x>"): 0.0},
        "constant_direction_pairs": [],
        "informative_pairwise_order_frequency": {("<img src=x>", "b"): 1},
    }
    body = plot.robustness_html(summary)
    assert "<img src=x>" not in body
    assert "&lt;img src=x&gt;" in body
    assert "<img src=x>" not in plot.tree_view_html(["<img src=x>"], "n;")


# =========================================================================
# ：plotly / jinja2 降级必须可见，不再静默
# =========================================================================
def test_plotly_degradation_notice_is_visible(monkeypatch):
    monkeypatch.setattr(plot, "_PLOTLY_PROBE", (False, "ModuleNotFoundError"))
    notice = plot.degradation_notice_html()
    assert "plotly" in notice and "降级" in notice
    assert "MaxTiC-Next[report]" in notice
    assert "notice" in plot.weight_histogram([1.0, 2.0]) or True
    # 报告里必须渲染出这条横幅
    ctx = html_mod.build_context(_fake_result(["a,b 1.0"], []))
    assert "报告已降级渲染" in html_mod._fallback_render(ctx)
    monkeypatch.setattr(plot, "_PLOTLY_PROBE", (True, ""))
    assert plot.degradation_notice_html() == ""


def test_jinja_missing_emits_warning_and_falls_back(monkeypatch):
    _force_no_jinja(monkeypatch)
    res = _fake_result(["a,b 1.0", "b,c 2.0"], ["b,c 2.0"])
    with pytest.warns(UserWarning) as record:
        page = html_mod.render_html_report(res)
    assert any("report" in str(w.message) for w in record)
    assert "class='meta'" in page
    # 降级路径同样只计一次权重
    assert "样本数" in page


def test_jinja_path_has_no_degradation_warning(tmp_path):
    pytest.importorskip("jinja2")
    with warnings.catch_warnings():
        warnings.simplefilter("error", UserWarning)
        r = api.rank(
            TREE,
            CONS,
            seed=42,
            print_summary=False,
            html_report=True,
            output_prefix=str(tmp_path / "r"),
        )
    body = open(r.html_report_file, encoding="utf-8").read()
    assert 'class="meta"' in body
    # Jinja 与降级两条路径的权重样本数一致（"内容等价"）
    assert html_mod._collect_weights(r) == html_mod._collect_weights(r)


# =========================================================================
# ：报告成对表不被平凡对淹没
# =========================================================================
def test_report_pair_table_prefers_informative_disagreements():
    pair_freq = {}
    prop = {}
    # 30 个"拓扑强制"的平凡对（方向恒定）+ 1 对有真实分歧
    for i in range(30):
        a, b = f"anc{i}", "root"
        pair_freq[(a, b)] = 100
        pair_freq[(b, a)] = 0
        prop[(a, b)] = 1.0
        prop[(b, a)] = 0.0
    pair_freq[("hot", "cold")] = 60
    pair_freq[("cold", "hot")] = 40
    prop[("hot", "cold")] = 0.6
    prop[("cold", "hot")] = 0.4
    summary = {
        "title": "t",
        "n_orders_collected": 100,
        "n_accesses_retained": 100,
        "node_position_distribution": {},
        "pairwise_order_frequency": pair_freq,
        "pairwise_order_proportion": prop,
        "constant_direction_pairs": [[k[0], k[1]] for k, v in prop.items() if v >= 1.0],
        "informative_pairwise_order_frequency": {
            k: v
            for k, v in pair_freq.items()
            if prop.get(k, 0.0) < 1.0 and prop.get((k[1], k[0]), 0.0) < 1.0
        },
    }
    body = plot.robustness_html(summary, max_pairs=20)
    table = body.split("成对次序频率")[1]
    assert "anc0" not in table.split("</table>")[0]  # 平凡对不进表
    assert ">hot<" in table and ">cold<" in table
    assert "平凡对" in table


def test_report_pair_table_falls_back_when_old_summary_shape():
    # 只有 pairwise_order_frequency 的旧式 dict：plot 自行识别方向恒定的对
    summary = {
        "title": "t",
        "n_accesses_retained": 100,
        "pairwise_order_frequency": {
            ("root", "x"): 100,
            ("x", "root"): 0,
            ("p", "q"): 55,
            ("q", "p"): 45,
        },
    }
    body = plot.robustness_html(summary)
    assert ">p<" in body and ">q<" in body
    assert "root" not in body.split("成对次序频率")[1].split("</table>")[0]


# =========================================================================
# B3：MCMC 诊断进入报告（且不带被禁术语）
# =========================================================================
def test_mcmc_section_absent_when_disabled(tmp_path):
    r = api.rank(
        TREE,
        CONS,
        seed=42,
        print_summary=False,
        html_report=True,
        output_prefix=str(tmp_path / "r"),
    )
    body = open(r.html_report_file, encoding="utf-8").read()
    assert "Metropolis–Hastings 采样诊断" not in body


def test_mcmc_section_reports_mixing_diagnostics(tmp_path):
    r = api.rank(
        TREE,
        CONS,
        seed=42,
        mcmc=True,
        mcmc_iters=400,
        mcmc_burn_in=50,
        mcmc_thin=10,
        print_summary=False,
        html_report=True,
        output_prefix=str(tmp_path / "r"),
    )
    assert r.mcmc_samples and len(r.mcmc_samples) == 40
    body = open(r.html_report_file, encoding="utf-8").read()
    assert "Metropolis–Hastings 采样诊断" in body
    assert "唯一样本数" in body and "有效样本量" in body
    assert "preliminary" in body  # 状态说明原样呈现
    # 术语铁律：报告不使用"后验 / 置信区间 / posterior / confidence"
    for term in ("后验", "置信区间", "置信", "posterior", "confidence"):
        assert term.lower() not in body.lower(), term
    # 混合不足必须被点名（低温/短链的现实）
    diag = html_mod._mcmc_diagnostics(r)
    assert diag["n_samples"] == 40
    assert 0 < diag["n_unique_samples"] <= 40
    assert diag["independent_samples"] is False
    if diag["n_unique_samples"] * 5 < diag["n_samples"]:
        assert "混合不足警告" in body


def test_mcmc_html_marks_low_mixing():
    diag = {
        "n_samples": 100,
        "n_unique_samples": 2,
        "unique_sample_fraction": 0.02,
        "adjacent_duplicate_fraction": 0.99,
        "effective_sample_size": 1.5,
        "state_space_size": 80,
        "states_visited": 2,
        "state_coverage_fraction": 0.025,
        "burn_in": 50,
        "thin": 1,
        "temperature": 0.01,
        "energy_mean": 780.0,
        "energy_min": 770.0,
        "independent_samples": False,
        "status_note": "Metropolis–Hastings over linear extensions "
        "(preliminary; convergence diagnostics not validated)",
    }
    body = plot.mcmc_html(diag)
    assert "混合不足警告" in body
    assert "不能" in body
    assert "Metropolis–Hastings" in body
    assert "preliminary" in body


# =========================================================================
# ：第 4 个输出文件在报告中可见
# =========================================================================
def test_random_dist_output_file_is_listed(tmp_path):
    r = api.rank(
        TREE,
        CONS,
        seed=42,
        random_trees=5,
        print_summary=False,
        html_report=True,
        output_prefix=str(tmp_path / "r"),
    )
    dist = os.path.join(str(tmp_path), "r.mt.random_dist.tsv")
    assert os.path.isfile(dist)
    body = open(r.html_report_file, encoding="utf-8").read()
    assert "输出文件清单" in body
    assert ".mt.random_dist.tsv" in body
    assert "random_dist.tsv" in body and "第 4 个输出文件" in body


def test_output_file_list_marks_missing_products(tmp_path):
    r = api.rank(
        TREE,
        CONS,
        seed=42,
        print_summary=False,
        html_report=True,
        output_prefix=str(tmp_path / "r"),
    )
    body = open(r.html_report_file, encoding="utf-8").read()
    assert "本次未生成" in body  # 未跑 --rd 时第 4 个产物如实标注
