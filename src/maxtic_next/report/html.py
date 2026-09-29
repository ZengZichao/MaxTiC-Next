"""交互式 HTML 报告构建。

优先使用 Jinja2 渲染 ``report/templates/report.html.j2``（含 partials）；若 Jinja2 未
安装（离线环境），降级为内置字符串模板渲染（等价内容、内联 SVG 图表），保证报告始终
可生成、可双击打开、无后端依赖。

降级是**有意为之的能力保留**，但不再静默：当报告被请求（默认开启）而
``jinja2`` / ``plotly`` 缺失时，本模块会通过 ``warnings.warn`` 显式提示
用户执行 ``pip install "MaxTiC-Next[report]"`` 以获得完整交互式报告，然后再降级渲染。
这样既保证离线可用，又让缺失依赖对用户可见（而非悄悄退化到功能较弱的模板）。

报告整合：
* 运行元数据（种子、参数、版本）；
* 冲突比例；
* 权重分布直方图（``plot.weight_histogram``）——**每条约束权重只计一次**（
  ``conflicting_lines`` 是 ``informative_lines`` 的子集，拼接求和会双计被违反约束）；
* 基于局部搜索访问解的稳健性/敏感性摘要（``robustness.sensitivity.NearOptimalCollector.summary``）；
* Metropolis–Hastings 链的混合诊断（``robustness.mcmc.mcmc_diagnostics``）；
* 输出文件清单（含第 4 个产物 ``*.mt.random_dist.tsv``）；
* 排序树可视化（``plot.tree_view_html``）。

两条渲染路径**都必须转义**用户可控内容：Jinja 路径靠 ``autoescape``，
降级路径靠 ``_e()``；依赖缺失（jinja2 / plotly）时都渲染可见的降级横幅。

术语铁律：报告全文仅使用"稳健性/敏感性摘要"，不出现"置信区间 /
后验概率 / 后验分布"等用语。
"""

import os
import warnings
from html import escape as _html_escape
from typing import Dict, List, Optional, Tuple

from maxtic_next.io.output import write_text_file_atomic
from maxtic_next.report import plot

_TEMPLATE_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "templates")


def _e(s) -> str:
    """降级模板用的 HTML 转义：仅 ``None`` 视为空，其余一律转义。"""
    if s is None:
        return ""
    return _html_escape(str(s))


def _parse_constraint_lines(lines) -> Tuple[Dict[str, float], List[str]]:
    """把 ``"donor,receptor weight"`` 行解析为 ``(edge, edge_keys)``（保序、去重）。"""
    edge: Dict[str, float] = {}
    keys: List[str] = []
    for line in lines or []:
        parts = str(line).split()
        if len(parts) < 2:
            continue
        key = parts[0]
        try:
            weight = float(parts[1])
        except ValueError:
            continue
        if key not in edge:
            edge[key] = weight
            keys.append(key)
    return edge, keys


def _collect_weights(result) -> List[float]:
    """从结果的信息性约束行抽取权重，用于权重分布直方图。

    ：``conflicting_lines`` 是 ``informative_lines`` 的**子集**（同一
    ``"key weight"`` 格式），旧实现把两个列表拼接取权重，导致**被违反的约束计两次**
    （实测 8 条约束 → 11 个权重样本，和 19.7 vs 真实 15.9，分布形状还被"是否违反"
    这一结果反向污染）。现按边键去重：每条约束的权重只计一次，且仅在它未出现过时
    才从 conflicting 补入（防御两侧不一致的构造结果）。
    """
    weights: List[float] = []
    seen = set()
    for source in (
        getattr(result, "informative_lines", []) or [],
        getattr(result, "conflicting_lines", []) or [],
    ):
        for line in source:
            parts = str(line).split()
            if len(parts) < 2:
                continue
            key = parts[0]
            if key in seen:
                continue
            try:
                weights.append(float(parts[1]))
            except ValueError:
                continue
            seen.add(key)
    return weights


def _mcmc_diagnostics(result) -> Optional[Dict]:
    """由 ``Result.mcmc_samples`` 组装 MH 链诊断（诊断必须随报告展示）。

    能量用交付口径重算：以 ``informative_lines`` 构造 ``(edge, edge_keys)`` 后调用
    ``ranking.value.value``，即"该序被违反的信息性约束权重和"（与 MaxTiC 目标函数
    在被违反集合上一致）。若上层（``Ranker``）已在 ``run_metadata`` 写入权威诊断，
    则以其为准并补齐缺失键。
    """
    samples = getattr(result, "mcmc_samples", None)
    meta = getattr(result, "run_metadata", {}) or {}
    if not samples:
        if not meta.get("mcmc"):
            return None
        return dict(meta.get("mcmc_diagnostics") or {}) or None
    from maxtic_next.robustness.mcmc import mcmc_diagnostics
    from maxtic_next.ranking.value import value as _value

    edge, keys = _parse_constraint_lines(getattr(result, "informative_lines", []))
    energies = [_value(list(s), edge, keys) for s in samples] if keys else None
    diag = mcmc_diagnostics(samples, energies)
    injected = meta.get("mcmc_diagnostics")
    if isinstance(injected, dict):
        for k, v in injected.items():
            diag[k] = v
    for k, label in (
        ("mcmc_burn_in", "burn_in"),
        ("mcmc_thin", "thin"),
        ("mcmc_temperature", "temperature"),
        ("mcmc_status", "status_note"),
    ):
        if k in meta and diag.get(label) in (None, ""):
            diag[label] = meta[k]
    diag.setdefault("n_samples", len(samples))
    diag.setdefault(
        "status_note",
        "Metropolis–Hastings over linear extensions "
        "(preliminary; convergence diagnostics not "
        "validated)",
    )
    return diag


def _output_files_html(result) -> str:
    """产物文件清单（含第 4 个输出文件 ``*.mt.random_dist.tsv``）。"""
    from maxtic_next.config import OUTPUT_SUFFIXES

    files: List[Tuple[str, str]] = [
        ("信息性约束（informative）", getattr(result, "informative_file", "")),
        ("与最优序冲突的约束（conflicts）", getattr(result, "conflicting_file", "")),
        ("偏序产物（partial order）", getattr(result, "partial_order_file", "")),
    ]
    random_stats = getattr(result, "random_stats", None) or {}
    files.append(
        (
            "随机序目标值分布（random trees，第 4 个输出文件；"
            f"short 风格后缀 {OUTPUT_SUFFIXES['short'][3]}，"
            f"legacy 风格后缀 {OUTPUT_SUFFIXES['legacy'][3]}；"
            "仅 --random-trees > 0 时生成）",
            random_stats.get("distribution_file", "") if random_stats else "",
        )
    )
    files.append(("交互式 HTML 报告", getattr(result, "html_report_file", "")))
    rows = []
    for label, path in files:
        rows.append(
            f"<tr><td>{_e(label)}</td>"
            f"<td><code>{_e(path) if path else '（本次未生成）'}</code></td></tr>"
        )
    return (
        "<table class='data'><thead><tr><th>产物</th><th>路径</th></tr></thead>"
        f"<tbody>{''.join(rows)}</tbody></table>"
    )


def _conflict_ratio(result) -> float:
    """冲突比例：与输入树冲突的偏序权重 / 偏序总权重。"""
    partial_total = getattr(result, "partial_total", 0.0) or 0.0
    conflict_with_input = getattr(result, "conflict_with_input", 0.0) or 0.0
    if partial_total == 0.0:
        return 0.0
    return conflict_with_input / partial_total


def build_context(
    result, sensitivity_summary: Optional[Dict] = None, weights: Optional[List[float]] = None
) -> Dict:
    """组装渲染上下文（Jinja2 与降级渲染共用）。"""
    from maxtic_next import __version__

    if weights is None:
        weights = _collect_weights(result)
    meta = getattr(result, "run_metadata", {}) or {}
    mcmc_diag = _mcmc_diagnostics(result)
    return {
        "version": __version__,
        "result": result,
        "run_metadata": meta,
        "seed": meta.get("seed", "?"),
        "params": meta.get("params", ""),
        "constraint_file": getattr(result, "constraint_file", ""),
        "internal_node_count": getattr(result, "internal_node_count", 0),
        "total_weight": getattr(result, "total_weight", 0.0),
        "best_source": getattr(result, "best_source", ""),
        "similarity_to_input": getattr(result, "similarity_to_input", 0.0),
        "conflict_ratio": _conflict_ratio(result),
        "values": getattr(result, "values", {}),
        "sensitivity_summary": sensitivity_summary,
        "weights": weights,
        "weight_histogram_html": plot.weight_histogram(weights),
        "robustness_html": plot.robustness_html(sensitivity_summary),
        "mcmc_diag": mcmc_diag,
        "mcmc_html": plot.mcmc_html(mcmc_diag),
        "degradation_html": plot.degradation_notice_html(),
        "output_files_html": _output_files_html(result),
        "n_weight_samples": len(weights),
        "tree_view_html": plot.tree_view_html(
            getattr(result, "best_order", []), getattr(result, "ranked_newick", "")
        ),
    }


def render_html_report(
    result, sensitivity_summary: Optional[Dict] = None, weights: Optional[List[float]] = None
) -> str:
    """渲染 HTML 报告字符串。

    优先 Jinja2；缺失则降级为内置字符串模板（内容等价）。
    """
    ctx = build_context(result, sensitivity_summary=sensitivity_summary, weights=weights)
    try:
        from jinja2 import Environment, FileSystemLoader, select_autoescape

        env = Environment(
            loader=FileSystemLoader(_TEMPLATE_DIR),
            autoescape=select_autoescape(default_for_string=True, default=True),
        )
        template = env.get_template("report.html.j2")
        return template.render(**ctx)
    except ImportError as exc:
        # 依赖缺失：显式告警后降级（不静默）。提示用户补齐 report 可选依赖。
        warnings.warn(
            f"HTML 报告所需的依赖缺失（{type(exc).__name__}: {exc}），已降级为内置模板。"
            f"如需完整交互式报告，请安装可选依赖："
            f'pip install "MaxTiC-Next[report]"（提供 jinja2 / plotly）。',
            stacklevel=2,
        )
        return _fallback_render(ctx)
    except Exception as exc:  # noqa: BLE001 - 模板加载/渲染错误，降级
        warnings.warn(
            f"Jinja2 模板渲染失败（{type(exc).__name__}: {exc}），已降级为内置模板。"
            f"如需完整交互式报告，请确认模板文件完整或重装 "
            f'pip install "MaxTiC-Next[report]"。',
            stacklevel=2,
        )
        return _fallback_render(ctx)


def _fallback_render(ctx: Dict) -> str:
    """Jinja2 缺失时的降级 HTML 渲染（内联样式，自包含）。

    ：Jinja 路径有 ``autoescape``，降级路径**必须**逐值转义，否则用户可控的
    约束文件路径 / 类群名可注入原始标签。这里除由 ``plot`` 生成并已完成转义的图表
    片段（``weight_histogram_html`` / ``robustness_html`` / ``mcmc_html`` /
    ``tree_view_html``）之外，所有插值一律经 :func:`_e` 转义。
    """
    style = """
    <style>
      body { font-family: -apple-system, 'Segoe UI', Arial, sans-serif; margin: 24px; color: #222; }
      h1 { border-bottom: 2px solid #4C72B0; padding-bottom: 6px; }
      h3, h4 { color: #4C72B0; }
      table.data { border-collapse: collapse; margin: 8px 0; font-size: 13px; }
      table.data th, table.data td { border: 1px solid #ccc; padding: 4px 8px; text-align: left; }
      .muted { color: #888; }
      .notice { background: #fff6e5; border: 1px solid #e6b800; padding: 8px 12px;
                border-radius: 6px; }
      .meta { background: #f5f7fa; padding: 10px 14px; border-radius: 6px; }
      pre.newick { background: #f5f7fa; padding: 8px; overflow-x: auto; }
      .card { border: 1px solid #e0e0e0; border-radius: 8px; padding: 12px 16px; margin: 12px 0; }
    </style>
    """
    meta_rows = [
        ("MaxTiC-Next 版本", ctx["version"]),
        ("种子", ctx["seed"]),
        ("参数", ctx["params"]),
        ("约束文件", ctx["constraint_file"]),
        ("内部节点数", ctx["internal_node_count"]),
        ("约束总权重", ctx["total_weight"]),
        ("最优序来源", ctx["best_source"]),
        ("与输入树相似度", ctx["similarity_to_input"]),
        ("冲突比例", f"{ctx['conflict_ratio']:.4f}"),
    ]
    meta_html = (
        "<div class='meta'>"
        + "".join(f"<b>{_e(k)}</b>：{_e(v)}<br/>" for k, v in meta_rows)
        + "</div>"
    )
    vals = ctx["values"]
    values_html = "<ul>"
    for k in ("input", "greedy", "mixing", "local_search", "best"):
        if k in vals:
            values_html += f"<li>{_e(k)}：{_e(vals[k])}</li>"
    values_html += "</ul>"

    mcmc_card = ""
    if ctx.get("mcmc_diag"):
        mcmc_card = f"<div class='card'>{ctx['mcmc_html']}</div>"
    degradation = ctx.get("degradation_html") or ""

    return f"""<!DOCTYPE html>
<html lang="zh-CN">
<head><meta charset="utf-8"/>
<meta name="viewport" content="width=device-width, initial-scale=1"/>
<title>MaxTiC-Next 报告 - {_e(ctx["constraint_file"])}</title>{style}</head>
<body>
<h1>MaxTiC-Next 排序报告</h1>
{degradation}
<div class='card'><h3>运行元数据</h3>{meta_html}</div>
<div class='card'><h3>目标值</h3>{values_html}</div>
<div class='card'><h3>约束权重分布</h3>{ctx["weight_histogram_html"]}
<p class='muted'>样本数（每条约束权重计一次，被违反的约束不重复计入）：
{_e(ctx.get("n_weight_samples"))}</p></div>
<div class='card'>{ctx["robustness_html"]}</div>
{mcmc_card}
<div class='card'><h3>排序树可视化</h3>{ctx["tree_view_html"]}</div>
<div class='card'><h3>输出文件清单</h3>{ctx["output_files_html"]}</div>
<hr/>
<p class='muted'>本报告为静态文件，无后端依赖，可双击打开。稳健性/敏感性摘要仅基于
局部搜索过程中访问到的解，反映各节点次序的稳健程度，不涉及概率推断。</p>
</body>
</html>"""


def write_html_report(
    prefix: str,
    result,
    sensitivity_summary: Optional[Dict] = None,
    weights: Optional[List[float]] = None,
    force: bool = False,
) -> str:
    """渲染并写出 HTML 报告，返回文件路径（``<prefix>.html``）。

    ``force`` 与三件套同语义：既有报告不会被静默覆盖。
    """
    html = render_html_report(result, sensitivity_summary=sensitivity_summary, weights=weights)
    path = prefix + ".html"
    return write_text_file_atomic(path, html, force=force)
