"""报告图表构建。

优先使用 Plotly 内联 JS 生成交互式图表；若 ``plotly`` 未安装（离线环境），自动降级为
**内联 SVG / 文本表格**，保证报告始终可生成、可双击打开。
降级**不再静默**：:func:`has_plotly` 探测依赖，缺失时报告顶部由
:func:`degradation_notice_html` 渲染一条可见提示，说明退化原因与恢复方法。

提供的图表：
* ``weight_histogram``：约束权重分布直方图；
* ``robustness_html``：基于局部搜索访问解的稳健性/敏感性摘要（节点位置分布 + 成对次序频率）；
* ``mcmc_html``：Metropolis–Hastings 链的混合诊断（必须与统计量一并展示）；
* ``tree_view_html``：排序树可视化（按 rank 的纵向阶梯 + Newick 文本）。

术语铁律：本模块**严禁**出现"置信区间 / 后验概率 / 后验分布"等用语。
"""

from typing import Dict, List, Optional, Tuple
from html import escape as _html_escape

#: 报告中"成对次序频率"表默认展示的行数
PAIR_TABLE_ROWS = 20


# ----------------------------------------------------------------------
# 报告文案（中 / 英）。所有面向读者的标签都从这里取，
# 图表与摘要函数经 ``lang`` 参数选择语言；非法值回退中文。
# 术语铁律同样约束英文：不出现 posterior / confidence 等用语。
# ----------------------------------------------------------------------
_LABELS: Dict[str, Dict[str, str]] = {
    "zh": {
        "report_degraded": "报告已降级渲染",
        "plotly_missing": (
            "plotly 未安装：交互式图表已降级为内联 SVG 静态图（信息不丢失，但不可缩放/悬停）。"
            '如需完整交互式报告，请执行 pip install "MaxTiC-Next[report]"。'
        ),
        "probe_info": "（探测信息：{}）",
        "hist_title": "约束权重分布",
        "hist_x": "约束权重",
        "hist_y": "约束数量",
        "hist_aria": "权重分布直方图",
        "no_weight_data": "无权重数据。",
        "robust_placeholder": (
            "未运行局部搜索（--local-search / --ls &gt; 0），"
            "无基于局部搜索访问解的稳健性/敏感性摘要。</p>"
        ),
        "robust_title": "基于局部搜索访问解的稳健性/敏感性摘要",
        "li_orders_dedup": "去重近优排序数（参与统计）：<b>{}</b>",
        "li_orders_retained": "保留的去重近优排序数：<b>{}</b>",
        "li_accesses_retained": "上述排序的<b>真实访问次数</b>合计：{}",
        "li_accesses_total": "局部搜索累计访问（含被容量裁剪淘汰者）：{}",
        "li_unique_total": "累计见过的不同排序：{}{}",
        "lb_lower_bound": "（下界）",
        "li_top_k": ("保留上限 top_k：{}（本次报告使用 {}）"),
        "li_best_value": "最优目标值：{}",
        "basis_fmt": (
            "口径：统计在“目标值最小的前 {k} 个去重近优排序”上进行，按每个序被局部"
            "搜索访问的真实次数加权；不是对输入数据重抽样的稳健度，也不是概率区间。"
        ),
        "robust_note": (
            "注意：本摘要不是对输入数据重抽样的稳健度，"
            "也不是概率区间；它只描述近优解集合内部各节点次序的敏感性。"
        ),
        "pos_dist_h": "节点位置分布（按 rank 位置的真实访问次数）",
        "th_node": "节点",
        "th_modal_pos": "众数位置",
        "th_pos_span": "位置跨度",
        "th_visits": "访问次数",
        "th_modal_share": "众数位置占比",
        "freq_note": "（占比按真实访问次数计算，非按去重排序数。）",
        "pair_h": "成对次序频率（a 在 b 之前的<b>真实访问次数</b>）",
        "pair_note_pre": (
            "下表已剔除方向恒定为 100% 的平凡对（多为物种树拓扑强制的"
            "祖先-后代对，其方向在<b>所有</b>合法排序中都相同、不含排序信息），"
            "以便看到真正有信息量的近优分歧。被剔除的平凡对共 "
        ),
        "pair_note_post": "对（{} 个有向键）。",
        "th_a": "a",
        "th_b": "b",
        "th_count": "访问次数",
        "th_freq": "频率",
        "pair_empty": "近优集合内所有成对次序方向均恒定（无可见的近优分歧）。",
        "trivial_details": "方向恒定的平凡对（前 12 个，点击展开）",
        "trivial_more": " …（另有 {} 个）",
        "mcmc_placeholder": (
            "未运行 Metropolis–Hastings 采样（--mcmc 未开启），无链混合诊断。</p>"
        ),
        "mcmc_h": "Metropolis–Hastings 采样诊断",
        "mcmc_status": "<b>状态</b>：{}",
        "mcmc_th_diag": "诊断量",
        "mcmc_th_value": "值",
        "lbl_n_samples": "记录样本数",
        "lbl_n_unique": "唯一样本数",
        "lbl_unique_frac": "唯一样本比例",
        "lbl_adj_dup": "相邻重复率",
        "lbl_max_run": "最长同态连续段占比",
        "lbl_ess": "有效样本量 (ESS)",
        "lbl_ess_frac": "ESS / 样本数",
        "lbl_state_space": "状态空间（线性扩展数）",
        "lbl_states_visited": "已访问状态数",
        "lbl_coverage": "状态访问比例",
        "lbl_burn_in": "burn-in 步数",
        "lbl_thin": "抽稀间隔",
        "lbl_temperature": "温度 T",
        "lbl_energy_mean": "能量均值",
        "lbl_energy_min": "能量最小值",
        "lbl_energy_sd": "能量标准差",
        "mix_warn_pre": "<b>混合不足警告</b>：{} 个记录样本中只有 {} 个不同"
        "排序，链实际上没有遍历状态空间；此时的均值/最小值等统计量<b>不能</b>"
        "用作排序不确定性的度量（提高温度、加长链或增大 burn-in 后重跑）。",
        "mcmc_corr_note": (
            "链上样本彼此相关，<b>不是独立样本</b>；任何均值/分位数解释"
            "必须连同有效样本量 (ESS) 与状态访问比例一起阅读。"
        ),
        "tree_no_order": "无排序结果。",
        "tree_title": "排序阶梯",
        "tree_newick_h": "排序树 Newick",
        "of_th_item": "产物",
        "of_th_path": "路径",
        "of_informative": "信息性约束（informative）",
        "of_conflicts": "与最优序冲突的约束（conflicts）",
        "of_partial": "偏序产物（partial order）",
        "of_random_pre": "随机序目标值分布（random trees，第 4 个输出文件；"
        "short 风格后缀 {}，legacy 风格后缀 {}；仅 --random-trees > 0 时生成）",
        "of_html_report": "交互式 HTML 报告",
        "of_not_generated": "（本次未生成）",
    },
    "en": {
        "report_degraded": "Report rendered in degraded mode",
        "plotly_missing": (
            "plotly is not installed: interactive charts are downgraded to inline "
            "static SVG (same information, but no zoom / hover). For the full "
            'interactive report run: pip install "MaxTiC-Next[report]".'
        ),
        "probe_info": " (probe: {})",
        "hist_title": "Constraint weight distribution",
        "hist_x": "Constraint weight",
        "hist_y": "Number of constraints",
        "hist_aria": "Weight distribution histogram",
        "no_weight_data": "No weight data.",
        "robust_placeholder": (
            "Local search was not run (--local-search / --ls &gt; 0), so no "
            "robustness/sensitivity summary based on locally searched solutions "
            "is available.</p>"
        ),
        "robust_title": "Robustness/sensitivity summary based on locally searched solutions",
        "li_orders_dedup": "Deduplicated near-optimal orders (counted): <b>{}</b>",
        "li_orders_retained": "Deduplicated near-optimal orders retained: <b>{}</b>",
        "li_accesses_retained": "Total <b>actual visits</b> of the retained orders: {}",
        "li_accesses_total": ("Total local-search visits (including capacity-evicted ones): {}"),
        "li_unique_total": "Distinct orders ever seen: {}{}",
        "lb_lower_bound": "(lower bound)",
        "li_top_k": ("Retention cap top_k: {} (this report uses {})"),
        "li_best_value": "Best objective value: {}",
        "basis_fmt": (
            "Basis: statistics are computed over the {k} deduplicated near-optimal "
            "orders with the smallest objective values, weighted by how often each "
            "order was actually visited by local search; this is not a resampling "
            "robustness of the input data and not a probabilistic interval."
        ),
        "robust_note": (
            "Note: this summary is not a resampling robustness of the input data "
            "and not a probabilistic interval; it only describes the sensitivity "
            "of node orders within the near-optimal set."
        ),
        "pos_dist_h": "Node position distribution (actual visits per rank position)",
        "th_node": "Node",
        "th_modal_pos": "Modal position",
        "th_pos_span": "Position span",
        "th_visits": "Visits",
        "th_modal_share": "Modal share",
        "freq_note": "(Shares are computed from actual visit counts, not from deduplicated orders.)",
        "pair_h": "Pairwise order frequency (actual visits with a before b)",
        "pair_note_pre": (
            "Pairs whose direction is constant at 100% are excluded (mostly "
            "ancestor-descendant pairs forced by the species-tree topology; their "
            "direction is identical in <b>every</b> legal order and carries no "
            "ranking information), so that genuinely informative near-optimal "
            "disagreements stand out. Excluded trivial pairs: "
        ),
        "pair_note_post": " ({} oriented keys).",
        "th_a": "a",
        "th_b": "b",
        "th_count": "Visits",
        "th_freq": "Frequency",
        "pair_empty": (
            "Every pairwise order direction is constant within the near-optimal "
            "set (no visible near-optimal disagreement)."
        ),
        "trivial_details": "Constant-direction trivial pairs (first 12, click to expand)",
        "trivial_more": " … ({} more)",
        "mcmc_placeholder": (
            "Metropolis–Hastings sampling was not run (--mcmc off); no chain "
            "mixing diagnostics.</p>"
        ),
        "mcmc_h": "Metropolis–Hastings sampling diagnostics",
        "mcmc_status": "<b>Status</b>: {}",
        "mcmc_th_diag": "Diagnostic",
        "mcmc_th_value": "Value",
        "lbl_n_samples": "Recorded samples",
        "lbl_n_unique": "Unique samples",
        "lbl_unique_frac": "Unique sample fraction",
        "lbl_adj_dup": "Adjacent duplicate rate",
        "lbl_max_run": "Longest same-state run fraction",
        "lbl_ess": "Effective sample size (ESS)",
        "lbl_ess_frac": "ESS / samples",
        "lbl_state_space": "State space (linear extensions)",
        "lbl_states_visited": "States visited",
        "lbl_coverage": "State coverage fraction",
        "lbl_burn_in": "Burn-in steps",
        "lbl_thin": "Thinning interval",
        "lbl_temperature": "Temperature T",
        "lbl_energy_mean": "Energy mean",
        "lbl_energy_min": "Energy minimum",
        "lbl_energy_sd": "Energy standard deviation",
        "mix_warn_pre": (
            "<b>Poor mixing warning</b>: only {} distinct orders among {} recorded "
            "samples — the chain has not actually traversed the state space; "
            "means/minima etc. must <b>not</b> be read as measures of ranking "
            "uncertainty (raise the temperature, lengthen the chain or increase "
            "burn-in and rerun)."
        ),
        "mcmc_corr_note": (
            "Chain samples are correlated, <b>not independent</b>; any mean/quantile "
            "interpretation must be read together with the effective sample size "
            "(ESS) and the state coverage fraction."
        ),
        "tree_no_order": "No ranking result.",
        "tree_title": "Ranking ladder",
        "tree_newick_h": "Ranking tree (Newick)",
        "of_th_item": "Artifact",
        "of_th_path": "Path",
        "of_informative": "Informative constraints (informative)",
        "of_conflicts": "Constraints conflicting with the best order (conflicts)",
        "of_partial": "Partial order (partial order)",
        "of_random_pre": (
            "Random-order objective distribution (random trees, 4th output file; "
            "short-style suffix {}, legacy-style suffix {}; generated only when "
            "--random-trees > 0)"
        ),
        "of_html_report": "Interactive HTML report",
        "of_not_generated": "(not generated in this run)",
    },
}


def _t(lang: str) -> Dict[str, str]:
    """取某语言的文案表；非法语言回退中文。"""
    return _LABELS.get(lang, _LABELS["zh"])


def _esc(s) -> str:
    """HTML/XML 转义（防止用户可控的标签注入 HTML/SVG）。

    ：旧实现写作 ``escape(str(s)) if s else ""``，把**整数 0**（以及
    任何假值）当成"没有内容"而渲染成空白单元格 —— 例如"众数位置 = 0"（根节点）
    在报告表格里就是空白。现仅对 ``None`` 返回空串，其余一律转义后输出。
    """
    if s is None:
        return ""
    return _html_escape(str(s))


# ----------------------------------------------------------------------
# 可选依赖探测（降级必须可见）
# ----------------------------------------------------------------------
_PLOTLY_PROBE: Optional[Tuple[bool, str]] = None

PLOTLY_MISSING_NOTE = (
    "plotly 未安装：交互式图表已降级为内联 SVG 静态图（信息不丢失，但不可缩放/悬停）。"
    '如需完整交互式报告，请执行 pip install "MaxTiC-Next[report]"。'
)


def has_plotly() -> bool:
    """当前环境能否使用 plotly（首次调用时探测并缓存）。"""
    return _probe_plotly()[0]


def plotly_missing_reason() -> str:
    """plotly 不可用时的原因文本（可用时为空串）。"""
    return _probe_plotly()[1]


def _probe_plotly() -> Tuple[bool, str]:
    global _PLOTLY_PROBE
    if _PLOTLY_PROBE is None:
        try:
            import plotly.graph_objects  # noqa: F401
            from plotly.offline import plot  # noqa: F401

            _PLOTLY_PROBE = (True, "")
        except Exception as exc:  # noqa: BLE001 - ImportError 之外还可能是版本不兼容
            _PLOTLY_PROBE = (False, f"{type(exc).__name__}: {exc}")
    return _PLOTLY_PROBE


def degradation_notice_html(lang: str = "zh") -> str:
    """依赖降级提示（可见横幅）；无降级时返回空串。"""
    t = _t(lang)
    notes: List[str] = []
    if not has_plotly():
        reason = plotly_missing_reason()
        notes.append(t["plotly_missing"] + (t["probe_info"].format(_esc(reason)) if reason else ""))
    if not notes:
        return ""
    items = "".join(f"<li>{_esc(t2)}</li>" for t2 in notes)
    return f"<div class='notice'><b>{_esc(t['report_degraded'])}</b><ul>{items}</ul></div>"


def weight_histogram(weights: List[float], bins: int = 12, lang: str = "zh") -> str:
    """约束权重分布直方图。Plotly 可用时返回交互式 div，否则返回内联 SVG。"""
    t = _t(lang)
    if has_plotly():
        try:
            import plotly.graph_objects as go
            from plotly.offline import plot as _plot

            fig = go.Figure(data=[go.Histogram(x=weights, nbinsx=bins)])
            fig.update_layout(
                title=t["hist_title"],
                xaxis_title=t["hist_x"],
                yaxis_title=t["hist_y"],
                margin=dict(l=40, r=20, t=40, b=40),
                height=320,
            )
            return _plot(fig, include_plotlyjs="inline", output_type="div")
        except Exception:  # noqa: BLE001 - 绘图期异常同样降级，不让报告失败
            return _weight_histogram_svg(weights, bins, t)
    return _weight_histogram_svg(weights, bins, t)


def _weight_histogram_svg(
    weights: List[float], bins: int = 12, t: Optional[Dict[str, str]] = None
) -> str:
    """权重分布的内联 SVG 直方图（Plotly 缺失时的降级实现）。"""
    if t is None:
        t = _t("zh")
    if not weights:
        return f"<p class='muted'>{t['no_weight_data']}</p>"
    lo, hi = min(weights), max(weights)
    if hi == lo:
        hi = lo + 1.0
    width, height = 560, 300
    pad_l, pad_b = 40, 30
    plot_w = width - pad_l - 10
    plot_h = height - pad_b - 10
    step = (hi - lo) / bins
    counts = [0] * bins
    for w in weights:
        idx = min(bins - 1, int((w - lo) / step)) if step > 0 else 0
        counts[idx] += 1
    max_count = max(counts) if counts else 1
    bar_w = plot_w / bins
    bars = []
    for i, c in enumerate(counts):
        if c == 0:
            continue
        x = pad_l + i * bar_w
        h = (c / max_count) * plot_h
        y = pad_b + (plot_h - h)
        bars.append(
            f'<rect x="{x:.1f}" y="{y:.1f}" width="{bar_w - 2:.1f}" height="{h:.1f}" '
            f'fill="#4C72B0"></rect>'
        )
    # 坐标轴
    axis = (
        f'<line x1="{pad_l}" y1="{pad_b}" x2="{pad_l}" y2="{pad_b + plot_h}" stroke="#333"/>'
        f'<line x1="{pad_l}" y1="{pad_b + plot_h}" x2="{pad_l + plot_w}" '
        f'y2="{pad_b + plot_h}" stroke="#333"/>'
        f'<text x="{pad_l}" y="{pad_b + plot_h + 18}" fill="#333" font-size="11">'
        f"{lo:.2f}</text>"
        f'<text x="{pad_l + plot_w - 30}" y="{pad_b + plot_h + 18}" fill="#333" '
        f'font-size="11">{hi:.2f}</text>'
    )
    return (
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" '
        f'viewBox="0 0 {width} {height}" role="img" aria-label="{_esc(t["hist_aria"])}">'
        f"{''.join(bars)}{axis}"
        f'<text x="{width // 2}" y="16" fill="#333" font-size="13" '
        f'text-anchor="middle">{_esc(t["hist_title"])}</text></svg>'
    )


def robustness_html(
    summary: Optional[Dict], max_pairs: int = PAIR_TABLE_ROWS, lang: str = "zh"
) -> str:
    """基于局部搜索访问解的稳健性/敏感性摘要的 HTML 表示。

    若未运行局部搜索（``summary is None``），返回占位说明。术语铁律：仅称
    "基于局部搜索访问解的稳健性/敏感性摘要"，杜绝后验/置信区间用语。

    ：同时展示**去重序数**与**真实访问次数/频率**（两者可以差几个数量级）。
    ：成对次序表默认只列"有信息量"的分歧对；树种系强制、方向恒定的平凡对
    （祖先-后代对在**所有**合法序中同向，计数必然等于全体访问数）被折叠为一行说明，
    不再挤占 top-20 名额。

    Args:
        summary: ``NearOptimalCollector.summary()`` 的返回值。
        max_pairs: 成对次序表最多展示的行数。
        lang: 报告语言（``"zh"`` / ``"en"``）。``summary`` 里的标题与口径说明
            由核心按中文生成，英文报告在此按同一模板重写，保证语义一致。
    """
    t = _t(lang)
    if summary is None:
        return f"<p class='muted'>{t['robust_placeholder']}</p>"
    # summary["title"] / ["basis"] 是核心按中文写死的；英文报告按同一语义改写，
    # 其余数值字段原样展示。
    title = t["robust_title"] if lang == "en" else summary.get("title", t["robust_title"])
    parts = [f"<h3>{_esc(title)}</h3>"]

    n_orders = summary.get("n_orders_collected", 0)
    acc_ret = summary.get("n_accesses_retained", None)
    acc_tot = summary.get("n_accesses_total", None)
    n_unique_total = summary.get("n_unique_orders_total", None)
    li = [
        f"<li>{t['li_orders_dedup'].format(_esc(n_orders))}</li>",
        f"<li>{t['li_orders_retained'].format(_esc(summary.get('n_orders_retained', n_orders)))}</li>",
    ]
    if acc_ret is not None:
        li.append(f"<li>{t['li_accesses_retained'].format(_esc(acc_ret))}</li>")
    if acc_tot is not None:
        li.append(f"<li>{t['li_accesses_total'].format(_esc(acc_tot))}</li>")
    if n_unique_total is not None:
        bound = t["lb_lower_bound"] if summary.get("n_unique_orders_total_is_lower_bound") else ""
        li.append(f"<li>{t['li_unique_total'].format(_esc(n_unique_total), bound)}</li>")
    li.append(
        f"<li>{t['li_top_k'].format(_esc(summary.get('top_k')), _esc(summary.get('top_k_reported', summary.get('top_k'))))}</li>"
    )
    li.append(f"<li>{t['li_best_value'].format(_esc(summary.get('best_value')))}</li>")
    parts.append("<ul>" + "".join(li) + "</ul>")
    if lang == "en":
        if summary.get("basis"):
            parts.append(
                f"<p class='muted'>{_esc(t['basis_fmt'].format(k=summary.get('top_k_reported', summary.get('top_k'))))}</p>"
            )
    elif summary.get("basis"):
        parts.append(f"<p class='muted'>{_esc(summary['basis'])}</p>")
    if summary.get("is_resampling_robustness") is False:
        parts.append(f"<p class='muted'>{_esc(t['robust_note'])}</p>")

    # 节点位置分布（紧凑表格：节点 / 众数位置 / 位置跨度 / 访问次数 / 频率）
    pos_dist = summary.get("node_position_distribution", {})
    pos_freq = summary.get("node_position_frequency", {})
    if pos_dist:
        parts.append(f"<h4>{_esc(t['pos_dist_h'])}</h4>")
        parts.append(
            "<table class='data'><thead><tr>"
            f"<th>{_esc(t['th_node'])}</th><th>{_esc(t['th_modal_pos'])}</th>"
            f"<th>{_esc(t['th_pos_span'])}</th>"
            f"<th>{_esc(t['th_visits'])}</th><th>{_esc(t['th_modal_share'])}</th>"
            "</tr></thead><tbody>"
        )
        rows = []
        for node, dist in pos_dist.items():
            if not dist:
                continue
            modal = max(dist, key=lambda p: (dist[p], -p))
            span = (max(dist) - min(dist)) + 1
            total = sum(dist.values())
            modal_share = dist[modal] / total if total else 0.0
            rows.append((node, modal, span, total, modal_share))
        rows.sort(key=lambda r: r[1])
        for node, modal, span, total, share in rows:
            parts.append(
                f"<tr><td>{_esc(node)}</td><td>{_esc(modal)}</td>"
                f"<td>{_esc(span)}</td><td>{_esc(total)}</td>"
                f"<td>{_esc(round(share, 4))}</td></tr>"
            )
        parts.append("</tbody></table>")
        if pos_freq:
            parts.append(f"<p class='muted'>{_esc(t['freq_note'])}</p>")

    # 成对次序：剔除"方向恒定"的平凡对（多数为树种系强制的祖先-后代对）
    pair_freq = summary.get("pairwise_order_frequency", {})
    if pair_freq:
        informative = summary.get("informative_pairwise_order_frequency")
        trivial = summary.get("constant_direction_pairs")
        if informative is None or trivial is None:
            informative, trivial = _split_trivial_pairs(pair_freq, summary)
        denom = summary.get("n_accesses_retained") or sum(pair_freq.values()) or 0
        parts.append(f"<h4>{t['pair_h']}</h4>")
        parts.append(
            "<p class='muted'>"
            + t["pair_note_pre"]
            + f"{_esc(len(trivial) // 2 if trivial else 0)}"
            + t["pair_note_post"].format(_esc(len(trivial or [])))
            + "</p>"
        )
        parts.append(
            "<table class='data'><thead><tr>"
            f"<th>{_esc(t['th_a'])}</th><th>{_esc(t['th_b'])}</th>"
            f"<th>{_esc(t['th_count'])}</th><th>{_esc(t['th_freq'])}</th>"
            "</tr></thead><tbody>"
        )
        ranked = sorted(informative.items(), key=lambda kv: kv[1], reverse=True)[:max_pairs]
        if not ranked:
            parts.append(f"<tr><td colspan='4' class='muted'>{_esc(t['pair_empty'])}</td></tr>")
        for (a, b), cnt in ranked:
            prop = (
                (cnt / denom)
                if denom
                else summary.get("pairwise_order_proportion", {}).get((a, b), 0.0)
            )
            parts.append(
                f"<tr><td>{_esc(a)}</td><td>{_esc(b)}</td>"
                f"<td>{_esc(cnt)}</td><td>{_esc(round(prop, 4))}</td></tr>"
            )
        parts.append("</tbody></table>")
        if trivial:
            shown = [
                f"{_esc(a)}&nbsp;≺&nbsp;{_esc(b)}"
                for a, b in sorted(trivial, key=lambda p: (str(p[0]), str(p[1])))[:12]
            ]
            more = "" if len(trivial) <= 12 else t["trivial_more"].format(len(trivial) - 12)
            parts.append(
                f"<details><summary>{_esc(t['trivial_details'])}</summary>"
                f"<p class='muted'>{'; '.join(shown)}{more}</p></details>"
            )

    return "\n".join(parts)


def _split_trivial_pairs(
    pair_freq: Dict[Tuple[str, str], int], summary: Dict
) -> Tuple[Dict[Tuple[str, str], int], List[List[str]]]:
    """在没有 ``informative_pairwise_order_frequency`` 时按计数重建平凡/信息量划分。

    兼容精简摘要 dict（只有 ``pairwise_order_frequency``）：若某对的正反两向计数
    中一个等于总访问数、另一个为 0，则该对方向恒定（平凡）。
    """
    prop = summary.get("pairwise_order_proportion") or {}
    total = summary.get("n_accesses_retained") or sum(pair_freq.values()) or 0
    trivial: List[List[str]] = []
    trivial_keys = set()
    for (a, b), p in prop.items():
        if p >= 1.0 - 1e-12:
            trivial.append([a, b])
            trivial_keys.add((a, b))
            trivial_keys.add((b, a))
    if not trivial and total:
        for (a, b), c in pair_freq.items():
            if c >= total and pair_freq.get((b, a), 0) == 0:
                trivial.append([a, b])
                trivial_keys.add((a, b))
                trivial_keys.add((b, a))
    informative = {k: v for k, v in pair_freq.items() if k not in trivial_keys}
    return informative, trivial


def mcmc_html(diag: Optional[Dict], max_rows: int = 12, lang: str = "zh") -> str:
    """Metropolis–Hastings 链的混合诊断（必须与统计量一并展示）。

    Args:
        diag: ``MCMCSampler.diagnostics()`` 的返回值（或同结构的 dict）；
            ``None`` / 空时返回占位说明。
        max_rows: 诊断表最多行数。
        lang: 报告语言（``"zh"`` / ``"en"``）。``status_note`` 由核心按中文生成时
            原样保留（它是运行时事实说明）；表格与说明文字按 ``lang`` 输出。
    """
    t = _t(lang)
    if not diag:
        return f"<p class='muted'>{t['mcmc_placeholder']}</p>"
    parts = [f"<h3>{_esc(t['mcmc_h'])}</h3>"]
    note = diag.get("status_note")
    if note:
        parts.append(f"<p class='notice'>{t['mcmc_status'].format(_esc(note))}</p>")
    order = [
        ("n_samples", t["lbl_n_samples"]),
        ("n_unique_samples", t["lbl_n_unique"]),
        ("unique_sample_fraction", t["lbl_unique_frac"]),
        ("adjacent_duplicate_fraction", t["lbl_adj_dup"]),
        ("max_chain_run_fraction", t["lbl_max_run"]),
        ("effective_sample_size", t["lbl_ess"]),
        ("ess_fraction", t["lbl_ess_frac"]),
        ("state_space_size", t["lbl_state_space"]),
        ("states_visited", t["lbl_states_visited"]),
        ("state_coverage_fraction", t["lbl_coverage"]),
        ("burn_in", t["lbl_burn_in"]),
        ("thin", t["lbl_thin"]),
        ("temperature", t["lbl_temperature"]),
        ("energy_mean", t["lbl_energy_mean"]),
        ("energy_min", t["lbl_energy_min"]),
        ("energy_sd", t["lbl_energy_sd"]),
    ]
    rows = []
    for key, label in order:
        if key not in diag or diag[key] is None:
            continue
        val = diag[key]
        if isinstance(val, float):
            val = f"{val:.6g}"
        rows.append(f"<tr><td>{_esc(label)}</td><td>{_esc(val)}</td></tr>")
    parts.append(
        "<table class='data'><thead><tr>"
        f"<th>{_esc(t['mcmc_th_diag'])}</th><th>{_esc(t['mcmc_th_value'])}</th>"
        f"</thead><tbody>{''.join(rows[: max_rows * 2])}</tbody></table>"
    )
    n = diag.get("n_samples") or 0
    uniq = diag.get("n_unique_samples") or 0
    if n and uniq * 5 < n:
        parts.append(f"<p class='notice'>{t['mix_warn_pre'].format(n, uniq)}</p>")
    parts.append(f"<p class='muted'>{_esc(t['mcmc_corr_note'])}</p>")
    return "\n".join(parts)


def tree_view_html(order: List[str], newick: str, lang: str = "zh") -> str:
    """排序树可视化：按 rank 的纵向阶梯（SVG）+ Newick 文本。

    说明：本可视化是排序位置的**占位可视化**（纵向阶梯表示 rank 1..n），完整
    系统发生树布局由外部工具渲染；此处保证报告自包含、可离线打开。
    """
    t = _t(lang)
    if not order:
        return f"<p class='muted'>{t['tree_no_order']}</p>"
    n = len(order)
    row_h = 22
    height = max(120, n * row_h + 40)
    width = 420
    rows = []
    for i, node in enumerate(order):
        y = 30 + i * row_h
        rows.append(
            f'<text x="20" y="{y}" fill="#333" font-size="12">{_esc(f"{i + 1}. {node}")}</text>'
            f'<line x1="20" y1="{y - 5}" x2="380" y2="{y - 5}" stroke="#ccc"/>'
        )
    svg = (
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" '
        f'viewBox="0 0 {width} {height}">'
        f"<title>{_esc(t['tree_title'])}</title>"
        f"{''.join(rows)}</svg>"
    )
    newick_block = f"<h4>{_esc(t['tree_newick_h'])}</h4><pre class='newick'>{_esc(newick)}</pre>"
    return svg + "\n" + newick_block
