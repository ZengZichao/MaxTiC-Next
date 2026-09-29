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


def degradation_notice_html() -> str:
    """依赖降级提示（可见横幅）；无降级时返回空串。"""
    notes: List[str] = []
    if not has_plotly():
        reason = plotly_missing_reason()
        notes.append(PLOTLY_MISSING_NOTE + (f"（探测信息：{_esc(reason)}）" if reason else ""))
    if not notes:
        return ""
    items = "".join(f"<li>{_esc(t)}</li>" for t in notes)
    return f"<div class='notice'><b>报告已降级渲染</b><ul>{items}</ul></div>"


def weight_histogram(weights: List[float], bins: int = 12) -> str:
    """约束权重分布直方图。Plotly 可用时返回交互式 div，否则返回内联 SVG。"""
    if has_plotly():
        try:
            import plotly.graph_objects as go
            from plotly.offline import plot as _plot

            fig = go.Figure(data=[go.Histogram(x=weights, nbinsx=bins)])
            fig.update_layout(
                title="约束权重分布",
                xaxis_title="约束权重",
                yaxis_title="约束数量",
                margin=dict(l=40, r=20, t=40, b=40),
                height=320,
            )
            return _plot(fig, include_plotlyjs="inline", output_type="div")
        except Exception:  # noqa: BLE001 - 绘图期异常同样降级，不让报告失败
            return _weight_histogram_svg(weights, bins)
    return _weight_histogram_svg(weights, bins)


def _weight_histogram_svg(weights: List[float], bins: int = 12) -> str:
    """权重分布的内联 SVG 直方图（Plotly 缺失时的降级实现）。"""
    if not weights:
        return "<p class='muted'>无权重数据。</p>"
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
        f'viewBox="0 0 {width} {height}" role="img" aria-label="权重分布直方图">'
        f"{''.join(bars)}{axis}"
        f'<text x="{width // 2}" y="16" fill="#333" font-size="13" '
        f'text-anchor="middle">约束权重分布</text></svg>'
    )


def robustness_html(summary: Optional[Dict], max_pairs: int = PAIR_TABLE_ROWS) -> str:
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
    """
    if summary is None:
        return (
            "<p class='muted'>未运行局部搜索（--local-search / --ls &gt; 0），"
            "无基于局部搜索访问解的稳健性/敏感性摘要。</p>"
        )
    title = summary.get("title", "基于局部搜索访问解的稳健性/敏感性摘要")
    parts = [f"<h3>{_esc(title)}</h3>"]

    n_orders = summary.get("n_orders_collected", 0)
    acc_ret = summary.get("n_accesses_retained", None)
    acc_tot = summary.get("n_accesses_total", None)
    n_unique_total = summary.get("n_unique_orders_total", None)
    li = [
        f"<li>去重近优排序数（参与统计）：<b>{_esc(n_orders)}</b></li>",
        f"<li>保留的去重近优排序数：<b>{_esc(summary.get('n_orders_retained', n_orders))}</b></li>",
    ]
    if acc_ret is not None:
        li.append(f"<li>上述排序的<b>真实访问次数</b>合计：{_esc(acc_ret)}</li>")
    if acc_tot is not None:
        li.append(f"<li>局部搜索累计访问（含被容量裁剪淘汰者）：{_esc(acc_tot)}</li>")
    if n_unique_total is not None:
        bound = "（下界）" if summary.get("n_unique_orders_total_is_lower_bound") else ""
        li.append(f"<li>累计见过的不同排序：{_esc(n_unique_total)}{bound}</li>")
    li.append(
        f"<li>保留上限 top_k：{_esc(summary.get('top_k'))}"
        f"（本次报告使用 {_esc(summary.get('top_k_reported', summary.get('top_k')))}）"
        "</li>"
    )
    li.append(f"<li>最优目标值：{_esc(summary.get('best_value'))}</li>")
    parts.append("<ul>" + "".join(li) + "</ul>")
    if summary.get("basis"):
        parts.append(f"<p class='muted'>{_esc(summary['basis'])}</p>")
    if summary.get("is_resampling_robustness") is False:
        parts.append(
            "<p class='muted'>注意：本摘要不是对输入数据重抽样的稳健度，"
            "也不是概率区间；它只描述近优解集合内部各节点次序的敏感性。</p>"
        )

    # 节点位置分布（紧凑表格：节点 / 众数位置 / 位置跨度 / 访问次数 / 频率）
    pos_dist = summary.get("node_position_distribution", {})
    pos_freq = summary.get("node_position_frequency", {})
    if pos_dist:
        parts.append("<h4>节点位置分布（按 rank 位置的真实访问次数）</h4>")
        parts.append(
            "<table class='data'><thead><tr>"
            "<th>节点</th><th>众数位置</th><th>位置跨度</th>"
            "<th>访问次数</th><th>众数位置占比</th></tr></thead><tbody>"
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
            parts.append("<p class='muted'>（占比按真实访问次数计算，非按去重排序数。）</p>")

    # 成对次序：剔除"方向恒定"的平凡对（多数为树种系强制的祖先-后代对）
    pair_freq = summary.get("pairwise_order_frequency", {})
    if pair_freq:
        informative = summary.get("informative_pairwise_order_frequency")
        trivial = summary.get("constant_direction_pairs")
        if informative is None or trivial is None:
            informative, trivial = _split_trivial_pairs(pair_freq, summary)
        denom = summary.get("n_accesses_retained") or sum(pair_freq.values()) or 0
        parts.append("<h4>成对次序频率（a 在 b 之前的<b>真实访问次数</b>）</h4>")
        parts.append(
            "<p class='muted'>下表已剔除方向恒定为 100% 的平凡对（多为物种树拓扑强制的"
            f"祖先-后代对，其方向在<b>所有</b>合法排序中都相同、不含排序信息），"
            f"以便看到真正有信息量的近优分歧。被剔除的平凡对共 "
            f"{_esc(len(trivial) // 2 if trivial else 0)} 对（{_esc(len(trivial or []))} "
            "个有向键）。</p>"
        )
        parts.append(
            "<table class='data'><thead><tr>"
            "<th>a</th><th>b</th><th>访问次数</th><th>频率</th></tr></thead><tbody>"
        )
        ranked = sorted(informative.items(), key=lambda kv: kv[1], reverse=True)[:max_pairs]
        if not ranked:
            parts.append(
                "<tr><td colspan='4' class='muted'>近优集合内所有成对次序方向"
                "均恒定（无可见的近优分歧）。</td></tr>"
            )
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
            more = "" if len(trivial) <= 12 else f" …（另有 {len(trivial) - 12} 个）"
            parts.append(
                "<details><summary>方向恒定的平凡对（前 12 个，点击展开）</summary>"
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


def mcmc_html(diag: Optional[Dict], max_rows: int = 12) -> str:
    """Metropolis–Hastings 链的混合诊断（必须与统计量一并展示）。

    Args:
        diag: ``MCMCSampler.diagnostics()`` 的返回值（或同结构的 dict）；
            ``None`` / 空时返回占位说明。
        max_rows: 诊断表最多行数。
    """
    if not diag:
        return (
            "<p class='muted'>未运行 Metropolis–Hastings 采样（--mcmc 未开启），无链混合诊断。</p>"
        )
    parts = ["<h3>Metropolis–Hastings 采样诊断</h3>"]
    note = diag.get("status_note")
    if note:
        parts.append(f"<p class='notice'><b>状态</b>：{_esc(note)}</p>")
    order = [
        ("n_samples", "记录样本数"),
        ("n_unique_samples", "唯一样本数"),
        ("unique_sample_fraction", "唯一样本比例"),
        ("adjacent_duplicate_fraction", "相邻重复率"),
        ("max_chain_run_fraction", "最长同态连续段占比"),
        ("effective_sample_size", "有效样本量 (ESS)"),
        ("ess_fraction", "ESS / 样本数"),
        ("state_space_size", "状态空间（线性扩展数）"),
        ("states_visited", "已访问状态数"),
        ("state_coverage_fraction", "状态访问比例"),
        ("burn_in", "burn-in 步数"),
        ("thin", "抽稀间隔"),
        ("temperature", "温度 T"),
        ("energy_mean", "能量均值"),
        ("energy_min", "能量最小值"),
        ("energy_sd", "能量标准差"),
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
        "<table class='data'><thead><tr><th>诊断量</th><th>值</th></tr>"
        f"</thead><tbody>{''.join(rows[: max_rows * 2])}</tbody></table>"
    )
    n = diag.get("n_samples") or 0
    uniq = diag.get("n_unique_samples") or 0
    if n and uniq * 5 < n:
        parts.append(
            f"<p class='notice'><b>混合不足警告</b>：{n} 个记录样本中只有 {uniq} 个不同"
            "排序，链实际上没有遍历状态空间；此时的均值/最小值等统计量<b>不能</b>"
            "用作排序不确定性的度量（提高温度、加长链或增大 burn-in 后重跑）。</p>"
        )
    parts.append(
        "<p class='muted'>链上样本彼此相关，<b>不是独立样本</b>；任何均值/分位数解释"
        "必须连同有效样本量 (ESS) 与状态访问比例一起阅读。</p>"
    )
    return "\n".join(parts)


def tree_view_html(order: List[str], newick: str) -> str:
    """排序树可视化：按 rank 的纵向阶梯（SVG）+ Newick 文本。

    说明：本可视化是排序位置的**占位可视化**（纵向阶梯表示 rank 1..n），完整
    系统发生树布局由外部工具渲染；此处保证报告自包含、可离线打开。
    """
    if not order:
        return "<p class='muted'>无排序结果。</p>"
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
        f"<title>排序阶梯</title>"
        f"{''.join(rows)}</svg>"
    )
    newick_block = f"<h4>排序树 Newick</h4><pre class='newick'>{_esc(newick)}</pre>"
    return svg + "\n" + newick_block
