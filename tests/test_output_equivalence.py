"""逐字节输出等价性基线。

**为什么需要这个文件**："不改动数值顺序" 是硬约束——`mix` 的平局要消耗 `rng.random()`
（见 ``ranking/mixing.py``），浮点末位一旦抖动就会改变排序结果，进而破坏可复现性契约。
单条 CLI 冒烟路径的覆盖面不足，不足以给 ``edgeweights`` 集合化与 ``--incremental`` 放行。

本模块对一组**固定用例**（示例数据 + 确定性合成树 × 多种子 × 多开关）产出规范化摘要，
与 ``tests/baselines/output_equivalence.json`` 比对：

    pytest tests/test_output_equivalence.py              # 校验模式（默认，CI 用）
    MAXTIC_UPDATE_OUTPUT_BASELINE=1 pytest tests/test_output_equivalence.py
                                                          # 重新生成

规范化处理：仓库根、临时目录与输出前缀分别替换为 ``<ROOT>`` / ``<TMP>`` / ``<OUT>``——
摘要首行和 ``run_metadata["constraint_files_merged"]`` 都会回显输入文件的绝对路径，
不归一化仓库根就会让基线绑死在生成它的那个目录上，换一台机器（或 CI）必然误报漂移；
``run_metadata`` 里的墙钟类字段（``time_for_search`` / ``local_search_elapsed_sec`` /
``*_seconds`` / ``*_sec``）被排除，因为它们天生不稳定。

浮点一律归一到 **12 位有效数字** 再取摘要，与 ``tests/test_equivalence.py`` 的
``_normalize_floats`` 同一口径：本基线看守的是"固定 seed 下算法输出可复现"，而不是
某个 CPython 版本的浮点末位表示。同一份代码在 3.11 上 ``total_weight`` 的 ``repr`` 是
``2218.3999999999996``、在 3.14 上是 ``2218.4``（相对 ~1e-16 的末位噪声），若按全精度
``repr`` 取摘要，CI 的 3.9–3.14 矩阵就只能在其中某一个版本上通过——那是版本脆弱性而非
回归。整数不参与归一（浮点正则要求小数点或指数），故节点计数与结构变化仍逐字节捕获，
真实算法回归（变化远大于 1e-9）同样会失败。
"""

from __future__ import annotations

import hashlib
import json
import os
import random
import re
import sys
from contextlib import redirect_stdout
from decimal import Decimal, ROUND_HALF_EVEN
from io import StringIO
from typing import Dict, List, Tuple

import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
SRC = os.path.join(ROOT, "src")
if SRC not in sys.path:
    sys.path.insert(0, SRC)

BASELINE = os.path.join(HERE, "baselines", "output_equivalence.json")
UPDATE = os.environ.get("MAXTIC_UPDATE_OUTPUT_BASELINE") == "1"

from maxtic_next import api  # noqa: E402

EXAMPLE_TREE = os.path.join(ROOT, "examples", "minitree.tree")
EXAMPLE_CONS = os.path.join(ROOT, "examples", "Cyano_CUTConstraints.tsv")

# 天生不稳定的墙钟类元数据字段：不进摘要（否则每次运行都"漂移"）。
# 实测键名包括 time_for_search 与 local_search_elapsed_sec —— 后者是 --incremental
# 与全量重算之间**唯一**的差异字段，正因如此必须排除，否则会误判为数值分歧。
VOLATILE_META = {
    "time_for_search",
    "elapsed_seconds",
    "wall_seconds",
    "local_search_elapsed_sec",
    "elapsed_sec",
}


def _is_volatile(key: str) -> bool:
    return (
        key in VOLATILE_META
        or key.endswith("_seconds")
        or key.endswith("_elapsed")
        or key.endswith("_elapsed_sec")
        or key.endswith("_sec")
    )


# 与 tests/test_equivalence.py 同口径的浮点归一化：只匹配真正的小数/科学计数（整数不匹配，
# 因此节点计数与结构差异仍逐字节捕获）。
_FLOAT_RE = re.compile(r"-?\d+\.\d+(?:[eE][-+]?\d+)?|-?\d+[eE][-+]?\d+")
_SIG_FIGS = 12


def _round_sig(token: str) -> str:
    """把单个浮点数字符串归一到固定有效数字（Decimal 精确舍入，避免 float 末位噪声）。"""
    d = Decimal(token)
    if d == 0:
        return "0.0"
    q = _SIG_FIGS - (d.adjusted() + 1)
    quantized = d.quantize(Decimal(1).scaleb(-q), rounding=ROUND_HALF_EVEN)
    return format(quantized, "f")


def _normalize_floats(text: str) -> str:
    """吸收不同 CPython / numpy 版本的浮点末位表示差异（见模块文档）。"""
    return _FLOAT_RE.sub(lambda m: _round_sig(m.group(0)), text)


# --------------------------------------------------------------------------- #
# 确定性合成输入（生成器见本文件 make_inputs；规模取小以保证 CI 快）
# --------------------------------------------------------------------------- #
def make_inputs(n_leaves: int, density: int, seed: int, outdir: str) -> Tuple[str, str]:
    """生成平衡二叉物种树 + 随机约束；路径只由入参决定（同参数复用同一份文件）。"""
    tag = f"syn-n{n_leaves}-d{density}-s{seed}"
    tp, cp = os.path.join(outdir, tag + ".tree"), os.path.join(outdir, tag + ".tsv")
    if os.path.isfile(tp) and os.path.isfile(cp):
        return tp, cp
    rnd = random.Random(seed)
    labels = [str(i) for i in range(1, n_leaves + 1)]
    nodes = [x + ":1" for x in labels]
    nxt, internals = n_leaves + 1, []
    while len(nodes) > 1:
        new = []
        for i in range(0, len(nodes) - 1, 2):
            name = str(nxt)
            nxt += 1
            internals.append(name)
            new.append(f"({nodes[i]},{nodes[i + 1]}){name}:1")
        if len(nodes) % 2:
            new.append(nodes[-1])
        nodes = new
    all_nodes = labels + internals
    cons = [
        f"{a}\t{b}\t{rnd.uniform(0.1, 30):.2f}"
        for a, b in (rnd.sample(all_nodes, 2) for _ in range(int(len(all_nodes) * density)))
    ]
    with open(tp, "w", encoding="utf-8") as fh:
        fh.write(nodes[0] + ";\n")
    with open(cp, "w", encoding="utf-8") as fh:
        fh.write("\n".join(cons) + "\n")
    return tp, cp


# (用例名, 输入规格, 传给 api.rank 的额外参数)；输入规格 = ("example",) 或 (叶数, 密度, 种子)
CASES: List[Tuple[str, tuple, dict]] = [
    ("example_default", ("example",), dict(seed=42)),
    ("example_seed7", ("example",), dict(seed=7)),
    (
        "example_local_search",
        ("example",),
        dict(seed=42, local_search=5.0, temperature=0.001, local_search_max_iterations=3000),
    ),
    (
        "example_ls_incremental",
        ("example",),
        dict(
            seed=42,
            local_search=5.0,
            temperature=0.001,
            local_search_max_iterations=3000,
            incremental=True,
        ),
    ),
    ("example_threshold", ("example",), dict(seed=42, threshold_constraints=0.15)),
    ("example_random_trees", ("example",), dict(seed=3, random_trees=5)),
    ("example_near_optimal", ("example",), dict(seed=11, top_k=5)),
    ("example_legacy_style", ("example",), dict(seed=42, output_style="legacy")),
    ("syn_small_s1", (60, 4, 1), dict(seed=1)),
    ("syn_small_s42", (60, 4, 42), dict(seed=42)),
    ("syn_small_s7", (60, 4, 7), dict(seed=7)),
    ("syn_mid_s1", (120, 4, 1), dict(seed=1)),
    (
        "syn_mid_ls",
        (120, 4, 1),
        dict(seed=1, local_search=5.0, temperature=0.001, local_search_max_iterations=4000),
    ),
    (
        "syn_mid_ls_incremental",
        (120, 4, 1),
        dict(
            seed=1,
            local_search=5.0,
            temperature=0.001,
            local_search_max_iterations=4000,
            incremental=True,
        ),
    ),
    ("syn_mid_mcmc", (120, 4, 5), dict(seed=5, mcmc=True, mcmc_iters=3000, mcmc_burn_in=500)),
    ("syn_mid_random_type1", (120, 4, 9), dict(seed=9, random_type=1)),
    # 同一输入、同一种子，只有 incremental 不同 → 这两条摘要必须逐字节相同
    (
        "pair_ls_full",
        (120, 4, 1),
        dict(
            seed=1,
            local_search=5.0,
            temperature=0.001,
            local_search_max_iterations=4000,
            incremental=False,
        ),
    ),
    (
        "pair_ls_incremental",
        (120, 4, 1),
        dict(
            seed=1,
            local_search=5.0,
            temperature=0.001,
            local_search_max_iterations=4000,
            incremental=True,
        ),
    ),
]


def _norm(text: str, tmpdir: str, prefix: str) -> str:
    """把环境相关字符串归一化，其余逐字节保留。"""
    if prefix:
        text = text.replace(prefix, "<OUT>")
    for token in {tmpdir, os.path.realpath(tmpdir)}:
        text = text.replace(token, "<TMP>")
    # 示例输入是仓库内的绝对路径，会出现在摘要首行与 constraint_files_merged 中；
    # 连软链解析形式一起归一化，基线才与仓库 clone 到何处无关。
    for token in {ROOT, os.path.realpath(ROOT)}:
        text = text.replace(token, "<ROOT>")
    return text


def _run(name: str, spec: tuple, kwargs: dict, tmpdir: str) -> Dict[str, str]:
    if spec[0] == "example":
        tree, cons = EXAMPLE_TREE, [EXAMPLE_CONS]
    else:
        n_leaves, density, seed = spec
        tree, cpath = make_inputs(n_leaves, density, seed, tmpdir)
        cons = [cpath]
    prefix = os.path.join(tmpdir, "out", name)
    os.makedirs(os.path.dirname(prefix), exist_ok=True)

    buf = StringIO()
    with redirect_stdout(buf):
        r = api.rank(
            tree, cons, output_prefix=prefix, print_summary=True, html_report=False, **kwargs
        )
    summary = buf.getvalue()
    assert not r.html_report_file, f"{name}: 意外产出 HTML 报告"

    def sha(s: str) -> str:
        # 归一浮点末位再取摘要：同一份代码在 3.11 与 3.14 上必须得到相同的基线。
        return hashlib.sha256(_normalize_floats(s).encode("utf-8")).hexdigest()[:20]

    files = {}
    for label, path in (
        ("informative", r.informative_file),
        ("conflicts", r.conflicting_file),
        ("partial", r.partial_order_file),
    ):
        with open(path, encoding="utf-8") as fh:
            files[label] = sha(_norm(fh.read(), tmpdir, prefix))

    meta = {k: v for k, v in sorted((r.run_metadata or {}).items()) if not _is_volatile(k)}
    payload = {
        "files": files,
        "stdout": sha(_norm(summary, tmpdir, prefix)),
        "ranked_newick": sha(_norm(r.ranked_newick, tmpdir, prefix)),
        "best_order": sha(json.dumps(list(r.best_order), ensure_ascii=False)),
        "values": sha(json.dumps({k: repr(v) for k, v in sorted(r.values.items())})),
        "counts": sha(
            json.dumps(
                {
                    "internal_node_count": r.internal_node_count,
                    "total_weight": repr(r.total_weight),
                    "uninformative": {k: repr(v) for k, v in sorted(r.uninformative.items())},
                    "trivial_conflict": repr(r.trivial_conflict),
                    "informative_count": repr(r.informative_count),
                    "similarity_to_input": repr(r.similarity_to_input),
                    "best_source": r.best_source,
                },
                sort_keys=True,
            )
        ),
        "warnings": sha(json.dumps([_norm(str(w), tmpdir, prefix) for w in (r.warnings or [])])),
        # ensure_ascii=False：默认转义会把非 ASCII 的仓库路径变成 \uXXXX，
        # 令 _norm 的字符串替换失效，基线就会绑死在生成它的目录上。
        "metadata": sha(
            _norm(
                json.dumps(meta, sort_keys=True, default=repr, ensure_ascii=False), tmpdir, prefix
            )
        ),
    }
    return payload


@pytest.fixture(scope="module")
def tmp_root(tmp_path_factory) -> str:
    return str(tmp_path_factory.mktemp("output_equivalence"))


def _collect(tmpdir: str) -> Dict[str, Dict[str, str]]:
    return {name: _run(name, spec, kwargs, tmpdir) for name, spec, kwargs in CASES}


def test_update_or_verify_output_equivalence(tmp_root):
    """改动前后的产物必须逐字节一致；基线缺失时明确失败而不是静默通过。"""
    current = _collect(tmp_root)
    if UPDATE:
        os.makedirs(os.path.dirname(BASELINE), exist_ok=True)
        payload = {
            "_comment": (
                "由 MAXTIC_UPDATE_OUTPUT_BASELINE=1 生成。任何代码改动都必须让本"
                "测试继续通过；确需改变数值行为时，须逐条说明偏离。"
            ),
            "cases": current,
        }
        with open(BASELINE, "w", encoding="utf-8") as fh:
            json.dump(payload, fh, ensure_ascii=False, indent=2, sort_keys=True)
        assert len(current) == len(CASES)
        return
    if not os.path.isfile(BASELINE):
        pytest.fail(
            f"缺少基线文件 {BASELINE}；请在**未改动代码**时运行 "
            "MAXTIC_UPDATE_OUTPUT_BASELINE=1 pytest tests/test_output_equivalence.py"
        )
    with open(BASELINE, encoding="utf-8") as fh:
        baseline = json.load(fh)["cases"]
    diffs = []
    for name, exp in sorted(baseline.items()):
        got = current.get(name)
        if got is None:
            diffs.append(f"{name}: 用例被删除")
            continue
        diffs += [
            f"{name}.{k}: {exp[k]} -> {got.get(k)}" for k in sorted(exp) if got.get(k) != exp[k]
        ]
    diffs += [f"{name}: 新增用例未登记基线" for name in sorted(set(current) - set(baseline))]
    assert not diffs, "输出逐字节漂移：\n" + "\n".join(diffs)


def test_incremental_must_be_byte_identical_to_full_recomputation(tmp_root):
    """判据：同输入同种子下 ``--incremental`` 与全量重算**逐字节一致**。

    刻意独立于基线文件——直接对比两条路径，即使有人重生成基线也无法掩盖两路径分歧。
    """
    common = dict(seed=1, local_search=5.0, temperature=0.001, local_search_max_iterations=4000)
    a = _run("pair_ls_full", (120, 4, 1), dict(common, incremental=False), tmp_root)
    b = _run("pair_ls_incremental", (120, 4, 1), dict(common, incremental=True), tmp_root)
    diff = [k for k in sorted(a) if a[k] != b[k]]
    assert not diff, (
        "incremental 与全量重算产出不同，差异字段："
        + ", ".join(diff)
        + "\n→ 分歧修好之前 `--incremental` 不得转为默认开启"
    )


def test_baseline_keeps_both_incremental_paths():
    """基线里必须同时存在 incremental 开/关的成对用例，否则这条安全网形同虚设。"""
    ls_cases = [(n, k) for n, _, k in CASES if k.get("local_search")]
    assert len(ls_cases) >= 4, "局部搜索用例太少，无法看守 --incremental"
    assert any(k.get("incremental") is True for _, k in ls_cases)
    assert any(k.get("incremental") in (None, False) for _, k in ls_cases)
