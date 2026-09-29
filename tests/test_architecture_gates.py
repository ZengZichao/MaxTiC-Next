"""架构与文档一致性的可执行门禁。

设计取向：**"要不要拆环"这个昂贵问题，换成"别让环变多"这个零成本问题。**
本文件不修架构，只把当前已核实、且各有明确理由的边界冻结成断言：

* ：包级双向依赖必须**恰好等于**已登记的四组；新增一对即失败。
* / ：numpy 位集分支与纯 Python 回退必须产出**逐字节相同**的结果。
* 护栏：核心 `dependencies` 必须保持为空数组（这是该工具在生信分发生态里的资产）。

全部只用标准库，不依赖 `../架构可视化/` 的脚本。
"""

from __future__ import annotations

import ast
import hashlib
import os
import re
from pathlib import Path

# 本仓库承诺支持 Python 3.9+（CI 矩阵 3.9–3.14），而 tomllib 是 3.11 才进标准库的；
# 缺它时整个模块会在收集阶段 ImportError，把 3.9/3.10 的全部用例一起打死。
try:
    import tomllib
except ModuleNotFoundError:  # pragma: no cover — Python < 3.11
    import tomli as tomllib

import pytest

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
SRC = ROOT / "src" / "maxtic_next"

# 改前由 架构可视化/data/import-graph.json 的 51 条包级边算出，且逐条核对过动机
# （类型标注接缝 / 报告副作用 / MCMC 复用能量函数 / 压缩读取归属）。
ALLOWED_PACKAGE_CYCLES: frozenset[frozenset[str]] = frozenset(
    {
        frozenset({"io", "constraints"}),
        frozenset({"io", "ranking"}),
        frozenset({"ranking", "report"}),
        frozenset({"ranking", "robustness"}),
    }
)


def _module_id(path: Path) -> str:
    """`src/maxtic_next/io/parsing.py` -> `maxtic_next.io.parsing`（含包 __init__ 归并）。"""
    rel = os.path.relpath(path, SRC.parent).replace(os.sep, ".")[:-3]
    return rel[: -len(".__init__")] if rel.endswith(".__init__") else rel


def _pkg_of(module: str) -> str:
    """与 架构可视化/data/import-graph.json 同口径：`maxtic_next.io.parsing` -> `io`，
    顶层模块（`maxtic_next.api`）自成一格，包自身（`maxtic_next`）记作 `(root)`。"""
    parts = module.split(".")
    return parts[1] if len(parts) >= 2 else "(root)"


def _package_edges() -> set[tuple[str, str]]:
    """AST 扫描 src/，返回包级 import 边集合（含函数内的惰性 import）。"""
    edges: set[tuple[str, str]] = set()
    for path in sorted(SRC.rglob("*.py")):
        if "__pycache__" in path.parts:
            continue
        src_mod = _pkg_of(_module_id(path))
        tree = ast.parse(path.read_text(encoding="utf-8"))
        for node in ast.walk(tree):
            targets: list[str] = []
            if isinstance(node, ast.ImportFrom) and node.module and node.level == 0:
                targets = [node.module]
            elif isinstance(node, ast.Import):
                targets = [a.name for a in node.names]
            for t in targets:
                if not t.startswith("maxtic_next"):
                    continue
                edges.add((src_mod, _pkg_of(t)))
    return {(a, b) for a, b in edges if a != b}


def test_no_new_package_cycles():
    """：包级双向依赖必须恰好是已登记的四组（新增即失败，减少则要求同步更新此常量）。"""
    edges = _package_edges()
    found = {frozenset((a, b)) for (a, b) in edges if (b, a) in edges}
    extra = found - ALLOWED_PACKAGE_CYCLES
    gone = ALLOWED_PACKAGE_CYCLES - found
    assert not extra, (
        "新增了包级双向依赖（架构约束要求冻结）："
        + ", ".join("↔".join(sorted(c)) for c in extra)
        + "\n→ 若确有必要，请连同 ALLOWED_PACKAGE_CYCLES "
        "一起更新，并在 PR 里说明接缝动机。"
    )
    assert not gone, (
        "包级双向依赖减少了但仍被登记为允许："
        + ", ".join("↔".join(sorted(c)) for c in gone)
        + " → 请把相关结论与 ALLOWED_PACKAGE_CYCLES 一起收紧。"
    )


def test_module_level_scc_members_are_known():
    """8 模块强连通分量的成员冻结（改这一簇要意识到其余 7 处都被牵动）。"""
    graph: dict[str, set[str]] = {}
    for path in sorted(SRC.rglob("*.py")):
        if "__pycache__" in path.parts:
            continue
        me = _module_id(path)
        tree = ast.parse(path.read_text(encoding="utf-8"))
        for node in ast.walk(tree):
            mods: list[str] = []
            if isinstance(node, ast.ImportFrom) and node.module and node.level == 0:
                # 只取被 import 的**模块**本身；`from pkg import name` 里的 name 可能是
                # 函数/类，展开它会把不存在的节点塞进图里（与 import-graph 口径保持一致）。
                mods = [node.module]
            elif isinstance(node, ast.Import):
                mods = [a.name for a in node.names]
            for m in mods:
                if m.startswith("maxtic_next") and m != me:
                    graph.setdefault(me, set()).add(m)
    known_scc = {
        "maxtic_next",
        "maxtic_next.api",
        "maxtic_next.dry_run",
        "maxtic_next.io",
        "maxtic_next.io.output",
        "maxtic_next.io.parsing",
        "maxtic_next.ranking.ranker",
        "maxtic_next.report.html",
    }
    sccs = _strong_components(graph, {n for n in graph})
    biggest = max((s for s in sccs if len(s) > 1), key=len, default=set())
    assert biggest <= known_scc, (
        f"import 环里出现了新的模块：{sorted(biggest - known_scc)}（架构约束已登记的成员见测试内注释）"
    )


def _strong_components(adj: dict[str, set[str]], nodes) -> list[set[str]]:
    index, low, stack, onstack, out = {}, {}, [], set(), []
    counter = [0]
    # 刻意不调 sys.setrecursionlimit：本实现是显式栈迭代，而且改全局递归上限会
    # 泄漏到同会话后续用例（会让 tests/test_correctness_guards.py 的深树用例误判）。

    def visit(start: str) -> None:
        work: list[tuple[str, iter]] = [(start, iter(sorted(adj.get(start, ()))))]
        index[start] = low[start] = counter[0]
        counter[0] += 1
        stack.append(start)
        onstack.add(start)
        while work:
            v, it = work[-1]
            pushed = False
            for w in it:
                if w not in index:
                    index[w] = low[w] = counter[0]
                    counter[0] += 1
                    stack.append(w)
                    onstack.add(w)
                    work.append((w, iter(sorted(adj.get(w, ())))))
                    pushed = True
                    break
                if w in onstack:
                    low[v] = min(low[v], index[w])
            if pushed:
                continue
            work.pop()
            if work:
                low[work[-1][0]] = min(low[work[-1][0]], low[v])
            if low[v] == index[v]:
                comp = set()
                while True:
                    w = stack.pop()
                    onstack.discard(w)
                    comp.add(w)
                    if w == v:
                        break
                out.append(comp)

    for n in sorted(nodes):
        if n not in index:
            visit(n)
    return out


# --------------------------------------------------------------------------- #
# / ：numpy 位集分支 vs 纯 Python 回退必须逐字节等价
# --------------------------------------------------------------------------- #
def _digest_run(tmp_path: Path, monkeypatch, *, fake_numpy: bool, seed: int) -> dict:
    """在指定"numpy 在场/缺席"视图下跑一次排序，返回产物摘要。"""
    from maxtic_next import api
    from maxtic_next.ranking import reachability as reach

    monkeypatch.setattr(reach, "_HAVE_NUMPY", fake_numpy, raising=True)
    tree = ROOT / "examples" / "minitree.tree"
    cons = ROOT / "examples" / "Cyano_CUTConstraints.tsv"
    prefix = str(tmp_path / f"numpy{int(fake_numpy)}-seed{seed}")
    r = api.rank(
        str(tree),
        [str(cons)],
        seed=seed,
        print_summary=False,
        html_report=False,
        output_prefix=prefix,
        force=True,
    )
    sha = lambda s: hashlib.sha256(s.encode("utf-8")).hexdigest()[:20]  # noqa: E731
    files = {}
    for label, p in (
        ("informative", r.informative_file),
        ("conflicts", r.conflicting_file),
        ("partial", r.partial_order_file),
    ):
        files[label] = sha(Path(p).read_text(encoding="utf-8"))
    return {
        "files": files,
        "order": sha(",".join(r.best_order)),
        "values": sha(repr(sorted((k, repr(v)) for k, v in r.values.items()))),
        "ranked_newick": sha(r.ranked_newick),
    }


@pytest.mark.parametrize("seed", [7, 42, 1234])
def test_numpy_bitset_branch_matches_pure_python(tmp_path, monkeypatch, seed) -> None:
    """两条可达性实现必须给出一模一样的产物。

    ``fake_numpy=True`` 这一侧要真的执行位集实现，因此**必须有 numpy 在场**：
    ``reachability`` 在缺 numpy 时把模块名 ``np`` 置为 ``None``，只翻转 ``_HAVE_NUMPY``
    会让 ``np.zeros`` 直接 ``AttributeError``——那是被测环境的缺陷，不是分支分歧。
    所以缺 numpy 的作业里本用例跳过，由 CI 矩阵的 ``numpy=true`` 腿覆盖两条分支；
    ``numpy=false`` 腿里 ``fake_numpy=False`` 一侧（纯 Python 回退）本就是该环境的默认路径，
    全套件都在跑它。
    """
    from maxtic_next.ranking import reachability as reach

    if not reach._HAVE_NUMPY:
        pytest.skip("numpy 未安装：位集分支无法真实执行，交由 CI 的 numpy 作业覆盖")

    a = _digest_run(tmp_path, monkeypatch, fake_numpy=False, seed=seed)
    b = _digest_run(tmp_path, monkeypatch, fake_numpy=True, seed=seed)
    assert a == b, f"numpy 位集分支与纯 Python 回退结果不一致（seed={seed}）：{a} != {b}"


# --------------------------------------------------------------------------- #
# 护栏：核心零第三方依赖
# --------------------------------------------------------------------------- #
def test_core_runtime_dependencies_stay_empty() -> None:
    cfg = tomllib.loads((ROOT / "pyproject.toml").read_text(encoding="utf-8"))
    assert cfg["project"]["dependencies"] == [], (
        "核心 `dependencies` 必须保持为空 —— 零依赖是对外的承诺，改动需同步更新文档。"
    )


# --------------------------------------------------------------------------- #
# 交付面一致性：版本号必须五处齐平
# --------------------------------------------------------------------------- #
def _version_bearing_texts() -> dict[str, str]:
    return {
        "pyproject.toml": (ROOT / "pyproject.toml").read_text(encoding="utf-8"),
        "src/maxtic_next/__init__.py": (SRC / "__init__.py").read_text(encoding="utf-8"),
        "bioconda/meta.yaml": (ROOT / "bioconda" / "meta.yaml").read_text(encoding="utf-8"),
        "Singularity.def": (ROOT / "Singularity.def").read_text(encoding="utf-8"),
        "Dockerfile": (ROOT / "Dockerfile").read_text(encoding="utf-8"),
    }


def test_version_is_declared_consistently() -> None:
    """版本号在 5 处声明处必须相同；Docker 侧必须可用 `ARG MAXTIC_VERSION` 校验。

    改前：`0.1.0` 在四处一致，但 Dockerfile 只在注释里出现 tag（无 ARG/LABEL），
    镜像版本也就无从校验。
    """
    texts = _version_bearing_texts()
    declared = {
        "pyproject.toml": re.search(r'^version = "([^"]+)"', texts["pyproject.toml"], re.M).group(
            1
        ),
        "src/maxtic_next/__init__.py": re.search(
            r'__version__ = "([^"]+)"', texts["src/maxtic_next/__init__.py"]
        ).group(1),
        # 配方用 jinja 变量声明版本（`{% set version = "0.1.0" %}`），
        # `version: {{ version }}` 只是引用它。
        "bioconda/meta.yaml": re.search(
            r'\{%\s*set version\s*=\s*"([\w.]+)"\s*%\}', texts["bioconda/meta.yaml"]
        ).group(1),
        "Singularity.def": re.search(
            r"^\s*Version\s+([\w.]+)", texts["Singularity.def"], re.M
        ).group(1),
        "Dockerfile": re.search(r"ARG MAXTIC_VERSION=([\w.]+)", texts["Dockerfile"]).group(1),
    }
    assert len(set(declared.values())) == 1, f"版本号不一致：{declared}"
    assert "org.opencontainers.image.version=" in texts["Dockerfile"], (
        "镜像缺少 OCI 版本标签，无法校验'镜像里是哪个版本'"
    )
