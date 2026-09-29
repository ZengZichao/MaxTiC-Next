#!/usr/bin/env python3
"""仓库边界回归测试：本仓库与桌面端仓库 MaxTiC-Next-Studio 的关系怎么描述。

桌面端是**独立的 GitHub 仓库**（MaxTiC-Next-Studio），两者只通过包依赖相连。
文档里若用本地目录名或相对路径（``../MaxTiC-Next-Studio-项目代码`` 之类）描述它，
换一台机器、换一个克隆布局就全是错的 —— 所以这里把它变成可执行的约束。
"""

from __future__ import annotations

from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]

STUDIO_REPO_URL = "https://github.com/ZengZichao/MaxTiC-Next-Studio"

# 只允许出现在"描述与另一仓库的关系"的文档里。
_LOCAL_PATH_MARKERS = (
    "../MaxTiC",
    "MaxTiC-Next-Studio-项目代码",
    "同级目录",
    "同级仓库",
    "同级独立仓库",
    "sibling directory",
    "sibling repository",
)

FRONT_END_DOCS = (
    "README.md",
    "README.en.md",
    "docs/manual/zh/00_index.md",
    "docs/manual/zh/03_cli_reference.md",
    "docs/manual/zh/07_workflows_deployment.md",
    "docs/manual/en/00_index.md",
    "docs/manual/en/03_cli_reference.md",
    "docs/manual/en/07_workflows_deployment.md",
)


@pytest.mark.parametrize("rel", FRONT_END_DOCS)
def test_docs_never_point_at_a_local_folder(rel):
    text = (ROOT / rel).read_text(encoding="utf-8")
    hits = [m for m in _LOCAL_PATH_MARKERS if m in text]
    assert not hits, f"{rel} 用本地路径/本地目录名描述了桌面端仓库：{hits}"


def test_readme_links_the_desktop_repository():
    for rel in ("README.md", "README.en.md"):
        text = (ROOT / rel).read_text(encoding="utf-8")
        assert STUDIO_REPO_URL in text, f"{rel} 没给出桌面端仓库地址"


def test_this_repository_holds_no_gui_code_or_dependency():
    """边界是双向的：核心仓库既不含 GUI 代码，也不声明 GUI 依赖。

    查的是 pyproject 的**实际字段**而不是全文 —— 注释里提到图形前端位于 Studio 仓库
    时出现 PySide6 字样是正确且有用的信息，不该被判为越界。
    """
    try:
        import tomllib
    except ModuleNotFoundError:  # Python < 3.11
        import tomli as tomllib

    assert not (ROOT / "src" / "maxtic_next" / "gui").exists()
    assert not (ROOT / "packaging").exists()
    with open(ROOT / "pyproject.toml", "rb") as fh:
        cfg = tomllib.load(fh)
    project = cfg["project"]
    deps = [d.lower() for d in project.get("dependencies", [])]
    extras = {
        name: [d.lower() for d in listed]
        for name, listed in project.get("optional-dependencies", {}).items()
    }
    all_deps = deps + [d for listed in extras.values() for d in listed]
    for banned in ("pyside6", "qt", "matplotlib"):
        assert not any(banned in d for d in all_deps), f"核心仓库声明了 GUI 依赖：{banned}"
    assert "studio" not in extras, "核心仓库仍有 studio extras"
    assert "maxtic-studio" not in project.get("scripts", {}), "核心仓库仍注册 GUI 入口"
    assert project["dependencies"] == [], "核心包承诺零运行依赖，请勿引入"
