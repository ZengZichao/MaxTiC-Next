"""pytest 公共配置：确保 ``src`` 在 ``sys.path`` 中，并提供示例数据路径 fixtures。"""

import os
import sys

import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC = os.path.join(ROOT, "src")
if SRC not in sys.path:
    sys.path.insert(0, SRC)

DATA_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")


def pytest_configure(config):
    """确保 src 目录在模块搜索路径中（兼容无 pip install 的运行）。"""
    if SRC not in sys.path:
        sys.path.insert(0, SRC)


@pytest.fixture(scope="session", autouse=True)
def _clean_maxtic_outputs():
    """会话级自动清理：移除 ``tests/data`` 顶层由测试产生的 MaxTiC 输出文件。

    部分测试在共享的 ``tests/data`` 约束文件上调用 ``rank()`` 且不传
    ``output_prefix``，会把三输出文件 / distribution / HTML 报告写入 ``tests/data``。
    本 fixture 在会话开始与结束各扫除一次，仅删除**输出后缀**的文件（legacy +
    short 两种命名风格及 ``.html``），输入文件（``*.tsv`` 约束 / ``*.tree``）不受影响，
    保持数据目录整洁（与 ``.gitignore`` 一致）。
    """
    from maxtic_next.config import OUTPUT_SUFFIXES

    markers = tuple(s.lower() for group in OUTPUT_SUFFIXES.values() for s in group) + (".html",)

    def _sweep():
        try:
            names = os.listdir(DATA_DIR)
        except OSError:
            return
        for n in names:
            p = os.path.join(DATA_DIR, n)
            if os.path.isfile(p) and n.lower().endswith(markers):
                try:
                    os.remove(p)
                except OSError:
                    pass

    _sweep()  # 会话开始：清理上次运行留下的陈旧输出
    yield
    _sweep()  # 会话结束：清理本次产生的输出


@pytest.fixture
def data_dir():
    return DATA_DIR


@pytest.fixture
def minitree_path():
    return os.path.join(DATA_DIR, "minitree.tree")


@pytest.fixture
def cyano_path():
    return os.path.join(DATA_DIR, "Cyano_CUTConstraints.tsv")
