"""移植等价自洽回归测试（ported-equivalence self-consistency）。

⚠️ 重要声明（非"机器级字节比对原版 Py2"）：
本测试**不是**把本项目输出与原版 Python 2.7 ``MaxTiC.py`` 做机器级逐字节比对。
本环境无 Python 2.7，无法运行原版做字节级对比。这里采用**移植等价自洽
（ported-equivalence self-consistency）** 验收标准：

1. 忠实逐函数移植、逻辑与源码一一对应（见各模块实现）；
2. 在固定 ``--seed 42`` 下，``minitree.tree + Cyano_CUTConstraints.tsv`` 跑通，
   生成三个输出文件 + stdout 摘要且无报错；
3. 排序树为合法拓扑序（每个内部节点在其后代之前）；
4. 将本次新代码在 seed=42 下的 stdout 与三文件记录为 ``tests/baselines/`` 快照，
   作为后续复现基线；首次运行生成基线，之后运行做精确比对。

也就是说：基线记录的是"本移植版本自身的确定性输出快照"，用于防止回归
（同一份代码在固定 seed 下必须可复现）。**真正的字节级等价验证**应由 CI 在
容器化 Python 2.7 环境中运行原版 ``MaxTiC.py`` 与本版本比对，
不在本测试范围内。

运行方式：``pytest tests/test_equivalence.py``（会自动生成基线）。

⚠️ 基线与原版 ``MaxTiC.py`` 的**已裁定偏离**：
``tests/baselines/`` 记录的是"本移植版本自身在 ``seed=42`` 下的确定性输出"，
用于防回归；它**有意不同于**原版输出，差异全部来自下列口径裁定：

* ``uninformative`` 百分比的分母改为 ``total_weight``（原版把无信息权重
  双重计入分母）：同一实例上原版打印 ``12%``，本基线打印 ``13.9%``；
* 偏序产物的哨兵边判定用 ``MAX_NUMBER``（1e10）而不是魔法数 ``100000``：
  权重 >= 1e5 的真实信息性约束不再被静默从 ``partial_order`` 里删掉；
* 信息性文件在"删除 0 权重边"**之后**写出（原版执行顺序；本基线的
  ``informative`` 中不含 0 权重行）；
* 局部搜索改进最优解后回写 ``values`` / ``best_source`` 并补打原版
  ``after local search … rejected`` / ``best found solution …`` 两行（本实例
  ``--ls 0``，故基线中不出现这两行，但 ``values``/摘要口径已随之改变）；
* ``--rd`` 的 p 值用 ``(k+1)/(n+1)`` 校正且检验**交付序**（本实例
  ``--rd 0``，基线不含 p 值行）；
* 谱系硬约束 NOTE 的措辞："…并已计入**上方** uninformative 统计（在
  total_weight 分母中只计一次，不重复计权）"。该 NOTE 打在摘要末尾，统计行在其
  上方，旧文案"计入下方"与打印顺序相反；原版则根本没有这一行（原版
  ``MaxTiC.py:444`` 的逐条命中打印在其执行路径上不可达）。

因此**任何"与原版逐字节等价"的表述都必须限定为**："除上述已裁定偏离外逐字节等价"。
真正的字节级对照须在容器化 Python 2.7 环境跑原版（见
``tests/test_equivalence_stats.py`` 的定位逻辑）。
"""

import os
import re
from decimal import Decimal, ROUND_HALF_EVEN

from maxtic_next.api import rank
from maxtic_next.tree.tree import Tree

DATA_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")
BASELINE_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "baselines")

# 项目根目录（tests/ 的上一级）；用于把基线中的绝对路径归一化，
# 使基线对目录重命名 / 迁移不敏感（仅算法输出进入字节级比对）。
PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def _normalize(text):
    """把项目根目录的绝对路径替换为稳定占位符 ``<PROJECT_ROOT>``。"""
    return text.replace(PROJECT_ROOT, "<PROJECT_ROOT>")


# 仅匹配真正的浮点数（含小数点或科学计数法）；整数（节点计数、百分号前的数字等）
# 不匹配，保持逐字节精确比对。
_FLOAT_RE = re.compile(r"-?\d+\.\d+(?:[eE][-+]?\d+)?|-?\d+[eE][-+]?\d+")
# 保留 12 位有效数字（对 ~1e3 量级容差约 1e-9）。足以吸收不同 numpy / Python 版本带来的
# 浮点末位噪声（相对 ~1e-13，如 ``2218.4`` 与 ``2218.3999999999996``），同时仍能捕获
# 真实算法回归（这类变化远大于 1e-9）。
_SIG_FIGS = 12


def _round_sig(token: str) -> str:
    """把单个浮点数字符串归一到固定有效数字（Decimal 精确舍入，避免 float 末位噪声）。"""
    d = Decimal(token)
    if d == 0:
        return "0.0"
    q = _SIG_FIGS - (d.adjusted() + 1)  # 保留 _SIG_FIGS 位有效数字所需的小数位数
    quantized = d.quantize(Decimal(1).scaleb(-q), rounding=ROUND_HALF_EVEN)
    return format(quantized, "f")


def _normalize_floats(text: str) -> str:
    """将文本中的浮点数字归一到固定有效数字，吸收末位浮点噪声。

    等价性基线用于捕获「算法回归」（固定 seed 下输出必须可复现），而非锁定特定
    numpy / Python 版本的浮点末位表示。不同环境下 ``2218.4`` 与 ``2218.3999999999996``
    这类末位噪声应被视为等价。整数不匹配浮点正则，保持精确比对，确保结构 / 计数变化仍
    能被捕获。
    """
    return _FLOAT_RE.sub(lambda m: _round_sig(m.group(0)), text)


def _read(path):
    with open(path, "r", encoding="utf-8") as fh:
        return fh.read()


def _check_or_create(name, content):
    """首次运行创建基线；之后精确比对。

    比对前将项目根目录的绝对路径归一化为 ``<PROJECT_ROOT>`` 占位符，
    使基线对目录重命名 / 迁移不敏感（仅算法输出进入比对）。
    """
    path = os.path.join(BASELINE_DIR, name)
    os.makedirs(BASELINE_DIR, exist_ok=True)
    content = _normalize_floats(_normalize(content))
    if not os.path.exists(path):
        with open(path, "w", encoding="utf-8") as fh:
            fh.write(content)
        return
    existing = _normalize_floats(_normalize(_read(path)))
    assert existing == content, f"基线不一致：{name}"


def _assert_topological(order, tree, internal_labels):
    """断言 order 是合法拓扑序：每个内部节点在其所有后代之前。"""
    index = {label: i for i, label in enumerate(order)}
    for label in internal_labels:
        nid = tree.label_to_node_id(label)
        # 叶子不参与 order，故只检查内部后代（见下方 stack 遍历）
        # 检查内部后代
        stack = list(tree.get_children(nid))
        while stack:
            child = stack.pop()
            if not tree.is_leaf(child):
                child_label = tree.get_bootstrap(child)
                assert index[label] < index[child_label], (
                    f"拓扑序违反：{label} 应早于其后裔 {child_label}"
                )
                stack.extend(tree.get_children(child))


def test_equivalence_minictree_cyano(tmp_path, capsys):
    tree_path = os.path.join(DATA_DIR, "minitree.tree")
    cons_path = os.path.join(DATA_DIR, "Cyano_CUTConstraints.tsv")
    prefix = str(tmp_path / "Cyano_CUTConstraints.tsv")

    # 载入树用于拓扑校验
    tree = Tree()
    with open(tree_path) as fh:
        tree.read_newick(fh.readline())
    internal_labels = tree.internal_node_labels()

    # 等价性校验以原版命名（legacy）为准，忠实复刻原 MaxTiC.py 的输出文件名
    result = rank(
        tree_path,
        cons_path,
        seed=42,
        output_prefix=prefix,
        print_summary=True,
        output_style="legacy",
    )
    captured = capsys.readouterr().out

    # 三文件应存在（legacy 长名）
    assert os.path.exists(result.informative_file)
    assert os.path.exists(result.conflicting_file)
    assert os.path.exists(result.partial_order_file)
    assert result.informative_file.endswith(
        "_MT_output_filtered_list_of_weighted_informative_constraints"
    )

    # 合法性校验
    assert result.best_source in ("greedy heuristic", "mixing heuristic")
    _assert_topological(result.best_order, tree, internal_labels)

    # 记录 / 比对基线
    def _join(lines):
        return "\n".join(lines) + ("\n" if lines else "")

    _check_or_create("stdout.txt", captured)
    _check_or_create("informative.txt", _join(result.informative_lines))
    _check_or_create("conflicting.txt", _join(result.conflicting_lines))
    _check_or_create("partial.txt", _join(result.partial_lines))
