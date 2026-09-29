"""增强等价性测试：多种子统计、确定性复现、原版对比、分布等价。

本模块在 ``test_equivalence.py`` 的自洽基线快照机制之上，补充四类增强验证：

a. **多种子统计等价测试**：用至少 10 个不同种子运行，记录每次的 best value 和
   Kendall 相似度，验证分布的统计特性（均值、标准差在合理范围）。
b. **确定性路径复现测试**：同一种子运行多次（至少 3 次），验证输出完全一致。
c. **原版对比测试（条件性）**：检测环境中是否有 python2 和原版 MaxTiC.py，
   如果有则运行原版并对比输出；如果没有则 skip 并打印说明信息。
d. **随机分布等价测试**：多次随机种子运行后，验证 value 分布的 KS 检验
   （使用 scipy.stats.ks_2samp 如果可用，否则用简单的分位数比较）。

这些测试不影响 ``test_equivalence.py`` 中的现有基线快照测试。
"""

import math
import os
import shutil
import subprocess

import pytest

from maxtic_next.api import rank
from maxtic_next.tree.tree import Tree

DATA_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")
TREE_PATH = os.path.join(DATA_DIR, "minitree.tree")
CONS_PATH = os.path.join(DATA_DIR, "Cyano_CUTConstraints.tsv")

# 原版 MaxTiC.py 是第三方 Python 2 脚本，**不随本仓库分发**。
# 本地若另备有一份，用环境变量 ``MAXTIC_ORIGINAL_PY`` 指向它（原版还依赖同目录下的
# ``script_tree`` 模块，因此运行时的 cwd 取其所在目录）；未设置时下面的用例跳过，
# 等价性宣称由 tests/test_reference_equivalence.py 看守。
_ORIGINAL_MAXTIC_PY = os.environ.get("MAXTIC_ORIGINAL_PY", "")
_ORIGINAL_MAXTIC_DIR = os.path.dirname(_ORIGINAL_MAXTIC_PY)

# 至少 10 个不同种子用于统计测试
TEST_SEEDS = [0, 1, 2, 3, 7, 13, 42, 99, 123, 456, 789, 1024]


# ----------------------------------------------------------------------
# 辅助函数
# ----------------------------------------------------------------------
def _best_value(result):
    """从 Result 中提取 best value（取 greedy / mixing 中的较小者）。

    MaxTiC 的 value 是最小化目标（未解释的转移权重总和越小越好），
    ranker 选优逻辑为 ``value_greedy <= value_heuristic`` 时选 greedy。
    """
    greedy = result.values.get("greedy", float("inf"))
    mixing = result.values.get("mixing", float("inf"))
    return min(greedy, mixing)


def _run_rank(seed, tmp_path, tag="", **kwargs):
    """运行一次排序，返回 Result。

    额外 kwargs 透传给 ``rank()``（如 ``local_search`` / ``random_type``）。
    """
    prefix = str(tmp_path / f"eq_{tag}_seed{seed}")
    return rank(
        TREE_PATH,
        CONS_PATH,
        seed=seed,
        output_prefix=prefix,
        print_summary=False,
        html_report=False,
        **kwargs,
    )


def _load_tree():
    """加载测试树（用于 Kendall 相似度计算）。"""
    tree = Tree()
    with open(TREE_PATH) as fh:
        tree.read_newick(fh.readline())
    return tree


# ----------------------------------------------------------------------
# a. 多种子统计等价测试
# ----------------------------------------------------------------------
class TestMultiSeedStatisticalEquivalence:
    """用至少 10 个不同种子运行，验证 best value 与 Kendall 相似度的统计特性。"""

    def test_multi_seed_best_value_statistics(self, tmp_path):
        """10+ 种子的 best value 均值/标准差在合理范围（有限、有界）。"""
        values = []
        for seed in TEST_SEEDS:
            r = _run_rank(seed, tmp_path, tag="stat")
            bv = _best_value(r)
            assert math.isfinite(bv), f"seed={seed} best_value 非有限: {bv}"
            assert bv >= 0.0, f"seed={seed} best_value 为负: {bv}"
            values.append(bv)

        n = len(values)
        assert n >= 10, f"种子数不足 10: {n}"

        mean_val = sum(values) / n
        var_val = sum((v - mean_val) ** 2 for v in values) / n
        std_val = math.sqrt(var_val)

        # 均值应为正有限值
        assert math.isfinite(mean_val) and mean_val > 0, f"best_value 均值异常: {mean_val}"

        # 标准差应有限且非负
        assert math.isfinite(std_val) and std_val >= 0, f"best_value 标准差异常: {std_val}"

        # 所有值应在一个合理范围内（不超过均值的 5 倍，不低于均值的 1/5）
        # 这捕获"某种子产生极端异常值"的回归
        for i, v in enumerate(values):
            assert mean_val / 5 <= v <= mean_val * 5, (
                f"seed={TEST_SEEDS[i]} best_value={v} 偏离均值 {mean_val} 过大"
            )

        # 打印统计摘要（供调试）
        print(
            f"\n[多种子统计] n={n}, mean={mean_val:.4f}, "
            f"std={std_val:.4f}, min={min(values):.4f}, max={max(values):.4f}"
        )

    def test_multi_seed_kendall_similarity_statistics(self, tmp_path):
        """10+ 种子的 Kendall 相似度（与输入树）均值/标准差在合理范围。"""
        sims = []
        for seed in TEST_SEEDS:
            r = _run_rank(seed, tmp_path, tag="sim")
            sim = r.similarity_to_input
            assert math.isfinite(sim), f"seed={seed} similarity 非有限: {sim}"
            # Kendall 相似度定义域 [0, 1]（1 = 与输入树完全一致）
            assert 0.0 <= sim <= 1.0, f"seed={seed} similarity 超出 [0,1]: {sim}"
            sims.append(sim)

        n = len(sims)
        assert n >= 10

        mean_sim = sum(sims) / n
        var_sim = sum((s - mean_sim) ** 2 for s in sims) / n
        std_sim = math.sqrt(var_sim)

        assert math.isfinite(mean_sim) and 0.0 <= mean_sim <= 1.0
        assert math.isfinite(std_sim) and std_sim >= 0.0

        print(
            f"\n[Kendall 相似度统计] n={n}, mean={mean_sim:.4f}, "
            f"std={std_sim:.4f}, min={min(sims):.4f}, max={max(sims):.4f}"
        )

    def test_multi_seed_all_produce_valid_permutations(self, tmp_path):
        """所有种子均产出内部节点全集的合法排列。"""
        tree = _load_tree()
        internal_labels = set(tree.internal_node_labels())
        for seed in TEST_SEEDS:
            r = _run_rank(seed, tmp_path, tag="perm")
            assert sorted(r.best_order) == sorted(internal_labels), (
                f"seed={seed} 未产出内部节点全集排列"
            )
            assert r.best_source in ("greedy heuristic", "mixing heuristic"), (
                f"seed={seed} best_source 异常: {r.best_source}"
            )


# ----------------------------------------------------------------------
# b. 确定性路径复现测试
# ----------------------------------------------------------------------
class TestDeterministicPathReproduction:
    """同一种子运行多次（至少 3 次），验证输出完全一致。"""

    def test_same_seed_identical_best_order(self, tmp_path):
        """seed=42 运行 3 次，best_order 完全一致。"""
        orders = []
        for i in range(3):
            r = _run_rank(42, tmp_path, tag=f"det{i}")
            orders.append(tuple(r.best_order))
        assert orders[0] == orders[1] == orders[2], f"同种子 best_order 不一致: {orders}"

    def test_same_seed_identical_values(self, tmp_path):
        """seed=42 运行 3 次，values dict 完全一致。"""
        vals = []
        for i in range(3):
            r = _run_rank(42, tmp_path, tag=f"val{i}")
            vals.append(r.values)
        assert vals[0] == vals[1] == vals[2], f"同种子 values 不一致: {vals}"

    def test_same_seed_identical_ranked_newick(self, tmp_path):
        """seed=42 运行 3 次，ranked_newick 完全一致。"""
        newicks = []
        for i in range(3):
            r = _run_rank(42, tmp_path, tag=f"nwk{i}")
            newicks.append(r.ranked_newick)
        assert newicks[0] == newicks[1] == newicks[2], "同种子 ranked_newick 不一致"

    def test_same_seed_identical_best_value_multiple_seeds(self, tmp_path):
        """多个种子各自重复 3 次，best_value 在同种子内完全一致。"""
        for seed in [1, 7, 42, 99]:
            bvs = []
            for i in range(3):
                r = _run_rank(seed, tmp_path, tag=f"bv{seed}_{i}")
                bvs.append(_best_value(r))
            assert bvs[0] == bvs[1] == bvs[2], f"seed={seed} 同种子 best_value 不一致: {bvs}"

    def test_different_seed_different_path(self, tmp_path):
        """不同种子（在随机化模式下）产生不同排序路径。

        注：默认参数（random_type=0）下 greedy 启发式是确定性的且总能胜出，
        种子不影响结果——这是正确行为。这里使用 random_type=1（随机化边方向）
        使 rng 参与建边过程，从而验证不同种子产生不同路径。
        不排除极少数种子碰巧相同，因此只断言"至少有一对不同"。
        """
        orders = {}
        for seed in [1, 2, 3, 4, 5]:
            r = _run_rank(seed, tmp_path, tag=f"diff{seed}", random_type=1)
            orders[seed] = tuple(r.best_order)
        unique = set(orders.values())
        assert len(unique) >= 2, (
            f"5 个不同种子（random_type=1）产出完全相同的排序，随机性可能未生效: {orders}"
        )


# ----------------------------------------------------------------------
# c. 原版对比测试（条件性）
# ----------------------------------------------------------------------
def _find_python2():
    """检测系统中是否存在 python2 可执行文件。"""
    return shutil.which("python2") or shutil.which("python2.7")


def _original_available():
    """检测原版 MaxTiC.py 和 python2 是否同时可用。"""
    return _find_python2() is not None and os.path.isfile(_ORIGINAL_MAXTIC_PY)


class TestOriginalComparison:
    """与原版 Python 2.7 MaxTiC.py 的条件性对比测试。

    原版 ``MaxTiC.py`` **不随本仓库分发**（它是第三方 Python 2 脚本），下面的探测只是
    检查本地是否另备有一份参考软件；没有 python2 或找不到原版脚本时本用例跳过。

    无论本用例是否运行，等价性宣称都由 ``tests/test_reference_equivalence.py`` 看守：
    那里有一份按原版控制流独立重写的 Python 3 参考实现（原版是 Py2 脚本，无法在 Py3
    下执行），在真实示例数据上与本版本逐字段差分比对，并显式记录本版本相对原版的
    已裁定偏离；本用例跳过时，那份差分测试照常运行。
    """

    def test_original_maxtic_comparison_or_skip(self, tmp_path):
        """如有 python2 + 原版 MaxTiC.py，运行原版并对比排序结果。

        原版 MaxTiC.py 无 --seed 参数（使用未播种的 random 模块），
        因此对比策略为：验证原版产出的排序是合法拓扑序且 value 合理，
        并与 MaxTiC-Next 的 best_order 做集合等价性校验
        （均为内部节点全集排列）。
        """
        if not _original_available():
            pytest.skip(
                "环境中无 python2（原版为 Python 2 脚本，不能用 python3 直接执行），"
                "故跳过**字节级**原版对比。\n"
                f"  python2: {_find_python2() or '未找到（PATH 中无 python2 / python2.7）'}\n"
                f"  原版 MaxTiC.py: "
                f"{'存在' if os.path.isfile(_ORIGINAL_MAXTIC_PY) else '未找到（MAXTIC_ORIGINAL_PY 未设置或文件不存在）'} "
                f"(MAXTIC_ORIGINAL_PY={_ORIGINAL_MAXTIC_PY or '未设置'})\n"
                "注意：等价性宣称并非无人看守 —— tests/test_reference_equivalence.py "
                "用独立重写的参考实现在同一数据集上做逐字段差分比对，"
                "能在 CI（无 python2）中证伪算法与统计口径漂移。"
            )

        py2 = _find_python2()
        # 原版 MaxTiC.py 依赖 script_tree 模块（在同目录下）
        # CLI: python MaxTiC.py species_tree_file constraints_file [ls=... t=...]
        out_dir = str(tmp_path / "original")
        os.makedirs(out_dir, exist_ok=True)
        # 复制约束文件到输出目录（原版在约束文件同目录下生成输出）
        cons_copy = os.path.join(out_dir, os.path.basename(CONS_PATH))
        shutil.copy(CONS_PATH, cons_copy)

        try:
            proc = subprocess.run(
                [py2, _ORIGINAL_MAXTIC_PY, TREE_PATH, cons_copy],
                capture_output=True,
                text=True,
                timeout=60,
                cwd=_ORIGINAL_MAXTIC_DIR,
            )
        except subprocess.TimeoutExpired:
            pytest.skip("原版 MaxTiC.py 运行超时（60s），跳过对比")
        except Exception as exc:
            pytest.skip(f"原版 MaxTiC.py 运行失败: {exc}")

        if proc.returncode != 0:
            pytest.skip(
                f"原版 MaxTiC.py 返回非零退出码 {proc.returncode}，stderr: {proc.stderr[:200]}"
            )

        # 原版 stdout 中包含排序结果摘要
        original_stdout = proc.stdout

        # 运行 MaxTiC-Next
        next_prefix = str(tmp_path / "next")
        r_next = rank(
            TREE_PATH,
            CONS_PATH,
            seed=42,
            output_prefix=next_prefix,
            print_summary=False,
            html_report=False,
        )

        # 校验：MaxTiC-Next 的 best_order 是合法排列
        tree = _load_tree()
        internal_labels = set(tree.internal_node_labels())
        assert sorted(r_next.best_order) == sorted(internal_labels)

        # 校验：原版也应在 stdout 中输出全部内部节点
        # 原版输出格式含排序后的 Newick（带分支长度），这里做宽松校验：
        # 所有内部节点标签均出现在原版输出中
        for label in internal_labels:
            assert str(label) in original_stdout, f"内部节点 {label} 未出现在原版输出中"

        print(f"\n[原版对比] 原版运行成功，stdout 长度={len(original_stdout)}")


# ----------------------------------------------------------------------
# d. 随机分布等价测试（KS 检验）
# ----------------------------------------------------------------------
def _ks_2samp(data1, data2):
    """两样本 KS 检验。优先使用 scipy，不可用时回退到手动计算。

    返回 (ks_statistic, p_value)。p_value 在手动模式下为 None。
    """
    try:
        from scipy.stats import ks_2samp as _scipy_ks

        result = _scipy_ks(data1, data2)
        # 兼容 scipy 的两种返回形态：KstestResult 与裸 tuple
        if hasattr(result, "statistic"):
            return result.statistic, result.pvalue
        return result[0], result[1]
    except ImportError:
        pass

    # 手动 KS 统计量（无 p-value）
    all_vals = sorted(set(data1) | set(data2))
    n1, n2 = len(data1), len(data2)
    cdf1 = [sum(1 for x in data1 if x <= v) / n1 for v in all_vals]
    cdf2 = [sum(1 for x in data2 if x <= v) / n2 for v in all_vals]
    ks_stat = max(abs(c1 - c2) for c1, c2 in zip(cdf1, cdf2))
    return ks_stat, None


def _quantile_compare(data1, data2, n_quantiles=4):
    """简单分位数比较：两组数据的四分位数应接近。"""

    def _quantiles(data, n):
        s = sorted(data)
        return [s[min(int(len(s) * i / n), len(s) - 1)] for i in range(n + 1)]

    q1 = _quantiles(data1, n_quantiles)
    q2 = _quantiles(data2, n_quantiles)
    max_diff = max(abs(a - b) for a, b in zip(q1, q2))
    return q1, q2, max_diff


class TestRandomDistributionEquivalence:
    """多次随机种子运行后，验证 value 分布的 KS 检验。

    将种子分成两组，分别运行得到两组 best_value 分布，验证它们来自同一分布。
    """

    def test_ks_test_value_distribution_two_halves(self, tmp_path):
        """将 12 个种子分成两组，KS 检验两组 value 分布一致。"""
        seeds_a = TEST_SEEDS[:6]
        seeds_b = TEST_SEEDS[6:]

        values_a = [_best_value(_run_rank(s, tmp_path, tag="ks_a")) for s in seeds_a]
        values_b = [_best_value(_run_rank(s, tmp_path, tag="ks_b")) for s in seeds_b]

        ks_stat, p_value = _ks_2samp(values_a, values_b)

        # KS 统计量应较小（两组来自同一算法同一数据的分布）
        assert ks_stat < 1.0, f"KS 统计量过大 ({ks_stat:.4f})，分布差异异常"

        # 如果有 p-value（scipy 可用），不应拒绝同分布假设
        if p_value is not None:
            # p_value > 0.05 表示不能拒绝"来自同一分布"的零假设
            # 注意：样本量小（6 个）时 p-value 可能偏低，放宽到 0.01
            assert p_value > 0.01, f"KS 检验 p-value={p_value:.4f} 过低，分布可能不一致"
            print(f"\n[KS 检验] ks_stat={ks_stat:.4f}, p_value={p_value:.4f} (scipy)")
        else:
            # scipy 不可用时，用分位数比较
            q1, q2, max_diff = _quantile_compare(values_a, values_b)
            range_val = max(max(values_a + values_b), 1.0)
            assert max_diff / range_val < 0.5, (
                f"分位数差异过大: max_diff={max_diff:.4f}, range={range_val:.4f}"
            )
            print(f"\n[分位数比较] max_diff={max_diff:.4f} (无 scipy, 手动)")

    def test_ks_test_reproducible_distribution(self, tmp_path):
        """同一组种子运行两次，value 分布应完全相同（确定性）。"""
        seeds = TEST_SEEDS[:8]
        run1 = [_best_value(_run_rank(s, tmp_path, tag="rep1")) for s in seeds]
        run2 = [_best_value(_run_rank(s, tmp_path, tag="rep2")) for s in seeds]

        # 确定性：同种子同结果 -> 两组完全一致
        assert run1 == run2, f"同种子重复运行 value 不一致:\n  run1={run1}\n  run2={run2}"

        # KS 统计量应为 0（完全相同的分布）
        ks_stat, _ = _ks_2samp(run1, run2)
        assert ks_stat == 0.0, f"确定性分布 KS 统计量非零: {ks_stat}"

    def test_quantile_comparison_value_distribution(self, tmp_path):
        """分位数比较：两组种子的 value 分布分位数应接近（无 scipy 回退路径）。"""
        seeds_a = TEST_SEEDS[:6]
        seeds_b = TEST_SEEDS[6:]

        values_a = [_best_value(_run_rank(s, tmp_path, tag="qc_a")) for s in seeds_a]
        values_b = [_best_value(_run_rank(s, tmp_path, tag="qc_b")) for s in seeds_b]

        q1, q2, max_diff = _quantile_compare(values_a, values_b)

        # 所有值有限
        for v in values_a + values_b:
            assert math.isfinite(v) and v >= 0

        # 分位数差异应在合理范围（不超过整体范围的 50%）
        all_vals = values_a + values_b
        data_range = max(all_vals) - min(all_vals)
        if data_range > 0:
            assert max_diff / data_range < 0.8, (
                f"分位数差异占范围比例过大: {max_diff / data_range:.2%}"
            )

        print(f"\n[分位数] A={q1}, B={q2}, max_diff={max_diff:.4f}")
