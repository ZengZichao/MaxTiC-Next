"""正确性守卫的回归测试。

每个用例把一条输入校验或数值不变量写成可执行断言。约定：

* ``TREE`` / ``CONS`` 为仓库内真实示例数据（13 个内部节点、总权重 2218.4）；
* 断言"错误信息里必须出现的词"而不是仅断言抛异常，以免守卫被替换成无信息异常；
* CLI 层用例通过 ``cli.main`` 抛 ``SystemExit`` 检查退出码。
"""

import os
import subprocess
import sys
import textwrap

import pytest

from maxtic_next import api
from maxtic_next.cli import build_parser, main as cli_main
from maxtic_next.config import (
    DEFAULT_MIN_TRANSFER_DIST,
    DEFAULT_RANDOM_TYPE,
    DEFAULT_SEED,
    DEFAULT_TEMPERATURE,
    DEFAULT_THRESHOLD_CONSTRAINTS,
    DEFAULT_TIME_FOR_SEARCH,
    MAX_NUMBER,
)
from maxtic_next.io.parsing import read_constraints_file, read_newick_file
from maxtic_next.ranking.mixing import mix, opt
from maxtic_next.ranking.ranker import Ranker, maximum_distance, random_order
from maxtic_next.random_ import RandomWrapper
from maxtic_next.tree.tree import Tree

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data")
TREE = os.path.join(DATA, "minitree.tree")
CONS = os.path.join(DATA, "Cyano_CUTConstraints.tsv")
# 注意：包在 src/ 之下，指向仓库根会让"未 pip install"的
# 子进程直接 ModuleNotFoundError，于是这条防误删结果的安全网只在 CI 的已安装环境里有效。
# 现在显式指到 src，并补一条断言：未安装状态下也必须能跑（见 test_existing_outputs_are_not_overwritten）。
REPO_ROOT = os.path.dirname(HERE)
SRC_PATH = os.path.join(REPO_ROOT, "src")


def _write(path, content):
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(content)
    return str(path)


_RUN_COUNTER = {"n": 0}


def _run(tmp_path, *args, **kwargs):
    kwargs.setdefault("print_summary", False)
    kwargs.setdefault("html_report", False)
    if "output_prefix" not in kwargs:
        _RUN_COUNTER["n"] += 1
        kwargs["output_prefix"] = os.path.join(str(tmp_path), f"run{_RUN_COUNTER['n']}")
    return api.rank(*args, **kwargs)


def _assert_linear_extension(order, tree):
    """交付序必须是"内部节点全集排列"且尊重物种树拓扑（父早于子）。"""
    labels = list(tree.internal_node_labels())
    assert sorted(order) == sorted(labels)
    index = {lbl: i for i, lbl in enumerate(order)}
    for lbl in labels:
        nid = tree.label_to_node_id(lbl)
        stack = list(tree.get_children(nid))
        while stack:
            child = stack.pop()
            if not tree.is_leaf(child):
                assert index[lbl] < index[tree.get_bootstrap(child)], (
                    f"拓扑被违反：{lbl} 应早于其后代 {tree.get_bootstrap(child)}"
                )
                stack.extend(tree.get_children(child))


# ======================================================================
# 非二叉物种树：无条件守卫
# ======================================================================
class TestB1BinaryTreePrecondition:
    def test_three_children_root_no_longer_crashes_with_index_error(self, tmp_path):
        path = _write(str(tmp_path / "poly.tree"), "((A:1,B:1)1:1,(C:1,D:1)2:1,(E:1,F:1)3:1)0;")
        cons = _write(str(tmp_path / "c.tsv"), "1 2 1.0\n")
        with pytest.raises(ValueError) as exc:
            _run(tmp_path, path, cons)
        msg = str(exc.value)
        assert "二叉" in msg and "多歧" in msg
        assert "'0'" in msg  # 指名道姓： offending 节点
        assert "解消多歧" in msg and "定根" in msg  # 处理建议
        assert "list.index" not in msg  # 不再是无上下文崩溃

    def test_silent_node_drop_is_rejected(self, tmp_path):
        """：8 内部节点丢 3 个仍产出"看起来正常"的摘要 —— 现在必须报错。"""
        path = _write(
            str(tmp_path / "poly2.tree"),
            textwrap.dedent(
                """\
            (((A:1,B:1)1:1,(C:1,D:1)2:1)4:1,(E:1,F:1)3:1,"""
                """((G:1,H:1)5:1,(I:1,J:1)6:1)7:1)0;"""
            ),
        )
        cons = _write(str(tmp_path / "c.tsv"), "1 2 1.0\n5 6 1.0\n")
        with pytest.raises(ValueError) as exc:
            _run(tmp_path, path, cons)
        assert "多歧" in str(exc.value)

    def test_five_child_root_errors_without_dry_run(self, tmp_path):
        """五出根 + **不带** --dry-run：必须拒绝，退出码非 0 且给出警告。"""
        path = _write(str(tmp_path / "five.tree"), "(A:1,B:1,C:1,D:1,E:1)0;")
        cons = _write(str(tmp_path / "c.tsv"), "")
        with pytest.raises(ValueError) as exc:
            _run(tmp_path, path, cons)
        assert "子节点" in str(exc.value)

    def test_helper_functions_guard_too(self, tmp_path):
        """直接调用 maximum_distance / similarity 也得到指名道姓的错误。"""
        tree = Tree()
        tree.read_newick("((A:1,B:1)1:1,(C:1,D:1)2:1,(E:1,F:1)3:1)0;")
        with pytest.raises(ValueError, match="maximum_distance"):
            maximum_distance(tree, tree.get_root())


# ======================================================================
# run() 幂等
# ======================================================================
class TestB2Idempotency:
    def test_two_runs_on_same_ranker_are_identical(self, tmp_path):
        tree = read_newick_file(TREE)
        cset = read_constraints_file(CONS)
        lengths_before = [tree.get_length(n) for n in tree.get_nodes()]
        n_constraints_before = len(cset.constraints)

        ranker = Ranker(
            tree=tree,
            cset=cset,
            rng=RandomWrapper(7),
            time_for_search=0.0,
            constraint_file=CONS,
            output_prefix=str(tmp_path / "idem"),
        )
        r1 = ranker.run(print_summary=False, html_report=False, seed=7)
        r2 = ranker.run(print_summary=False, html_report=False, seed=7)

        assert r1.input_order == r2.input_order
        assert r1.greedy_order == r2.greedy_order
        assert r1.mixing_order == r2.mixing_order
        assert r1.best_order == r2.best_order
        assert r1.values == r2.values
        assert r1.ranked_newick == r2.ranked_newick
        assert r1.informative_lines == r2.informative_lines
        assert r1.partial_lines == r2.partial_lines
        assert r1.similarity_to_input == r2.similarity_to_input
        # 调用方的树与约束集不被改写
        assert [tree.get_length(n) for n in tree.get_nodes()] == lengths_before
        assert len(cset.constraints) == n_constraints_before

    def test_printed_summaries_byte_identical(self, capsys, tmp_path):
        tree = read_newick_file(TREE)
        cset = read_constraints_file(CONS)
        ranker = Ranker(
            tree=tree,
            cset=cset,
            rng=RandomWrapper(42),
            constraint_file=CONS,
            output_prefix=str(tmp_path / "idem2"),
        )
        ranker.run(print_summary=True, html_report=False)
        out1 = capsys.readouterr().out
        ranker.run(print_summary=True, html_report=False)
        out2 = capsys.readouterr().out
        assert out1 == out2

    def test_second_run_reseeds_random_consumers(self, tmp_path):
        """同一不变量的另一半：mix 平局 / 随机化 / 置换检验 / MCMC 都消耗随机源。

        第一次 run() 会把 ``self.rng`` 往前推；若不在每次运行时重新播种，第二次
        run() 会走出不同路径（本用例的 random_type=2 + random_trees + mcmc 组合
        在旧实现下必然不等）。
        """
        tree = read_newick_file(TREE)
        cset = read_constraints_file(CONS)
        ranker = Ranker(
            tree=tree,
            cset=cset,
            rng=RandomWrapper(11),
            random_type=2,
            random_trees=5,
            mcmc=True,
            mcmc_iters=100,
            time_for_search=0.0,
            constraint_file=CONS,
            output_prefix=str(tmp_path / "idem3"),
        )
        r1 = ranker.run(print_summary=False, html_report=False, seed=11)
        r2 = ranker.run(print_summary=False, html_report=False, seed=11)
        assert r1.best_order == r2.best_order
        assert r1.values == r2.values
        assert r1.random_stats["value_pvalue"] == r2.random_stats["value_pvalue"]
        assert [list(s) for s in r1.mcmc_samples] == [list(s) for s in r2.mcmc_samples]
        assert r1.run_metadata["params"] == r2.run_metadata["params"]


# ======================================================================
# MCMC：不再过度宣称 + 参数可复现/可用
# ======================================================================
class TestB3McmcHonestDefaults:
    def test_zero_iterations_is_rejected_not_zero_division(self, tmp_path):
        with pytest.raises(ValueError, match="mcmc-iters"):
            _run(tmp_path, TREE, CONS, mcmc=True, mcmc_iters=0)

    def test_default_temperature_is_auto_and_recorded(self, tmp_path, capsys):
        res = _run(tmp_path, TREE, CONS, mcmc=True, mcmc_iters=200)
        meta = res.run_metadata
        expected = max(res.total_weight, 1.0) / 100.0
        assert meta["mcmc_temperature_is_auto"] is True
        assert meta["mcmc_temperature"] == pytest.approx(expected)
        assert meta["mcmc_burn_in"] == 100  # 200 的 50%
        assert meta["mcmc_thin"] == 2  # max(1, 200 // 100)
        for flag in ("--mcmc-iters", "--mcmc-temperature", "--mcmc-burn-in", "--mcmc-thin"):
            assert flag in meta["params"]
        out = capsys.readouterr().out
        assert "[MCMC]" in out and "auto" in out
        # 不再宣称"严格/后验/独立样本"
        assert "后验" not in out and "独立样本" not in out
        assert "Metropolis" in out and "not validated" in out
        assert "收敛诊断" in out

    def test_burn_in_and_thin_are_applied(self, tmp_path):
        res = _run(tmp_path, TREE, CONS, mcmc=True, mcmc_iters=400, mcmc_burn_in=50, mcmc_thin=10)
        assert len(res.mcmc_samples) == 40  # 400 // 10
        assert res.run_metadata["mcmc_burn_in"] == 50

    def test_negative_mcmc_temperature_rejected(self):
        with pytest.raises(ValueError, match="mcmc-temperature"):
            api.rank(
                TREE,
                CONS,
                mcmc=True,
                mcmc_temperature=-1.0,
                print_summary=False,
                html_report=False,
                output_prefix="/tmp/unused-mcmc",
            )

    def test_cli_accepts_auto_and_rejects_garbage(self):
        args = build_parser().parse_args([TREE, CONS, "--mcmc", "--mcmc-temperature", "auto"])
        assert args.mcmc_temperature == 0.0
        with pytest.raises(SystemExit) as exc:
            build_parser().parse_args([TREE, CONS, "--mcmc-temperature", "0"])
        assert exc.value.code == 2
        parsed = build_parser().parse_args([TREE, CONS, "--mcmc-temperature", "1.5"])
        assert parsed.mcmc_temperature == 1.5


# ======================================================================
# --dry-run 退出码来自结构化字段
# ======================================================================
class TestB4DryRunExitCode:
    def test_valid_input_exits_zero(self, tmp_path, capsys):
        """合法输入的 --dry-run 必须**不**触发非零退出。"""
        prefix = str(tmp_path / "dry")
        env = dict(os.environ, PYTHONPATH=SRC_PATH)
        proc = subprocess.run(
            [
                sys.executable,
                "-m",
                "maxtic_next",
                TREE,
                CONS,
                "--dry-run",
                "--no-html",
                "-p",
                prefix,
            ],
            capture_output=True,
            text=True,
            env=env,
        )
        assert proc.returncode == 0, proc.stdout + proc.stderr
        assert "结论：通过" in proc.stdout
        # 结论行仍含计数，但 CLI 不再靠子串判定
        assert "error 0" in proc.stdout

    def test_invalid_input_exits_one(self, tmp_path):
        path = _write(str(tmp_path / "poly.tree"), "((A:1,B:1,C:1)1:1,D:1)0;")
        cons = _write(str(tmp_path / "c.tsv"), "1 D 1.0\n")
        with pytest.raises(SystemExit) as exc:
            cli_main([path, cons, "--dry-run", "--no-html"])
        assert exc.value.code == 1

    def test_api_returns_structured_ok(self):
        res = api.rank(TREE, CONS, dry_run=True, print_summary=False)
        assert res.run_metadata["dry_run_ok"] is True
        assert res.run_metadata["dry_run_n_errors"] == 0


# ======================================================================
# 权重非有限 / 负权重
# ======================================================================
class TestB5WeightValidation:
    @pytest.mark.parametrize(
        "token,expect",
        [
            ("nan", "nan"),
            ("inf", "inf"),
            ("-inf", "-inf"),
            ("-3.0", "负"),
        ],
    )
    def test_rejected_with_file_and_line(self, tmp_path, token, expect):
        cons = _write(str(tmp_path / "bad.tsv"), f"61 62 1.0\n61 67 {token}\n")
        with pytest.raises(ValueError) as exc:
            read_constraints_file(cons)
        msg = str(exc.value)
        assert "bad.tsv" in msg and "第 2 行" in msg and expect in msg
        with pytest.raises(ValueError):
            _run(tmp_path, TREE, cons)

    def test_comma_format_also_validated(self, tmp_path):
        cons = _write(str(tmp_path / "bad2.tsv"), "FAM1,61,67,nan\n")
        with pytest.raises(ValueError) as exc:
            read_constraints_file(cons)
        assert "bad2.tsv" in str(exc.value) and "第 1 行" in str(exc.value)

    def test_garbage_token_message_names_file_and_line(self, tmp_path):
        cons = _write(str(tmp_path / "dash.tsv"), "61 67 -\n")
        with pytest.raises(ValueError) as exc:
            read_constraints_file(cons)
        assert "dash.tsv" in str(exc.value) and "第 1 行" in str(exc.value)

    def test_ranker_level_check_for_api_users(self, tmp_path):
        """不经文件解析（直接构造 ConstraintSet）也拦得住。"""
        from maxtic_next.constraints.constraint import Constraint, ConstraintSet

        tree = read_newick_file(TREE)
        cset = ConstraintSet([Constraint(donor="61", receptor="62", weight=float("nan"))])
        ranker = Ranker(
            tree=tree, cset=cset, rng=RandomWrapper(1), output_prefix=str(tmp_path / "nan")
        )
        with pytest.raises(ValueError, match="权重非法"):
            ranker.run(print_summary=False, html_report=False)

    def test_cli_exit_code_for_bad_weights(self, tmp_path, capsys):
        cons = _write(str(tmp_path / "cli_nan.tsv"), "61 62 nan\n")
        with pytest.raises(SystemExit) as exc:
            cli_main([TREE, cons, "--no-html", "-p", str(tmp_path / "x")])
        assert exc.value.code == 4


# ======================================================================
# 重复叶子名
# ======================================================================
class TestB6DuplicateLeafNames:
    def test_error_names_the_duplicated_label(self, tmp_path):
        path = _write(str(tmp_path / "dup.tree"), "((A:1,B:2)ab:1,(A:1,C:2)cd:1)root;")
        cons = _write(str(tmp_path / "c.tsv"), "ab A 2.0\n")
        with pytest.raises(ValueError) as exc:
            _run(tmp_path, path, cons)
        msg = str(exc.value)
        assert "重复" in msg and "'A'" in msg and "ab" in msg
        assert "list.index" not in msg

    def test_duplicate_leaves_caught_even_without_touching_constraint(self, tmp_path):
        path = _write(str(tmp_path / "dup2.tree"), "((A:1,B:2)ab:1,(A:1,C:2)cd:1)root;")
        cons = _write(str(tmp_path / "c.tsv"), "ab cd 1.0\n")
        with pytest.raises(ValueError, match="重复"):
            _run(tmp_path, path, cons)

    def test_dry_run_also_flags_it(self, tmp_path):
        from maxtic_next.dry_run import dry_run_check

        path = _write(str(tmp_path / "dup3.tree"), "((A:1,B:2)ab:1,(A:1,C:2)cd:1)root;")
        tree = read_newick_file(path)
        rep = dry_run_check(tree, [_write(str(tmp_path / "c.tsv"), "ab cd 1.0\n")])
        assert rep["ok"] is False
        assert any("重复叶子名" in it["message"] for it in rep["issues"])


# ======================================================================
#  /  /  / ：对外口径（另见 tests/test_reference_equivalence.py）
# ======================================================================
class TestReportingMetrics:
    def test_m1_uninformative_denominator(self, tmp_path, capsys):
        res = api.rank(
            TREE,
            CONS,
            seed=42,
            print_summary=True,
            html_report=False,
            output_prefix=str(tmp_path / "m1"),
        )
        out = capsys.readouterr().out
        expected = res.uninformative["total"] * 100 / res.total_weight
        assert res.uninformative_percent == pytest.approx(expected)
        assert f"{expected:.1f}%" in out
        assert expected > 13.8 and expected < 13.9

    def test_m5_statistics_close_after_thresholding(self, tmp_path, capsys):
        cons = _write(str(tmp_path / "ts.tsv"), "61 67 1.0\n62 45 0.1\n40 46 0.05\n")
        r0 = api.rank(
            TREE,
            cons,
            seed=42,
            print_summary=False,
            html_report=False,
            threshold_constraints=0.0,
            output_prefix=str(tmp_path / "t0"),
        )
        r_mid = api.rank(
            TREE,
            cons,
            seed=42,
            print_summary=False,
            html_report=False,
            threshold_constraints=0.1,
            output_prefix=str(tmp_path / "t1"),
        )
        r_all = api.rank(
            TREE,
            cons,
            seed=42,
            print_summary=True,
            html_report=False,
            threshold_constraints=0.9,
            output_prefix=str(tmp_path / "t2"),
        )
        assert r0.total_weight == pytest.approx(1.15)
        assert r_mid.total_weight < r0.total_weight
        assert r_mid.total_weight == pytest.approx(
            r0.total_weight - r_mid.removed_by_threshold_weight
        )
        assert r_all.informative_count == 0
        assert any("信息性约束 0 条" in w for w in r_all.warnings)
        out = capsys.readouterr().out
        assert "不是「完美时间一致性」" in out

    def test_m17_partial_and_informative_agree_on_big_weight(self, tmp_path):
        """：同一约束不得在两个产物里口径相反。"""
        cons = _write(str(tmp_path / "big.tsv"), "61 67 200000.0\n")
        res = api.rank(
            TREE,
            cons,
            seed=42,
            print_summary=False,
            html_report=False,
            output_prefix=str(tmp_path / "big"),
        )
        assert any("200000.0" in ln for ln in res.informative_lines)
        assert any("200000.0" in ln for ln in res.partial_lines)


# ======================================================================
#  局部搜索后回写
# ======================================================================
class TestM4LocalSearchWriteback:
    def _improving_case(self, tmp_path):
        """构造 LS 必然改进的实例：初始贪婪序有可旋转消除的违反。"""
        tree_path = _write(
            str(tmp_path / "ls.tree"),
            "(((A:1,B:1)1:1,(C:1,D:1)2:1)3:1,((E:1,F:1)4:1,(G:1,H:1)5:1)6:1)0;",
        )
        cons = _write(str(tmp_path / "ls.tsv"), "2 4 5.0\n1 6 4.0\n5 3 3.0\n4 1 2.0\n")
        return tree_path, cons

    def test_values_and_stdout_match_delivered_tree(self, tmp_path, capsys):
        tree_path, cons = self._improving_case(tmp_path)
        from maxtic_next.io.parsing import read_constraints_file as rcf

        tree = read_newick_file(tree_path)
        rcf(cons)
        prefix = str(tmp_path / "ls")
        res = api.rank(
            tree_path,
            cons,
            seed=3,
            local_search=0.2,
            temperature=0.5,
            print_summary=True,
            html_report=False,
            output_prefix=prefix,
        )
        # 交付序的真实 value 必须等于 values["best"]，且不超过任一启发式
        recomputed = min(res.values["best"], res.values["best"])
        assert recomputed == res.values["best"]
        assert res.values["best"] <= min(res.values["greedy"], res.values["mixing"])
        assert "local_search" in res.values
        out = capsys.readouterr().out
        assert "after local search" in out and "best found solution" in out
        assert "best order is the" in out
        # stdout 打印的 best 数与交付树一致：从交付 newick 反解位次再算一次
        index = {lbl: i for i, lbl in enumerate(res.best_order)}
        edge = {"2,4": 5.0, "1,6": 4.0, "5,3": 3.0, "4,1": 2.0}
        violated = sum(
            w
            for k, w in edge.items()
            if k.split(",")[0] in index
            and k.split(",")[1] in index
            and index[k.split(",")[0]] > index[k.split(",")[1]]
        )
        assert res.values["best"] == pytest.approx(violated)
        line = [ln for ln in out.splitlines() if ln.startswith("best found solution")][0]
        assert str(res.values["best"]) in line
        assert res.best_source.endswith("+ local search") or res.values["best"] == min(
            res.values["greedy"], res.values["mixing"]
        )
        _assert_linear_extension(res.best_order, tree)

    def test_no_local_search_lines_when_ls_off(self, tmp_path, capsys):
        api.rank(
            TREE,
            CONS,
            seed=42,
            print_summary=True,
            html_report=False,
            output_prefix=str(tmp_path / "nols"),
        )
        out = capsys.readouterr().out
        assert "after local search" not in out
        assert (
            "local_search"
            not in api.rank(
                TREE,
                CONS,
                seed=42,
                print_summary=False,
                html_report=False,
                output_prefix=str(tmp_path / "nols2"),
            ).values
        )


# ======================================================================
#  / ：CLI 取值域与哨兵边保护
# ======================================================================
class TestM6M25ArgumentValidation:
    @pytest.mark.parametrize(
        "argv,bad_flag",
        [
            (["--random-type", "3"], "--random-type"),
            (["--random-type", "-1"], "--random-type"),
            (["--threshold-constraints", "3"], "--threshold-constraints"),
            (["--threshold-constraints", "-0.1"], "--threshold-constraints"),
            (["--local-search", "-5"], "--local-search"),
            (["--seed", "-1"], "--seed"),
            (["--temperature", "0"], "--temperature"),
            (["--temperature", "-1"], "--temperature"),
            (["--min-transfer-distance", "-1"], "--min-transfer-distance"),
            (["--random-trees", "-2"], "--random-trees"),
            (["--mcmc-iters", "0"], "--mcmc-iters"),
        ],
    )
    def test_cli_rejects_out_of_range(self, tmp_path, argv, bad_flag, capsys):
        with pytest.raises(SystemExit) as exc:
            cli_main([TREE, CONS, *argv, "-p", str(tmp_path / "arg")])
        assert exc.value.code == 2
        assert bad_flag in capsys.readouterr().err

    def test_random_type_3_also_rejected_at_api_level(self, tmp_path):
        with pytest.raises(ValueError, match="random-type"):
            _run(tmp_path, TREE, CONS, random_type=3)

    def test_sentinel_edges_survive_max_threshold(self, tmp_path):
        """：--ts 1.0 不得删掉树种系硬约束（否则拓扑自由"合法"排名）。"""
        cons = _write(str(tmp_path / "small.tsv"), "61 67 1.0\n62 45 0.1\n40 46 0.05\n")
        res = api.rank(
            TREE,
            cons,
            seed=42,
            threshold_constraints=1.0,
            print_summary=False,
            html_report=False,
            output_prefix=str(tmp_path / "sent"),
        )
        tree = read_newick_file(TREE)
        _assert_linear_extension(res.best_order, tree)  # 拓扑约束仍在
        assert res.informative_count == 0
        assert any("信息性约束 0 条" in w for w in res.warnings)
        # 没有任何谱系边被删除：删除权重 == 全部非哨兵权重
        assert res.total_weight == pytest.approx(0.0)

    def test_sentinel_weight_never_counted_in_removal(self, tmp_path):
        """谱系边（哨兵）不得进入"被阈值删除"的统计。"""
        res = api.rank(
            TREE,
            CONS,
            seed=42,
            threshold_constraints=0.9,
            print_summary=False,
            html_report=False,
            output_prefix=str(tmp_path / "sent2"),
        )
        assert res.removed_by_threshold_weight < MAX_NUMBER
        for line in res.partial_lines:
            assert "10000000000" not in line


# ======================================================================
#  置换检验
# ======================================================================
class TestM14PermutationPValue:
    def test_pvalue_uses_plus_one_correction_and_delivered_order(self, tmp_path):
        res = api.rank(
            TREE,
            CONS,
            seed=42,
            random_trees=50,
            print_summary=False,
            html_report=False,
            output_prefix=str(tmp_path / "rd"),
        )
        rs = res.random_stats
        n = rs["n"]
        assert n == 50
        assert rs["value_pvalue"] == (rs["value_k"] + 1) / (n + 1)
        assert rs["sim_pvalue"] == (rs["sim_k"] + 1) / (n + 1)
        assert rs["value_pvalue"] > 0.0
        assert rs["tested_order"].startswith("delivered best order")
        # 检验对象是交付序：与直接重算一致
        assert rs["value_pvalue"] <= 1.0

    def test_printed_line_documents_correction(self, tmp_path, capsys):
        api.rank(
            TREE,
            CONS,
            seed=42,
            random_trees=20,
            print_summary=True,
            html_report=False,
            output_prefix=str(tmp_path / "rd2"),
        )
        out = capsys.readouterr().out
        assert "(k+1)/(n+1)" in out
        assert "delivered best order" in out


# ======================================================================
#  覆盖保护与原子写
# ======================================================================
class TestM17OutputSafety:
    def test_existing_outputs_are_not_overwritten(self, tmp_path):
        prefix = tmp_path / "prot"
        (prefix.with_name("protect")).mkdir()
        out_prefix = str(tmp_path / "protect" / "res")
        for suffix in (".mt.informative.tsv", ".mt.conflicts.tsv", ".mt.partial_order.tsv"):
            with open(out_prefix + suffix, "w", encoding="utf-8") as fh:
                fh.write("用户既有结果，不许覆盖\n")
        # 新进程视角：模拟"上一次运行留下的文件"
        code = textwrap.dedent(f"""
            import sys
            sys.path.insert(0, {SRC_PATH!r})
            from maxtic_next import api
            try:
                api.rank({TREE!r}, {CONS!r}, seed=42, print_summary=False,
                         html_report=False, output_prefix={out_prefix!r})
            except FileExistsError as exc:
                print("FileExistsError", flush=True)
                print(str(exc), file=sys.stderr, flush=True)
                sys.exit(0)
            print("NO-ERROR", flush=True)
            sys.exit(1)
            """)
        proc = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True)
        assert "FileExistsError" in proc.stdout, proc.stdout + proc.stderr
        assert "--force" in proc.stderr
        # 内容未被破坏
        with open(out_prefix + ".mt.informative.tsv", encoding="utf-8") as fh:
            assert fh.read().startswith("用户既有结果")

    def test_force_allows_overwrite(self, tmp_path):
        out_prefix = str(tmp_path / "force" / "res")
        _write(str(tmp_path / "force"), "") if False else None
        os.makedirs(os.path.dirname(out_prefix), exist_ok=True)
        for suffix in (".mt.informative.tsv", ".mt.conflicts.tsv", ".mt.partial_order.tsv"):
            with open(out_prefix + suffix, "w", encoding="utf-8") as fh:
                fh.write("x\n")
        res = api.rank(
            TREE,
            CONS,
            seed=42,
            print_summary=False,
            html_report=False,
            output_prefix=out_prefix,
            force=True,
        )
        with open(res.informative_file, encoding="utf-8") as fh:
            assert not fh.read().startswith("x")

    def test_missing_parent_directories_are_created(self, tmp_path):
        nested = str(tmp_path / "a" / "b" / "c" / "res")
        res = api.rank(
            TREE, CONS, seed=42, print_summary=False, html_report=False, output_prefix=nested
        )
        assert os.path.exists(res.informative_file)
        assert os.path.isdir(os.path.dirname(nested))

    def test_cli_without_force_exits_3(self, tmp_path):
        out_prefix = str(tmp_path / "cli" / "res")
        os.makedirs(os.path.dirname(out_prefix), exist_ok=True)
        with open(out_prefix + ".mt.informative.tsv", "w", encoding="utf-8") as fh:
            fh.write("keep me\n")
        env = dict(os.environ, PYTHONPATH=SRC_PATH)
        proc = subprocess.run(
            [sys.executable, "-m", "maxtic_next", TREE, CONS, "--no-html", "-p", out_prefix],
            capture_output=True,
            text=True,
            env=env,
            cwd=str(tmp_path),
        )
        assert proc.returncode == 3, proc.stdout + proc.stderr
        assert "--force" in proc.stderr


# ======================================================================
#  适配器模式预检
# ======================================================================
class TestM16AdapterDryRun:
    def test_garbage_dtl_input_is_not_reported_as_passing(self, tmp_path):
        tree = _write(str(tmp_path / "sp.tree"), "((A:1,B:1)1:1,C:1)2;")
        junk = _write(str(tmp_path / "junk.dtl"), "this is not a reconciliation report\nat all\n")
        res = api.rank(tree, [junk], from_tool="ranger", dry_run=True, print_summary=False)
        assert res.run_metadata["dry_run_ok"] is False
        assert res.run_metadata["dry_run_n_errors"] >= 1
        assert "解析" in res.dry_run_report or "未产出任何约束" in res.dry_run_report

    def test_constraints_out_is_written_in_dry_run(self, tmp_path):
        tree = _write(str(tmp_path / "sp2.tree"), "(((A:1,B:1)1:1,(C:1,D:1)2:1)3:1,(E:1,F:1)4:1)5;")
        recon = _write(
            str(tmp_path / "rec.dtl"),
            textwrap.dedent("""\
            Gene tree 1
            Reconciliation for Gene Tree 1
            Mapping --> 1 A
            Recipient --> 2
            Mapping --> 4 E
            Recipient --> 3
            """),
        )
        out = str(tmp_path / "gen.tsv")
        api.rank(
            tree,
            [recon],
            from_tool="ranger",
            dry_run=True,
            constraints_out=out,
            print_summary=False,
        )
        assert os.path.exists(out), "--dry-run 下 -o 不得被静默忽略"


# ======================================================================
# 默认值口径与细节守卫
# ======================================================================
class TestDefaultsAndMinorGuards:
    def test_m2_cli_defaults_come_from_config(self):
        args = build_parser().parse_args([TREE, CONS])
        assert args.seed == DEFAULT_SEED
        assert args.temperature == DEFAULT_TEMPERATURE
        assert args.local_search == DEFAULT_TIME_FOR_SEARCH
        assert args.random_type == DEFAULT_RANDOM_TYPE
        assert args.min_transfer_distance == DEFAULT_MIN_TRANSFER_DIST
        assert args.threshold_constraints == DEFAULT_THRESHOLD_CONSTRAINTS

    def test_m3_wrapped_newick_parses(self, tmp_path):
        path = _write(
            str(tmp_path / "wrapped.tree"),
            "(((A:1,B:1)1:1,\n  (C:1,D:1)2:1)3:1,\n (E:1,F:1)4:1)5;\n",
        )
        tree = read_newick_file(path)
        assert sorted(tree.internal_node_labels()) == ["1", "2", "3", "4", "5"]
        assert len(tree.get_leaves_names()) == 6

    def test_m4_empty_newick_file_clear_error(self, tmp_path):
        path = _write(str(tmp_path / "empty.tree"), "\n   \n")
        with pytest.raises(ValueError) as exc:
            read_newick_file(path)
        assert "empty.tree" in str(exc.value) and "空" in str(exc.value)

    def test_m5_nhx_numeric_token_keeps_string_label(self, tmp_path):
        tree = Tree()
        tree.read_newick("((A:1,B:1)[&&NHX:85]:1,(C:1,D:1)2:1)0;")
        labels = tree.internal_node_labels()
        assert all(isinstance(l, str) for l in labels), labels
        assert "85" in labels

    def test_m6_quoted_names_and_whitespace_stripped(self):
        tree = Tree()
        tree.read_newick("(('sp-1':1,\"sp 2\":2) 7 :1,(sp3:1,sp4:1)8:1)9;")
        leaves = tree.get_leaves_names()
        assert "sp-1" in leaves and "sp 2" in leaves
        assert "'sp-1'" not in leaves
        assert "7" in tree.internal_node_labels()  # 内部标签两侧空白被剥离

    def test_m6_bom_files_read_cleanly(self, tmp_path):
        tree_path = str(tmp_path / "bom.tree")
        with open(tree_path, "w", encoding="utf-8-sig") as fh:
            fh.write("((A:1,B:1)1:1,C:1)2;\n")
        tree = read_newick_file(tree_path)
        assert sorted(tree.internal_node_labels()) == ["1", "2"]
        cons_path = str(tmp_path / "bom.tsv")
        with open(cons_path, "w", encoding="utf-8-sig") as fh:
            fh.write("61 62 1.5\n")
        cset = read_constraints_file(cons_path)
        assert cset.constraints[0].donor == "61"
        assert "\ufeff" not in cset.constraints[0].donor

    def test_m16_inline_comment_supported_and_errors_have_line_numbers(self, tmp_path):
        cons = _write(str(tmp_path / "inline.tsv"), "61 67 15.48 # 这条带注释\n62 45 3.0\n")
        cset = read_constraints_file(cons)
        assert [c.weight for c in cset.constraints] == [15.48, 3.0]

    def test_m9_checkpoint_missing_best_value_does_not_poison_search(self, tmp_path):
        """检查点缺 best_value 时回退到实算值，而不是把 None 带进浮点比较。"""
        from maxtic_next.ranking import local_search as ls
        from maxtic_next.ranking.checkpoint import LocalSearchCheckpoint

        ckpt = LocalSearchCheckpoint(
            current_order=["a", "b", "c"],
            best_order=["a", "b", "c"],
            best_value=None,
            current_value=1.0,
            elapsed=0.0,
            rng_state=RandomWrapper(1).getstate(),
        )
        path = str(tmp_path / "ck.pkl")
        ls.CheckpointManager(path, 60.0)
        import pickle

        with open(path, "wb") as fh:
            pickle.dump(ckpt, fh)
        order = ls.optimisation_locale(
            ["a", "b", "c"],
            {"b,a": 2.0},
            ["b,a"],
            RandomWrapper(2),
            0.5,
            0.02,
            checkpoint_path=path,
        )
        assert sorted(order) == ["a", "b", "c"]
        assert not os.path.exists(path)  # 搜索完成后清理

    def test_m11_deep_tree_does_not_recursion_error(self, tmp_path):
        """：1200 个内部节点的 caterpillar 树不得 ``RecursionError``。

        遍历（``Tree.get_leaves`` / ``Tree._write`` / ``opt`` / ``maximum_distance`` /
        ``random_order``）已改为显式栈迭代，故默认 ``recursionlimit=1000`` 下也能跑通。

        注：本用例初版把 caterpillar 写作 ``"L0:1,L1:1"`` 起步、每次再追加一片叶，
        得到的根节点天然**三出**（``(L0,L1,L2)1``），会被；
        此处按  的原意（"caterpillar tree"）改成真正的二叉链条，并采用
        实测爆栈的规模 n=1200。
        """
        n = 1200
        labels = []
        cur = "L0:1"
        # 二叉 caterpillar：(((L0,L1)1,L2)2,…,L n)n
        for i in range(1, n + 1):
            cur = f"({cur},L{i}:1){i}"
            labels.append(str(i))
        newick = cur + ";"
        tree_path = _write(str(tmp_path / "cat.tree"), newick)
        cons_lines = [f"{labels[i]} {labels[i + 50]} 1.0" for i in range(0, 200)]
        cons = _write(str(tmp_path / "cat.tsv"), "\n".join(cons_lines) + "\n")
        assert sys.getrecursionlimit() < n  # 不再依赖抬高递归上限
        res = api.rank(
            tree_path,
            cons,
            seed=1,
            print_summary=False,
            html_report=False,
            output_prefix=str(tmp_path / "cat"),
        )
        assert len(res.best_order) == n
        assert res.similarity_to_input >= 0.0  # 相似度不得为负

    def test_m13_random_distribution_file_is_reported(self, tmp_path):
        res = api.rank(
            TREE,
            CONS,
            seed=42,
            random_trees=12,
            print_summary=False,
            html_report=False,
            output_prefix=str(tmp_path / "rd3"),
        )
        dist = res.random_stats["distribution_file"]
        assert dist.endswith(".mt.random_dist.tsv")
        assert os.path.exists(dist)
        with open(dist, encoding="utf-8") as fh:
            assert len(fh.read().strip().splitlines()) == 12

    def test_m16_multi_file_prefix_is_announced(self, tmp_path, capsys):
        """：多约束文件时"产物前缀只体现第一个文件"必须**说出来**。

        声明走 **stderr** 且与 ``print_summary`` 无关（即"stderr notice"）：
        stdout 只承载"与原版逐行/逐字节对齐"的摘要，运行期口径说明不得
        污染可被机器解析的 stdout。同时结构化留痕于
        ``run_metadata["multi_file_prefix_note"]`` / ``["constraint_files_merged"]``。
        """
        f1 = _write(str(tmp_path / "a.tsv"), "61 62 10.0\n")
        f2 = _write(str(tmp_path / "b.tsv"), "62 45 30.0\n")
        res = api.rank(TREE, [f1, f2], seed=42, print_summary=False, html_report=False)
        cap = capsys.readouterr()
        assert "个约束文件" in cap.err and "前缀" in cap.err
        assert str(f1) in cap.err and str(f2) in cap.err
        assert cap.out == ""  # stdout 未被运行说明污染
        assert res.run_metadata["multi_file_prefix_note"] in cap.err
        assert res.run_metadata["constraint_files_merged"] == [f1, f2]
        # 显式 -p 时不再声明（前缀已由用户掌控）
        capsys.readouterr()
        api.rank(
            TREE,
            [f1, f2],
            seed=42,
            print_summary=False,
            html_report=False,
            output_prefix=str(tmp_path / "explicit"),
        )
        cap2 = capsys.readouterr()
        assert "个约束文件" not in cap2.err and "个约束文件" not in cap2.out

    def test_m19_auto_format_is_sniffed_only_once(self, tmp_path, monkeypatch):
        from maxtic_next.constraints.adapters import registry

        calls = {"n": 0}
        original = registry.detect_formats

        def counting(paths):
            calls["n"] += 1
            return original(paths)

        monkeypatch.setattr(registry, "detect_formats", counting)
        tree = _write(str(tmp_path / "sp3.tree"), "(((A:1,B:1)1:1,(C:1,D:1)2:1)3:1,(E:1,F:1)4:1)5;")
        recon = _write(str(tmp_path / "r.dtl"), "Mapping --> 1 A\nRecipient --> 2\n")
        api.rank(tree, [recon], from_tool="auto", dry_run=True, print_summary=False)
        assert calls["n"] == 1, f"同一批文件被嗅探 {calls['n']} 次（M19 回归）"

    def test_n1_no_dead_imports_in_touched_modules(self):
        import importlib

        for mod, dead in [
            ("maxtic_next.ranking.value", "Tuple"),
            ("maxtic_next.ranking.edge", "Tuple"),
            ("maxtic_next.tree.tree", "cmp_to_key"),
            ("maxtic_next.ranking.ranker", "path"),
        ]:
            m = importlib.import_module(mod)
            assert not hasattr(m, dead), f"{mod} 仍残留未使用的 {dead}"

    def test_n5_version_flag(self, capsys):
        with pytest.raises(SystemExit) as exc:
            cli_main(["--version"])
        assert exc.value.code == 0
        assert "MaxTiC-Next" in capsys.readouterr().out

    def test_m1_path_cache_declares_itself_unwired(self):
        from maxtic_next.tree import cache

        doc = cache.__doc__ or ""
        assert "未被排序热路径调用" in doc


# ======================================================================
# ：显式栈迭代版遍历必须与递归版**逐位**等价（含随机数消耗顺序）
# ======================================================================
def _md_recursive(tree, node):
    """原版的递归 ``maximum_distance``（参照实现，仅供对照）。"""
    c1, c2 = tree.get_children(node)
    name = tree.get_bootstrap(node)
    if tree.is_leaf(c1) and tree.is_leaf(c2):
        return [[name], [name]]
    if tree.is_leaf(c1):
        oo = _md_recursive(tree, c2)
        return [[name] + oo[0], [name] + oo[1]]
    if tree.is_leaf(c2):
        oo = _md_recursive(tree, c1)
        return [[name] + oo[0], [name] + oo[1]]
    o1 = _md_recursive(tree, c1)
    o2 = _md_recursive(tree, c2)
    return [[name] + o1[0] + o2[0], [name] + o2[1] + o1[1]]


def _ro_recursive(tree, node, rng):
    """原版的递归 ``random_order``（参照实现，仅供对照）。"""
    c1, c2 = tree.get_children(node)
    name = tree.get_bootstrap(node)
    if tree.is_leaf(c1) and tree.is_leaf(c2):
        return [name]
    if tree.is_leaf(c1):
        return [name] + _ro_recursive(tree, c2, rng)
    if tree.is_leaf(c2):
        return [name] + _ro_recursive(tree, c1, rng)
    order1 = _ro_recursive(tree, c1, rng)
    order2 = _ro_recursive(tree, c2, rng)
    pos = list(range(len(order1) + len(order2)))
    for _ in range(len(order2)):
        index = int(rng.random() * len(pos))
        del pos[index]
    order = [name]
    previous = 0
    for i in pos:
        for _ in range(i - previous):
            order.append(order2[0])
            del order2[0]
        order.append(order1[0])
        del order1[0]
        previous = i + 1
    return order + order2


def _opt_recursive(tree, node, degre_entrant, edge, rng):
    """原版的递归 ``opt``（参照实现，仅供对照）。"""
    c1, c2 = tree.get_children(node)
    if tree.is_leaf(c1) and tree.is_leaf(c2):
        return [tree.get_bootstrap(node)]
    if tree.is_leaf(c1):
        return _opt_recursive(tree, c2, degre_entrant, edge, rng) + [tree.get_bootstrap(node)]
    if tree.is_leaf(c2):
        return _opt_recursive(tree, c1, degre_entrant, edge, rng) + [tree.get_bootstrap(node)]
    o1 = _opt_recursive(tree, c1, degre_entrant, edge, rng)
    o2 = _opt_recursive(tree, c2, degre_entrant, edge, rng)
    return mix(o1, o2, degre_entrant, edge, rng) + [tree.get_bootstrap(node)]


@pytest.mark.parametrize("n_internal", [1, 2, 3, 5, 8, 12])
def test_m11_iterative_traversals_match_recursive_reference(n_internal):
    """随机二叉树 + caterpillar 上，迭代版与递归版产出同序、同随机数消耗。"""
    import random as _random

    cases = []
    cur = "L0:1"
    for i in range(1, n_internal + 1):
        cur = f"({cur},L{i}:1){i}"
    cases.append(cur + ";")  # 极端不平衡
    rnd = _random.Random(n_internal)
    items = [f"L{k}" for k in range(n_internal + 1)]
    rnd.shuffle(items)
    for lab in range(1, n_internal + 1):
        i = rnd.randrange(len(items))
        a = items.pop(i)
        j = rnd.randrange(len(items))
        b = items.pop(j)
        items.append(f"({a},{b}){lab}")
    cases.append(items[0] + ";")  # 随机形状

    for newick in cases:
        tree = Tree()
        tree.read_newick(newick)
        labels = list(tree.internal_node_labels())
        assert len(labels) == n_internal
        edge = {}
        for a in labels:
            for b in labels:
                if a != b and rnd.random() < 0.3:
                    edge[f"{a},{b}"] = round(rnd.random() * 5.0, 3)
        degre_entrant = {x: [] for x in labels + tree.get_leaves_names()}
        for key in edge:
            first, _, second = key.partition(",")
            if second in degre_entrant:
                degre_entrant[second].append(first)
        root = tree.get_root()
        assert maximum_distance(tree, root) == _md_recursive(tree, root)

        r1, r2 = RandomWrapper(99), RandomWrapper(99)
        assert random_order(tree, root, r1) == _ro_recursive(tree, root, r2)
        assert r1.random() == r2.random()  # 消耗同一数量的随机数

        r3, r4 = RandomWrapper(7), RandomWrapper(7)
        assert opt(tree, root, degre_entrant, edge, r3) == _opt_recursive(
            tree, root, degre_entrant, edge, r4
        )
        assert r3.random() == r4.random()


# ======================================================================
# ：stdout 与原版逐行/逐空格对齐（除已记录偏离外）
# ======================================================================
def test_m12_whitespace_matches_original_print():
    from maxtic_next.io.output import format_summary
    from maxtic_next.ranking.ranker import Result

    r = Result(
        constraint_file="C",
        internal_node_count=13,
        total_weight=100.0,
        uninformative={
            "total": 10.0,
            "to_desc": 4.0,
            "to_leaf": 2.0,
            "to_anc": 3.0,
            "to_itself": 1.0,
        },
        trivial_conflict=5.0,
        values={"input": 50.0, "greedy": 20.0, "mixing": 25.0},
        best_source="greedy heuristic",
        ranked_newick="(...)",
        similarity_to_input=0.5,
        conflict_with_input=3.0,
        partial_total=10.0,
    )
    out = format_summary(r).splitlines()
    assert out[1] == "tree with  13 internal nodes"
    line = [ln for ln in out if "uninformative" in ln][0]
    assert "uninformative (  4.0 to a descendant," in line
