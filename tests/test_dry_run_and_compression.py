"""预检、压缩输入与元数据完整性的回归测试。

覆盖的不变量：

* ``--dry-run`` 的适配器分支里，异常路径必须真的让预检失败——打印 ``[error]``
  却仍判"通过"并以退出码 0 结束是不可接受的（计数与 ``ok`` 必须在**两条**路径上重算）。
* ``--min-endpoint-hit-rate 0`` 只关闭低命中率**告警**，绝不关闭"端点全部未命中"
  的硬错误，否则与 ``--help`` 的承诺相反。
* ``FRQ`` 汇总行按**首字段**判定；按整行子串匹配会把类群名恰含该子串的数据行静默删除。
* 适配器对自环转移写 ``distance = None``，不能写 ``0.0``——否则
  ``--min-transfer-distance`` 会连"供体 = 受体谱系"的事件一起吃掉。
* 上游输出常以 ``.gz`` / ``.tgz`` 分发，必须透明解包；多成员归档给出可执行提示。
* 近优解收集容量 ``top_k`` 由 API/CLI 参数化，不写死。
* 内核提供协作式停止钩子 ``stop_check``，取消在迭代边界真正生效。
* ``--random-trees`` 的进度反馈写 stderr，长检验期间不是零输出。
* ``run_metadata`` 不得有重复键（重复字面量键会静默丢数据）。

约定与其余测试一致：全部产物写入 ``tmp_path``，不改动仓库文件。

``maxtic-studio`` 入口点属图形前端，其回归测试在 Studio 仓库的
``tests/test_params_validation.py``。
"""

from __future__ import annotations

import gzip
import io as _io
import os
import tarfile
import time

import pytest

from maxtic_next import api
from maxtic_next.constraints.constraint import Constraint, ConstraintSet
from maxtic_next.constraints.adapters._diagnostics import (
    LabelMismatchError,
    check_endpoint_hit_rate,
)
from maxtic_next.constraints.parsers import parse_constraints
from maxtic_next.io.compression import CompressedArchiveError, read_text_lines
from maxtic_next.ranking.ranker import check_near_optimal_top_k

HERE = os.path.dirname(os.path.abspath(__file__))
TREE = os.path.join(HERE, "data", "minitree.tree")
CONS = os.path.join(HERE, "data", "Cyano_CUTConstraints.tsv")
ECC_SPECIES = os.path.join(HERE, "data", "eccetera_example", "species.tree")
ECC_REPORT = os.path.join(HERE, "data", "eccetera_example", "FAM1.recphyloxml")


def _run(tmp_path, name, constraints=TREE and CONS, **kwargs):
    kwargs.setdefault("print_summary", False)
    kwargs.setdefault("html_report", False)
    kwargs.setdefault("seed", 42)
    kwargs["output_prefix"] = str(tmp_path / name)
    return api.rank(TREE, constraints, force=True, **kwargs)


# ----------------------------------------------------------------------
# 适配器分支的错误必须真的让预检失败
# ----------------------------------------------------------------------
def test_dry_run_adapter_error_path_is_not_reported_as_pass(tmp_path):
    broken = tmp_path / "truncated.recphyloxml"
    with open(ECC_REPORT, encoding="utf-8") as fh:
        text = fh.read()
    broken.write_text(text[: len(text) // 2], encoding="utf-8")

    res = api.rank(
        ECC_SPECIES,
        [str(broken)],
        from_tool="eccetera",
        dry_run=True,
        print_summary=False,
        output_prefix=str(tmp_path / "dry"),
    )
    meta = res.run_metadata
    assert meta["dry_run"] is True
    # 若只在正常路径重算，这里会打印 [error] 却 dry_run_ok=True、n_errors=0（退出码 0）
    assert meta["dry_run_n_errors"] >= 1
    assert meta["dry_run_ok"] is False


# ----------------------------------------------------------------------
#  残留：阈值为 0 只关闭告警，0 命中仍是硬错误
# ----------------------------------------------------------------------
def test_zero_hit_rate_still_raises_when_guard_disabled():
    diag = {
        "transfers_seen": 5,
        "transfers_resolved": 0,
        "unresolved_endpoint_samples": ["3", "5"],
        "warnings": [],
    }
    with pytest.raises(LabelMismatchError):
        check_endpoint_hit_rate(diag, species_labels={"A", "B"}, tool="eccetera", threshold=0.0)


def test_threshold_zero_does_not_warn_on_partial_hit():
    diag = {
        "transfers_seen": 4,
        "transfers_resolved": 2,
        "unresolved_endpoint_samples": ["7"],
        "warnings": [],
    }
    rate = check_endpoint_hit_rate(diag, species_labels={"A", "B"}, tool="ranger", threshold=0.0)
    assert rate == pytest.approx(0.5)
    assert diag["warnings"] == []


# ----------------------------------------------------------------------
# FRQ 只在作为首字段时才是汇总行
# ----------------------------------------------------------------------
def test_taxon_name_containing_frq_is_not_dropped():
    cset = parse_constraints(["FRQX 2 5.0", "1 2 1.0", "XFRQ 3 2.0"])
    keys = {(c.donor, c.receptor) for c in cset.constraints}
    assert ("FRQX", "2") in keys
    assert ("XFRQ", "3") in keys
    assert len(cset.constraints) == 3


def test_real_frq_summary_line_is_skipped():
    assert len(parse_constraints(["FRQ gene1 61 67 0.5"]).constraints) == 0


# ----------------------------------------------------------------------
# ：自环事件的 distance 必须是 None（无距离），不是 0.0
# ----------------------------------------------------------------------
def test_self_loop_without_distance_survives_filter():
    cset = ConstraintSet(
        [
            Constraint(donor="X", receptor="X", weight=1.0, metadata={"distance": None}),
            Constraint(donor="Y", receptor="Y", weight=2.0, metadata={"distance": 0.0}),
        ]
    )
    report = cset.filter_by_distance(0.0)
    assert len(cset.constraints) == 1
    assert cset.constraints[0].donor == "X"
    assert report["self_loops_dropped"] == 1


# ----------------------------------------------------------------------
# ：gzip / tar.gz 输入
# ----------------------------------------------------------------------
def test_gzip_inputs_produce_identical_results(tmp_path):
    gz_tree = tmp_path / "tree.nwk.gz"
    gz_cons = tmp_path / "cons.tsv.gz"
    for src, dst in ((TREE, gz_tree), (CONS, gz_cons)):
        with open(src, "rb") as fin, gzip.open(dst, "wb") as fout:
            fout.write(fin.read())

    plain = _run(tmp_path, "plain")
    packed = _run(tmp_path, "gz", constraints=str(gz_cons))
    # 同一约束内容：结果必须逐字节一致（树不同则用 plain 的树）
    res_gz_tree = api.rank(
        str(gz_tree),
        str(gz_cons),
        force=True,
        print_summary=False,
        html_report=False,
        seed=42,
        output_prefix=str(tmp_path / "gz2"),
    )
    assert res_gz_tree.best_order == plain.best_order
    assert res_gz_tree.values == plain.values
    assert packed.best_order == plain.best_order


def test_multi_member_tar_gives_actionable_error(tmp_path):
    arc = tmp_path / "many.tar.gz"
    with tarfile.open(arc, "w:gz") as tar:
        for name, payload in (("a.txt", b"1 2 1.0\n"), ("b.txt", b"2 3 1.0\n")):
            raw = payload
            info = tarfile.TarInfo(name)
            info.size = len(raw)
            tar.addfile(info, _io.BytesIO(raw))
    with pytest.raises(CompressedArchiveError) as exc:
        read_text_lines(str(arc))
    message = str(exc.value)
    assert "tar" in message and "a.txt" in message


def test_single_member_tar_is_read_transparently(tmp_path):
    member = tmp_path / "cons.tsv"
    with open(CONS, "rb") as fin:
        payload = fin.read()
    member.write_bytes(payload)
    arc = tmp_path / "one.tar.gz"
    with tarfile.open(arc, "w:gz") as tar:
        tar.add(str(member), arcname="cons.tsv")
    plain_lines = [ln.rstrip("\n") for ln in open(CONS, encoding="utf-8").read().splitlines()]
    got = [ln.rstrip("\n") for ln in read_text_lines(str(arc))]
    assert got == plain_lines


# ----------------------------------------------------------------------
#  残留：top_k 三层可覆盖并被记录
# ----------------------------------------------------------------------
@pytest.mark.parametrize("bad", [0, -1, "many", 1.5, True, float("nan")])
def test_top_k_rejects_invalid_values(bad):
    with pytest.raises(ValueError):
        check_near_optimal_top_k(bad)


def test_top_k_is_threaded_and_recorded(tmp_path):
    res = _run(tmp_path, "topk", local_search=0.2, temperature=0.001, top_k=300)
    assert res.run_metadata["near_optimal_top_k"] == 300


def test_cli_exposes_near_optimal_top_k():
    from maxtic_next import cli

    parser = cli.build_parser()
    args = parser.parse_args([TREE, CONS, "--near-optimal-top-k", "250"])
    assert args.near_optimal_top_k == 250
    with pytest.raises(SystemExit):
        parser.parse_args([TREE, CONS, "--near-optimal-top-k", "0"])


# ----------------------------------------------------------------------
#  残留：取消必须真的打断长时长局部搜索
# ----------------------------------------------------------------------
def test_stop_check_interrupts_local_search(tmp_path):
    started = time.time()
    res = _run(tmp_path, "cancel", local_search=30.0, temperature=0.001, stop_check=lambda: True)
    elapsed = time.time() - started
    assert elapsed < 10.0, f"取消未生效：耗时 {elapsed:.1f}s"
    assert res.best_order


def test_mcmc_stop_check_reaches_the_chain(tmp_path):
    """取消钩子必须抵达 MCMC 链：否则 10 万步会一路跑完。

    链在产出第一个样本之前被取消时，本实现按"没有可交付样本"如实报错
    （而不是静默返回空列表让上层以为采到了样本），因此错误信息本身就是
    钩子已接通的证据。
    """
    started = time.time()
    with pytest.raises(RuntimeError, match="取消"):
        _run(tmp_path, "mcmcstop", mcmc=True, mcmc_iters=100000, stop_check=lambda: True)
    assert time.time() - started < 10.0


# ----------------------------------------------------------------------
# ：--random-trees 的进度走 stderr，stdout 仍是可机读摘要
# ----------------------------------------------------------------------
def test_random_trees_progress_goes_to_stderr(tmp_path, capsys):
    api.rank(
        TREE,
        CONS,
        random_trees=5,
        force=True,
        print_summary=True,
        html_report=False,
        output_prefix=str(tmp_path / "rd"),
    )
    captured = capsys.readouterr()
    assert "random" in captured.err.lower()
    assert "pvalue" in captured.out


# ----------------------------------------------------------------------
# run_metadata 不得有重复键（重复字面量键会静默丢数据）
# ----------------------------------------------------------------------
def test_no_duplicate_literal_keys_in_source():
    """字面量重复键会静默丢弃前一个值 —— ``run_metadata`` 有多个写入点，必须看守。

    用 AST 扫全源码，任何字典字面量里出现同一字符串键两次即失败。
    """
    import ast

    root = os.path.join(os.path.dirname(api.__file__))
    offenders = []
    for dirpath, dirnames, filenames in os.walk(root):
        dirnames[:] = [d for d in dirnames if d != "__pycache__"]
        for fn in filenames:
            if not fn.endswith(".py"):
                continue
            path = os.path.join(dirpath, fn)
            tree = ast.parse(open(path, encoding="utf-8").read(), str(path))
            for node in ast.walk(tree):
                if not isinstance(node, ast.Dict):
                    continue
                seen = set()
                for key in node.keys:
                    if isinstance(key, ast.Constant) and isinstance(key.value, str):
                        if key.value in seen:
                            offenders.append(f"{path}:{node.lineno}:{key.value}")
                        seen.add(key.value)
    assert offenders == [], f"字典字面量存在重复键：{offenders}"


def test_run_metadata_records_mcmc_wiring(tmp_path):
    res = _run(tmp_path, "meta", mcmc=True, mcmc_iters=200)
    meta = res.run_metadata
    params = meta.get("params")
    assert "mcmc" in str(params), str(params)[:200]
    assert "--near-optimal-top-k" in str(params)
    assert meta.get("mcmc_chain_diagnostics") or "mcmc_diagnostics" in str(params)
