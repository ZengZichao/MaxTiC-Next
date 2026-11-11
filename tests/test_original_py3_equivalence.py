"""Cross-version check against a Python 3 (2to3) translation of the original MaxTiC.

The original ``MaxTiC.py`` is a Python 2 script and cannot run under Python 3.
This test runs the vendored 2to3 translation in
``tests/data/original_maxitc/`` (see its README) on the two instances used in the
article and compares the three output tables and the objective values with
``maxtic_next``.  It is a cross-check of the original control flow, not a
substitute for running the original under Python 2.7.
"""

import os
import shutil
import subprocess
import sys

import pytest

from maxtic_next import api

HERE = os.path.dirname(os.path.abspath(__file__))
ORIG = os.path.join(HERE, "data", "original_maxitc")
MAIN = os.path.join(HERE, "data", "main_synthetic")
CYANO_TREE = os.path.join(HERE, "data", "minitree.tree")
CYANO_CONS = os.path.join(HERE, "data", "Cyano_CUTConstraints.tsv")

ORIG_SUFFIXES = (
    "_MT_output_filtered_list_of_weighted_informative_constraints",
    "_MT_output_list_of_constraints_conflicting_with_best_order",
    "_MT_output_partial_order",
)

pytestmark = pytest.mark.skipif(
    not os.path.isfile(os.path.join(ORIG, "MaxTiC.py")),
    reason="vendored 2to3 translation of the original MaxTiC is absent",
)


def _run_original(tmp_path, tree, cons):
    """Run the 2to3-translated original and return (stdout, three output paths)."""
    work = tmp_path / "original"
    work.mkdir()
    tree_name = "tree.nwk"
    cons_name = "constraints.tsv"
    shutil.copy(tree, work / tree_name)
    shutil.copy(cons, work / cons_name)
    env = dict(os.environ)
    env["PYTHONPATH"] = ORIG + os.pathsep + env.get("PYTHONPATH", "")
    proc = subprocess.run(
        [sys.executable, os.path.join(ORIG, "MaxTiC.py"), tree_name, cons_name],
        cwd=work,
        env=env,
        capture_output=True,
        text=True,
        timeout=300,
    )
    assert proc.returncode == 0, proc.stderr
    outs = [str(work / (cons_name + suffix)) for suffix in ORIG_SUFFIXES]
    for path in outs:
        assert os.path.isfile(path), path
    return proc.stdout, outs


def _parse_values(stdout):
    """Extract the input/greedy/mixing objective values from the original stdout."""
    values = {}
    for line in stdout.splitlines():
        if line.startswith("value of the order given by the input tree"):
            values["input"] = float(line.split()[9])
        elif line.startswith("value of the greedy heuristic"):
            values["greedy"] = float(line.split()[5])
        elif line.startswith("value of the mixing heuristic:"):
            values["mixing"] = float(line.split()[5])
    return values


def _compare(tmp_path, tree, cons):
    stdout, orig_outs = _run_original(tmp_path, tree, cons)
    ref = api.rank(
        tree,
        cons,
        seed=42,
        print_summary=False,
        html_report=False,
        output_prefix=str(tmp_path / "next"),
    )
    new_outs = [ref.informative_file, ref.conflicting_file, ref.partial_order_file]

    # The three tables must be identical (line order included: both write in the
    # same descending-weight order).  If a future change makes the order differ,
    # compare sorted lines instead and investigate before relaxing this.
    for orig, new in zip(orig_outs, new_outs):
        with open(orig, encoding="utf-8") as fh_o, open(new, encoding="utf-8") as fh_n:
            assert fh_o.read() == fh_n.read(), (orig, new)

    orig_values = _parse_values(stdout)
    assert orig_values["input"] == pytest.approx(ref.values["input"], abs=1e-9)
    assert orig_values["greedy"] == pytest.approx(ref.values["greedy"], abs=1e-9)
    assert orig_values["mixing"] == pytest.approx(ref.values["mixing"], abs=1e-9)


def test_original_py3_matches_on_main_synthetic(tmp_path):
    _compare(tmp_path, os.path.join(MAIN, "tree.nwk"), os.path.join(MAIN, "constraints.tsv"))


def test_original_py3_matches_on_cyanobacteria(tmp_path):
    _compare(tmp_path, CYANO_TREE, CYANO_CONS)
