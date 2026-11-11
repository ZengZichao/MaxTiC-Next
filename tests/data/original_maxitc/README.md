# Original MaxTiC — Python 3 (2to3) translation for cross-version tests

These two files are a **2to3 translation** of the original MaxTiC scripts:

- `MaxTiC.py`
- `script_tree.py`

from the original MaxTiC distribution by Eric Tannier (Inria), distributed
under the CeCILL licence. The original header and citation are retained in the
files. The translation is used only by `tests/test_original_py3_equivalence.py`
to run the original control flow under Python 3 in environments where a Python 2
interpreter is unavailable.

What was changed relative to the Python 2 source:

1. `2to3` converted `print` statements, `has_key`, `iteritems`, `xrange`, etc.
2. Tabs were expanded to 8 spaces so Python 3 accepts the indentation.
3. `sort(lambda x, y: cmp(...))` was replaced by
   `sort(key=cmp_to_key(lambda x, y: cmp(...)))`.
4. A two-line `cmp(a, b)` helper was added (Python 3 removed the built-in).

No algorithmic change was made. This is **not** a substitute for running the
original under Python 2.7; it is a reproducible cross-check of the original
control flow. The true Python 2 binary comparison remains future work.
