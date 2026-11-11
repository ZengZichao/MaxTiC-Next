# Real ARTra reconciliation report (adapter validation)

- `output.txt`: the reconciliation report distributed with ARTra (Mondal et al.;
  reference [9] of the manuscript).
- `species.nwk`: the species tree, taken from the first line of ARTra's
  `input.newick` for the same run.

The report describes a 199-taxon species tree and one gene tree with 13
transfers. It is used by `tests/test_adapter_registry.py` to validate the ARTra
adapter on a real report (with `min_support=0`; the default 0.05 would drop all
13 one-count transfers at 1/24 = 0.042). The files are redistributed unchanged
for testing under the terms of the original ARTra distribution.
