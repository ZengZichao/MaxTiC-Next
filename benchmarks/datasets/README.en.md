# MaxTiC-Next Real-Dataset Benchmark Directory

> 🌐 中文版：[README.md](README.md)

This directory holds real large-scale HGT (horizontal gene transfer) datasets for
performance stress testing via the `--real` mode of the benchmark suite module
`src/maxtic_next/benchmarks/suite.py` (entry point
`python -m maxtic_next.benchmarks.suite --real`).

## Directory Structure

```
benchmarks/datasets/
├── README.md            ← Chinese version
├── README.en.md         ← this file
├── <dataset_name>/
│   ├── <dataset_name>.tree          ← species tree (Newick format)
│   └── <dataset_name>.tsv           ← constraint file (TSV / plain-text dual format)
└── ...
```

> Directory scanning is **not recursive**, so `--data-dir` should point at the level that
> directly contains the tree and the constraints (e.g.
> `benchmarks/datasets/<dataset_name>/`). Scanning also skips MaxTiC's own output files
> (the three output files / distribution / HTML report suffixes, see
> `suite._is_maxtic_output`), so that ranking output is not re-ingested as input
> constraints.

## How to Add a Real Dataset

### 1. Prepare the Data Files

Each dataset consists of two files:

- **Species tree file**: Newick format, with extension `.tree` / `.nwk` / `.newick` / `.tre`
- **Constraint file**: MaxTiC constraint format (three columns: `donor\treceptor\tweight`),
  with extension `.tsv` / `.txt` / `.csv` / `.constraints`

### 2. Naming Convention

To enable automatic pairing, it is recommended that the tree file and the constraint
file share the same filename prefix (identical after stripping the extension):

```
my_dataset/
├── cyano.tree              ← species tree
└── cyano.tsv               ← constraints
```

If the directory contains only one tree file and one constraint file, they will be
paired automatically even if the names differ.

### 3. Running the Benchmark

> **This repository does not ship real datasets** (they are too large):
> `benchmarks/datasets/` currently holds only the README files, so you have to place
> your own data following the conventions above. The **default value of `--data-dir`
> is the `examples/` directory at the repository root**, which contains only the small
> demo pair `minitree.tree` + `Cyano_CUTConstraints.tsv`.

```bash
# Use the data in the default examples/ directory
python -m maxtic_next.benchmarks.suite --real

# Specify a custom dataset directory (replace the <path/to/datasets> placeholder with
# the directory where you actually put your data)
python -m maxtic_next.benchmarks.suite --real --data-dir <path/to/datasets>

# Specify the seed and the time budget
python -m maxtic_next.benchmarks.suite --real --data-dir <path/to/datasets> --seed 42 --time-budget 120
```

If the directory given to `--data-dir` contains no tree + constraints pair (e.g. the
still-empty `benchmarks/datasets/`), `suite.py` prints the Chinese warning
"警告：目录 <data_dir> 中未找到 tree+constraints 数据对。" verbatim (English gloss:
"Warning: no tree+constraints data pair found in directory <data_dir>") and skips all
datasets.

### 4. Output Metrics

For each dataset, the benchmark reports the following key metrics:

| Metric | Description |
|--------|-------------|
| `n_internal` | Number of internal nodes in the species tree |
| `n_constraints` | Number of valid lines in the constraint file (constraint count, excluding `#` comment lines and FRQ header lines) |
| `n_distinct_pairs` | Number of distinct `(donor, receptor)` pairs after deduplication |
| `cons_*` | Honest size breakdown of the constraint file: `cons_n_lines_nonempty` / `cons_n_constraints_parsed` / `cons_n_distinct_pairs` / `cons_n_comment_or_header_lines` |
| `cpu_time` | CPU time (seconds, median of `time.process_time`; wall-clock reported separately as `wall_time`) |
| `cpu_time_iqr` / `wall_time_iqr` | Inter-quartile range of `cpu_time` / `wall_time` |
| `wall_time` | Wall-clock time (seconds, median of `time.perf_counter`) |
| `io_time` / `kernel_cpu_time` | I/O cost of writing the output files (CPU median) / ranking-kernel CPU time after subtracting I/O |
| `peak_memory_mb` | Peak memory (MiB, `tracemalloc`; computed as `peak / (1024 * 1024)`, the key name in the result dict is still `peak_memory_mb`) |
| `best_value` | Objective value of the best ranking (the smaller of greedy / mixing) |
| `delivered_value` | Independently recomputed objective value of the delivered ranking |
| `best_source` | Provenance of the delivered ranking (`greedy heuristic` / `mixing heuristic`, possibly suffixed `+ local search`) |
| `correctness_ok` | Whether the recomputed delivered value passes the correctness assertions |
| `similarity_to_input` | Kendall similarity between the best ranking and the input tree |
| `exceeded` | Whether the time budget was exceeded |

## Data Format Details

### Species Tree (Newick)

Standard Newick format, with internal nodes labeled by bootstrap values:

```
(((((CYAP8:1,CYAP0:1)42:7,(CYAA5:7,UCYNA:7)41:1)59:2,...
```

### Constraint File (TSV)

One constraint per line, three tab-separated columns: `donor\treceptor\tweight`

```
61	62	15.04
61	67	15.48
61	46	28.73
```

- `donor`: transfer donor node (internal node label)
- `receptor`: transfer receptor node (internal node label)
- `weight`: constraint weight (floating point, usually the transfer support)

## Notes

- Large-scale datasets (more than 500 internal nodes) may take considerable time;
  set `--time-budget` appropriately.
- The benchmark measures peak memory with `tracemalloc`, which tracks only memory
  allocated by Python and excludes direct allocations made by C extensions.
- All output files (informative / conflicting / partial order) are written to a
  temporary directory by default and do not pollute the workspace.
