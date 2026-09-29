# MaxTiC-Next Real-Dataset Benchmark Directory

> 🌐 中文版：[README.md](README.md)

This directory holds real large-scale HGT (horizontal gene transfer) datasets for
performance stress testing via the `--real` mode of `benchmarks/suite.py`.

## Directory Structure

```
benchmarks/datasets/
├── README.md            ← this file
├── <dataset_name>/
│   ├── <dataset_name>.tree          ← species tree (Newick format)
│   └── <dataset_name>.tsv           ← constraint file (TSV / plain-text dual format)
└── ...
```

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

```bash
# Use the data in the default examples/ directory
python -m maxtic_next.benchmarks.suite --real

# Specify a custom dataset directory
python -m maxtic_next.benchmarks.suite --real --data-dir benchmarks/datasets/

# Specify the seed and the time budget
python -m maxtic_next.benchmarks.suite --real --data-dir benchmarks/datasets/ --seed 42 --time-budget 120
```

### 4. Output Metrics

For each dataset, the benchmark reports the following metrics:

| Metric | Description |
|--------|-------------|
| `n_internal` | Number of internal nodes in the species tree |
| `n_constraints` | Number of valid lines in the constraint file (constraint count) |
| `cpu_time` | CPU time (seconds, `time.perf_counter`) |
| `peak_memory_mb` | Peak memory (MB, `tracemalloc`) |
| `best_value` | Objective value of the best ranking (the smaller of greedy / mixing) |
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
