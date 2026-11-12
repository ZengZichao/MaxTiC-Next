# MaxTiC-Next 真实数据集基准目录

> 🌐 English version: [README.en.md](README.en.md)

本目录用于存放真实大规模 HGT（水平基因转移）数据集，供基准套件模块
`src/maxtic_next/benchmarks/suite.py` 的 `--real` 模式（入口
`python -m maxtic_next.benchmarks.suite --real`）进行性能压测。

## 目录结构

```
benchmarks/datasets/
├── README.md            ← 本文件
├── README.en.md         ← 英文版
├── <dataset_name>/
│   ├── <dataset_name>.tree          ← 物种树（Newick 格式）
│   └── <dataset_name>.tsv           ← 约束文件（TSV / 文本双格式）
└── ...
```

> 目录扫描**不递归子目录**，故 `--data-dir` 应指到直接包含树与约束的那一层
> （例如 `benchmarks/datasets/<dataset_name>/`）。扫描还会跳过 MaxTiC 自身产出的输出
> 文件（三文件 / distribution / HTML 报告等后缀，见 `suite._is_maxtic_output`），以免把
> 排序输出重新当作输入约束摄入。

## 如何添加真实数据集

### 1. 准备数据文件

每组数据集需要两个文件：

- **物种树文件**：Newick 格式，扩展名为 `.tree` / `.nwk` / `.newick` / `.tre`
- **约束文件**：MaxTiC 约束格式（三列：`donor\treceptor\tweight`），
  扩展名为 `.tsv` / `.txt` / `.csv` / `.constraints`

### 2. 命名约定

为便于自动配对，建议树文件与约束文件使用相同的文件名前缀（去掉扩展名后相同）：

```
my_dataset/
├── cyano.tree              ← 物种树
└── cyano.tsv               ← 约束
```

如果目录内只有一个树文件和一个约束文件，即使名称不同也会自动配对。

### 3. 运行基准测试

> **本仓库未随附真实数据集**（体积原因）：`benchmarks/datasets/` 目前只有说明文件，
> 需按上文约定自行放置数据。`--data-dir` 的**默认值是项目根目录下的 `examples/`**
> （只含 `minitree.tree` + `Cyano_CUTConstraints.tsv` 这一组小型演示数据）。

```bash
# 使用默认 examples/ 目录中的数据
python -m maxtic_next.benchmarks.suite --real

# 指定自定义数据集目录（占位符 <path/to/datasets> 请换成你实际放置数据的目录）
python -m maxtic_next.benchmarks.suite --real --data-dir <path/to/datasets>

# 指定种子和时间预算
python -m maxtic_next.benchmarks.suite --real --data-dir <path/to/datasets> --seed 42 --time-budget 120
```

若 `--data-dir` 指向的目录里没有 tree + constraints 配对（例如仍为空的
`benchmarks/datasets/`），`suite.py` 会打印「未找到 tree+constraints 数据对」的警告并
跳过全部数据组。

### 4. 输出指标

每组数据集的基准测试输出以下主要指标：

| 指标 | 说明 |
|------|------|
| `n_internal` | 物种树内部节点数 |
| `n_constraints` | 约束文件有效行数（约束数，不含 `#` 注释行与 FRQ 头行） |
| `n_distinct_pairs` | 去重后的不同 `(donor, receptor)` 对数 |
| `cons_*` | 约束文件的诚实规模明细：`cons_n_lines_nonempty` / `cons_n_constraints_parsed` / `cons_n_distinct_pairs` / `cons_n_comment_or_header_lines` |
| `cpu_time` | CPU 时间（秒，`time.process_time` 中位数；墙钟另报 `wall_time`） |
| `cpu_time_iqr` / `wall_time_iqr` | `cpu_time` / `wall_time` 的四分位距 |
| `wall_time` | 墙钟时间（秒，`time.perf_counter` 中位数） |
| `io_time` / `kernel_cpu_time` | 输出写入的 I/O 开销（CPU 中位数）/ 扣除 I/O 后的排序内核 CPU 时间 |
| `peak_memory_mb` | 峰值内存（MiB，`tracemalloc`；实现按 `peak / (1024 * 1024)` 计算，结果字典的键名仍为 `peak_memory_mb`） |
| `best_value` | 最优排序的目标值（greedy / mixing 中的较小者） |
| `delivered_value` | 交付排序目标值的独立重算结果 |
| `best_source` | 交付排序的来源（`greedy heuristic` / `mixing heuristic`，可能带 `+ local search`） |
| `correctness_ok` | 交付值独立重算是否通过正确性断言 |
| `similarity_to_input` | 最优排序与输入树的 Kendall 相似度 |
| `exceeded` | 是否超出时间预算 |

## 数据格式说明

### 物种树（Newick）

标准 Newick 格式，内部节点用 bootstrap 标签标注：

```
(((((CYAP8:1,CYAP0:1)42:7,(CYAA5:7,UCYNA:7)41:1)59:2,...
```

### 约束文件（TSV）

每行一条约束，三列以制表符分隔：`donor\treceptor\tweight`

```
61	62	15.04
61	67	15.48
61	46	28.73
```

- `donor`：转移供体节点（内部节点标签）
- `receptor`：转移受体节点（内部节点标签）
- `weight`：约束权重（浮点数，通常为转移支持度）

## 注意事项

- 大规模数据集（内部节点数 > 500）可能需要较长时间，请合理设置 `--time-budget`。
- 基准测试使用 `tracemalloc` 测量峰值内存，仅追踪 Python 分配的内存，
  不包含 C 扩展的直接内存分配。
- 所有输出文件（informative / conflicting / partial order）默认写入临时目录，
  不污染工作区。
