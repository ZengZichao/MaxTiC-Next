# MaxTiC-Next 真实数据集基准目录

本目录用于存放真实大规模 HGT（水平基因转移）数据集，供 `benchmarks/suite.py` 的
`--real` 模式进行性能压测。

## 目录结构

```
benchmarks/datasets/
├── README.md            ← 本文件
├── <dataset_name>/
│   ├── <dataset_name>.tree          ← 物种树（Newick 格式）
│   └── <dataset_name>.tsv           ← 约束文件（TSV / 文本双格式）
└── ...
```

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

```bash
# 使用默认 examples/ 目录中的数据
python -m maxtic_next.benchmarks.suite --real

# 指定自定义数据集目录
python -m maxtic_next.benchmarks.suite --real --data-dir benchmarks/datasets/

# 指定种子和时间预算
python -m maxtic_next.benchmarks.suite --real --data-dir benchmarks/datasets/ --seed 42 --time-budget 120
```

### 4. 输出指标

每组数据集的基准测试输出以下指标：

| 指标 | 说明 |
|------|------|
| `n_internal` | 物种树内部节点数 |
| `n_constraints` | 约束文件有效行数（约束数） |
| `cpu_time` | CPU 时间（秒，`time.perf_counter`） |
| `peak_memory_mb` | 峰值内存（MB，`tracemalloc`） |
| `best_value` | 最优排序的目标值（greedy / mixing 中的较小者） |
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
