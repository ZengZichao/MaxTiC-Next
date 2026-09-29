# 02 · 快速开始

本章用仓库自带的示例数据（蓝细菌小数据集）在 5 分钟内跑通一次完整排序。

> 本章所有 stdout / 文件清单均为**实际运行捕获**（Python 3.14 本机，命令逐字照抄即可复现）。

## 2.1 输入数据

- 物种树：`examples/minitree.tree`（13 个内部节点，Newick，内部标签在 bootstrap 字段，**单行**）
- 约束：`examples/Cyano_CUTConstraints.tsv`（空格格式 `donor receptor weight`，**123 行**，
  即原版随包分发的那份蓝细菌约束集）

物种树片段（内部节点标签为数字 `42/59/61/...`）：

```
(((((CYAP8:1,CYAP0:1)42:7,(CYAA5:7,UCYNA:7)41:1)59:2,...)69;
```

约束片段（`donor receptor weight`，制表符分隔）：

```
61	62	15.04
61	67	15.48
61	46	28.73
```

## 2.2 最小运行

```bash
maxtic-next examples/minitree.tree examples/Cyano_CUTConstraints.tsv --seed 42
```

未安装（或不想装）时，从仓库根目录用源码树直接跑：

```bash
PYTHONPATH=src python3 -m maxtic_next examples/minitree.tree \
    examples/Cyano_CUTConstraints.tsv --seed 42
```

两者完全等价（`maxtic-next` 就是 `maxtic_next.cli:main` 的控制台入口）。

## 2.3 期望的 stdout 摘要

以下为上述命令在仓库根目录的**实际输出**（HTML 报告同时生成；加 `--no-html` 只去掉
`*.html`，stdout 不变）。注意 `tree with` 后面是**两个空格**（与原版 `print` 一致）：

```
examples/Cyano_CUTConstraints.tsv
tree with  13 internal nodes
2218.4 total weight of constraints from transfers
307.45000000000005 uninformative (  307.45000000000005 to a descendant, 0.0 to a leaf 0.0 to an ancestor 0.0 to itself 13.9%)
676.0400000000001 trivially conflicting constraints (30%) (descendant to ancestor or trivial cycle)
value of the order given by the input tree 755.8800000000005 (34.0732059141724%)
value of the greedy heuristic 755.8800000000005 (34.0732059141724%)
value of the mixing heuristic: 755.8800000000005 (34.0732059141724%)
best order is the greedy heuristic
(((((CYAP8:1,CYAP0:1)42:7,(CYAA5:7,UCYNA:7)41:1)59:2,(MICAN:9,(CYAP2:3,CYAP7:3)43:6)45:1)61:2,(SYNY3:11,SYNP2:11)56:1)65:1,(TRIEI:6,(NOSA0:5,(NOSP7:4,(NOSS1:2,ANAVT:2)40:2)46:1)62:1)67:7)69;
Similarity of the best order compared with the input order 1.0
0.0 constraints in agreement with the best tree in conflict with input tree
Similarity of the order compared with the input order 1.0
NOTE: 12 条约束与物种树谱系边同向（权重合计 63.55），已被置为不可违反的拓扑硬约束，并已计入上方 uninformative 统计（在 total_weight 分母中只计一次，不重复计权）。
```

**如何读这份摘要**：

| 行 | 含义 |
|----|------|
| `examples/Cyano_CUTConstraints.tsv` | 本次产物前缀 = **第一个**约束文件的路径 |
| `tree with  13 internal nodes` | 待排序的内部节点数（两个空格为原版排版） |
| `2218.4 total weight of constraints from transfers` | 所有转移约束的权重之和（分母） |
| `uninformative (...)` | 被判为无信息（到后代 / 到叶子 / 到祖先 / 自环）的权重及占比；**百分比 = uninformative ÷ total_weight**（，原版分母把该项双重计入而报 `12%`，本版本报 `13.9%`） |
| `trivially conflicting` | 平凡冲突（后代→祖先或平凡环）权重及占比 |
| `value of the ... heuristic` | 三种排序各自“被违反约束权重”（越小越好）及占 total_weight 的百分比 |
| `best order is the greedy heuristic` | 最终选用的启发式；开局部搜索且确有改进时变为 `... + local search` |
| `after local search X rejected` / `best found solution Y (%)` | **仅 `--ls > 0`** 时追加两行，报告搜索后的真实值 |
| 一行 Newick | **排序树**：分支长度按排序位置重写 |
| `Similarity ...` / `constraints in agreement ...` | 与输入树排序的 Kendall 相似度（∈[0,1]，1 为完全一致）等 |
| `values from N random orders ...` 等 4 行 | **仅 `--random-trees > 0`** 时追加（见 08 章 8.7） |
| `NOTE: N 条约束与物种树谱系边同向…` / `WARNING: …` | 运行期口径声明与告警，恒打在摘要**末尾** |

> 末行的 NOTE 是**本版新增**（原版 `MaxTiC.py:444` 的逐条命中打印在其执行路径上不可达）。
> 它说明有 12 条约束与物种树本身的父子边同向，已被置为不可违反的拓扑硬约束。

## 2.4 生成的输出文件

默认 **short** 命名（前缀为**第一个**约束文件的路径）：

| 文件 | 内容 |
|------|------|
| `examples/Cyano_CUTConstraints.tsv.mt.informative.tsv` | 过滤后的加权信息性约束 |
| `examples/Cyano_CUTConstraints.tsv.mt.conflicts.tsv` | 与最优序冲突的约束 |
| `examples/Cyano_CUTConstraints.tsv.mt.partial_order.tsv` | 偏序（含 black/green 标记） |
| `examples/Cyano_CUTConstraints.tsv.html` | 交互式 HTML 报告（未加 `--no-html` 时） |

第 4 个数据文件只有开启随机树检验时才出现：

| 文件 | 出现条件 |
|------|----------|
| `<前缀>.mt.random_dist.tsv` | `--random-trees > 0`（每行 `value similarity`，共 N 行） |

> - 用 `--output-style legacy` 可切回原版长名（复现验证用），两种风格**内容逐字节一致**。
> - **默认拒绝覆盖**：若产物已存在，再次运行会以退出码 `3` 中止而不覆盖。
>   确认要覆盖请加 `-f/--force`，或用 `-p` 换个前缀。

## 2.5 打开 HTML 报告

```bash
# macOS
open examples/Cyano_CUTConstraints.tsv.html
# Linux
xdg-open examples/Cyano_CUTConstraints.tsv.html
```

报告为**纯静态**页面（内联/CDN 引入 Plotly），无需后端，可离线打开。
未安装 `jinja2` / `plotly` 时不会报错，而是**降级**为基础模板（见 09 章 FAQ）。

## 2.6 加一点“味道”

```bash
# 开局部搜索 10 秒（触发稳健性/敏感性摘要 + 两行 local search 输出）
maxtic-next examples/minitree.tree examples/Cyano_CUTConstraints.tsv --seed 42 --ls 10

# 随机树置换检验（追加 4 行 p 值输出 + 第 4 个产物文件）
maxtic-next examples/minitree.tree examples/Cyano_CUTConstraints.tsv --seed 42 --rd 50

# 从 RANGER-DTLx 调和报告直接排序（--ale-min-family-size 0 是样例夹具所需）
maxtic-next examples/minitree.tree examples/adapters/ranger/FAM1.dtl \
    --from ranger --ale-min-family-size 0

# 自动识别上游工具输出格式
maxtic-next examples/minitree.tree examples/adapters/eccetera/FAM1.recphyloxml \
    --from-auto --ale-min-family-size 0

# 只检查输入，不排序（合法输入退出码 0）
maxtic-next examples/minitree.tree examples/Cyano_CUTConstraints.tsv --dry-run
```

适配器模式会先向 **stderr** 打印一行口径声明与一行诊断摘要（；
不进 stdout，以免污染与原版对齐的摘要）：

```
[ranger] 警告：报告未声明样本数（无 ``Total number of optimal solutions: N`` 行）且只有 1 个调和块：权重按**整数计数**输出，min_support 阈值不具 [0,1] 支持度语义。
[ranger] donor 端点约定=donor_itself 权重口径=integer_count 分母=不可确立 端点命中率=100.0% | 文件=1 块=1 转移=2 保留=2 丢弃[标签不匹配=0, 支持度=0, 家族规模=0, 家族规模探测失败=0, 自环=0]
```

诊断量同时写入 `Result.run_metadata["adapter_diagnostics"]`；只想静默可加
`--quiet-adapters`（**只影响打印**，错误照旧抛出）。

## 2.7 Python 一行式

```python
from maxtic_next import rank
result = rank("examples/minitree.tree", "examples/Cyano_CUTConstraints.tsv", seed=42)
print(result.best_source, result.similarity_to_input)
# greedy heuristic 1.0
```

> 上一章：[01 · 安装](01_installation.md) ｜ 下一章：[03 · 命令行参考](03_cli_reference.md)
