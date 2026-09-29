# 04 · Python API

MaxTiC-Next 提供简洁的 Python API，便于嵌入脚本与流水线。核心导出在包顶层：

```python
from maxtic_next import rank, build_constraints
```

## 4.1 `rank(...)` —— 端到端排序

下面的签名逐字取自 `src/maxtic_next/api.py`（默认常量解析自 `config.py`）：

```python
def rank(species_tree_path, constraints_path,
         seed=42, local_search=0.0, temperature=0.001,
         random_type=0, min_transfer_distance=0,
         threshold_constraints=0.0, random_trees=0,
         output_prefix=None, print_summary=True,
         *, from_ale=False, from_tool=None,
         target_clade=None, ancestor_map=False, dry_run=False,
         html_report=True, ale_min_support=0.05,
         ale_min_family_size=5, ale_cache_dir=None,
         ale_source="trf", ale_parallel="process",
         artra_transfer_kind="all",
         adapter_min_endpoint_hit_rate=None, adapter_quiet=False,
         mcmc=False, mcmc_iters=1000, mcmc_temperature=0.0,
         mcmc_burn_in=None, mcmc_thin=None,
         incremental=True, checkpoint_path=None, checkpoint_interval=60.0,
         constraints_out=None, output_style="short",
         force=False, local_search_max_iterations=0,
         top_k=50, stop_check=None) -> Result
```

- `constraints_path` 可为 **str 或 list**；多个文件会全部解析并合并（边键聚合与文件顺序无关），
  产物前缀取**第一个**文件（多文件时 stdout 不变、**stderr** 打一条 NOTE，并写入
  `result.run_metadata["multi_file_prefix_note"]`）。
- 第 10 个位置参数之后（`from_ale` 起）全部是**关键字专用**参数。
- `mcmc_temperature=0.0` 是 **`auto` 哨兵**（`max(总权重,1)/100`），不是“零温度”；传 `0.0`
  与 CLI 传 `--mcmc-temperature auto` 等价。传负数或非正数会报错。
- `top_k=50` 是**近优解收集容量**（CLI `--near-optimal-top-k`）：稳健性/敏感性摘要最多考察多少个
  **去重**近优排序，须为 ≥ 1 的整数（`ranking/ranker.py:check_near_optimal_top_k` 校验），
  仅 `local_search > 0` 时有意义，生效值见 `run_metadata["near_optimal_top_k"]`（03 章 3.6b）。
- `stop_check=None` 是**协作式取消谓词**（`()`→`bool`）：在局部搜索与 MCMC 链的**迭代边界**
  打断计算并交付截至当下的结果，Studio 的「取消」按钮就用它（08 章 8.9）。
- `species_tree_path` / `constraints_path` 同样**透明接受** `.gz` / gzip 流与 tar 归档
  （05 章 5.7）：`api.rank` 会先把 `.tgz`/`.tar.gz`/`.tar` 安全解压到进程临时目录、把成员
  当作输入（`run_metadata["archives_expanded"]` 留痕），单文件 gzip 则由读取层就地解压。
- `threshold_constraints` 取值域 **[0, 1]**；`random_type` ∈ {0,1,2}；`temperature` > 0；
  `local_search` ≥ 0；`random_trees` ≥ 0；`mcmc_iters` > 0；`top_k` ≥ 1。这些都在
  `Ranker._validate_params` 里强制，违反即 `ValueError`（**包括 `dry_run=True` 路径**）。
- `adapter_min_endpoint_hit_rate=None` 表示沿用适配器层默认（`0.5`）。
- `force=False` 时，若产物已存在会抛 `FileExistsError`（CLI 转为退出码 3）。

### 参数与 CLI 的对应

| Python 参数 | CLI 参数 |
|-------------|----------|
| `species_tree_path` | 位置参数 `SPECIES_TREE` |
| `constraints_path`（str 或 list） | 位置参数 `CONSTRAINTS` |
| `seed` | `--seed` |
| `local_search` | `--local-search/--ls` |
| `temperature` | `--temperature/--t` |
| `random_type` | `--random-type/--r` |
| `min_transfer_distance` | `--min-transfer-distance/--d` |
| `threshold_constraints` | `--threshold-constraints/--ts` |
| `random_trees` | `--random-trees/--rd` |
| `from_tool` | `--from` |
| `from_ale` | `--from-ale` |
| `ale_min_support` | `--ale-min-support` |
| `ale_min_family_size` | `--ale-min-family-size` |
| `ale_cache_dir` | `--ale-cache-dir` |
| `ale_source`（默认 `"trf"`） | `--ale-source` |
| `ale_parallel` | `--ale-parallel` |
| `artra_transfer_kind` | `--artra-transfer-kind` |
| `adapter_min_endpoint_hit_rate` | `--min-endpoint-hit-rate` |
| `adapter_quiet` | `--quiet-adapters` |
| `target_clade` / `ancestor_map` | `--target-clade` / `--target-clade-ancestor-map` |
| `dry_run` | `--dry-run` |
| `html_report`（默认 True） | `--no-html`（取反） |
| `mcmc` / `mcmc_iters` / `mcmc_temperature` | `--mcmc` / `--mcmc-iters` / `--mcmc-temperature` |
| `mcmc_burn_in` / `mcmc_thin` | `--mcmc-burn-in` / `--mcmc-thin` |
| `incremental` / `checkpoint_path` / `checkpoint_interval` | `--incremental` / `--checkpoint` / `--checkpoint-interval` |
| `local_search_max_iterations` | `--local-search-max-iters` |
| `top_k` | `--near-optimal-top-k` |
| `constraints_out` | `-o/--constraints-out` |
| `output_style` | `--output-style` |
| `force` | `-f/--force` |
| `output_prefix` | （API 专属，等价 `-p`）输出文件前缀，默认取第一个约束文件路径 |
| `print_summary` | （API 专属）是否打印 stdout 摘要 |
| `stop_check` | （API/GUI 专属，**无对应 CLI 选项**）协作式取消谓词，见 4.4.6 |

> CLI 独有：`--from-auto`（等价 `from_tool="auto"`）、`--version`、`--help`。
> 所有 CLI 选项都抵达 `api.rank`，无一“死选项”；反向不成立——`stop_check` 只在
> Python / Studio 一侧存在，CLI 靠 `Ctrl-C` 之外没有中断开关。

## 4.2 `Result` 返回对象

`rank()` 返回 `maxtic_next.ranking.ranker.Result`（dataclass），关键字段：

| 字段 | 类型 | 含义 |
|------|------|------|
| `constraint_file` | str | 本次前缀来源（第一个约束文件） |
| `input_order` | list[str] | 输入树给出的内部节点排序 |
| `greedy_order` | list[str] | 贪婪启发式排序 |
| `mixing_order` | list[str] | 混合启发式排序 |
| `best_order` | list[str] | 最终最优排序 |
| `best_source` | str | `"greedy heuristic"` / `"mixing heuristic"`；`--ls > 0` 且有改进时为 `"… + local search"` |
| `values` | dict | `{"input":..,"greedy":..,"mixing":..}`；开局部搜索时另有 `"local_search"` 与 `"best"` |
| `uninformative` | dict | `{"total","to_desc","to_leaf","to_anc","to_itself"}` |
| `uninformative_percent` | float/None | **= total ÷ `total_weight` × 100**（，非原版双重计入的分母） |
| `from_leaf` | float | 供体来自叶子的约束权重（无法满足，忽略于打分） |
| `trivial_conflict` | float | 平凡冲突权重 |
| `total_weight` | float | 约束总权重（`--ts` 生效时为删边**之后**的口径） |
| `removed_by_threshold_weight` / `removed_by_threshold_count` | float / int | `--ts` 删除的权重与条数 |
| `informative_count` | int | 信息性约束条数 |
| `internal_node_count` | int | 内部节点数 |
| `ranked_newick` | str | 排序树 Newick |
| `similarity_to_input` | float | 与输入树的 Kendall 相似度 |
| `conflict_with_input` / `partial_total` | float | 摘要里“与最优树一致却与输入树冲突”比例的分子/分母 |
| `informative_lines` / `conflicting_lines` / `partial_lines` | list[str] | 三个文件的原始行 |
| `informative_file` / `conflicting_file` / `partial_order_file` | str | 三输出文件路径 |
| `html_report_file` | str | HTML 报告路径（生成时） |
| `sensitivity_summary` | dict/None | 稳健性/敏感性摘要（局部搜索开启时） |
| `mcmc_samples` | list/None | MCMC 线性扩展样本（`--mcmc` 时） |
| `random_stats` | dict/None | 随机树 p 值统计（`--random-trees>0` 时），含 `value_pvalue`、`sim_pvalue`、`correction="(k+1)/(n+1)"`、`tested_order="delivered best order (after local search)"`、`distribution_file`（第 4 个产物路径） |
| `warnings` | list[str] | 运行期口径声明/告警（与 stdout 末尾那些行同源） |
| `run_metadata` | dict | 运行元数据：`params`、`dry_run`/`dry_run_ok`/`dry_run_n_errors`、`adapter_diagnostics`、`constraint_files_merged`、`multi_file_prefix_note`、`mcmc_status`、`near_optimal_top_k`（`--near-optimal-top-k` 的生效值）、`cancelled` / `stop_check_requested` / `mcmc_stopped_early`（协作式取消）、`archives_expanded` / `archive_notes`（归档输入展开了哪些成员）等 |
| `dry_run_report` | str | 预检报告（`dry_run=True` 时） |

> **没有“逐节点 MTC 分数”这一量**：MaxTiC 的目标函数是整个排序的被违反权重之和
> （`ranking/value.py`），任何按节点展示的分数都只是输入约束的描述量，不解释排名。

## 4.3 `build_constraints(path)` —— 仅解析

```python
from maxtic_next import build_constraints
cset = build_constraints("constraints.tsv")   # 返回 ConstraintSet
for c in cset.constraints:
    print(c.donor, c.receptor, c.weight, c.metadata)
```

按距离过滤**不**在此处发生（由排序阶段的 `ConstraintSet.filter_by_distance` 执行）。
非法权重（`nan`/`inf`/负值）在解析层即报 `ValueError`，信息含文件名与行号。

## 4.4 典型用法

### 4.4.1 基本排序

```python
from maxtic_next import rank

r = rank("examples/minitree.tree", "examples/Cyano_CUTConstraints.tsv",
         seed=42, print_summary=False, html_report=False)
print("最优启发式:", r.best_source)          # greedy heuristic
print("与输入相似度:", r.similarity_to_input)  # 1.0
print("三文件:", r.informative_file, r.conflicting_file, r.partial_order_file)
```

### 4.4.2 从上游工具输出直接排序

```python
# RANGER-DTLx 调和报告
r = rank("species.tree", ["FAM1.dtl", "FAM2.dtl"], from_tool="ranger",
         ale_min_family_size=5, print_summary=False, html_report=False)

# 自动检测
r = rank("species.tree", ["out.recphyloxml"], from_tool="auto", html_report=False)

# 取适配器口径与丢弃计数
print(r.run_metadata["adapter_diagnostics"])
```

ALE 默认 `ale_source="trf"`（转移事件、ALE 官方 MaxTiC 集成口径）。要与 `rec` 口径对齐，
显式传 `ale_source="rec"` —— 两者是**不同的约束集**，不是等价写法。

### 4.4.3 两阶段：先产统一约束，后排序

```python
# 阶段 1：把 ALE .uml_rec 转为统一约束文件（不排序）
rank("species.tree", ["fam1.uml_rec", "fam2.uml_rec"], from_tool="ale",
   ale_cache_dir="cache/", constraints_out="constraints.tsv")

# 阶段 2：对约束文件排序
rank("species.tree", "constraints.tsv", local_search=180)
```

### 4.4.4 局部搜索 + 稳健性摘要

```python
r = rank("species.tree", "constraints.tsv", local_search=30,
         print_summary=False)
summ = r.sensitivity_summary   # dict，含 node_position_distribution 等
print(summ["title"])           # "基于局部搜索访问解的稳健性/敏感性摘要"
print("收集近优序数:", summ["n_orders_collected"])
print("上限 K:", summ["top_k"], r.run_metadata["near_optimal_top_k"])
print("链实际见过的不同序数:", summ["n_unique_orders_total"])
```

要跨机器复现这条路径，加 `local_search_max_iterations=2000000`。
摘要只统计目标值最小的 `top_k` 个**去重**近优序（默认 50），因此报告这些频率时
**必须**一并给出 K；调大它（`top_k=400`）只是拓宽支持集，**不会**让它变成后验
（见 03 章 3.6b / 08 章 8.1b）。

### 4.4.5 只要预检结论，不要排序

```python
r = rank("species.tree", "constraints.tsv", dry_run=True)
ok = r.run_metadata["dry_run_ok"]        # bool —— CLI 退出码的唯一依据
print(r.dry_run_report)                  # 人类可读报告
print(r.run_metadata["dry_run_n_errors"])
```

适配器模式下同样成立：上游输出**解析失败**时 `dry_run_ok` 为 `False`、
`dry_run_n_errors ≥ 1`（CLI 因此退出 1），而不是像旧实现那样判“通过”。

```python
r = rank("examples/adapters/ale/species.tree", "truncated.recphyloxml",
         from_tool="eccetera", dry_run=True, print_summary=False)
r.run_metadata["dry_run_ok"]          # False
r.run_metadata["dry_run_n_errors"]    # 1
```

### 4.4.6 协作式取消（`stop_check`）

```python
import itertools
from maxtic_next import rank

calls = itertools.count()
r = rank("species.tree", "constraints.tsv", local_search=30.0, temperature=60.0,
         stop_check=lambda: next(calls) > 50)      # 第 50 次检查之后转真

r.run_metadata["cancelled"]                     # True
r.run_metadata["local_search_iterations"]       # 51 —— 计划 30 秒，在迭代边界停下
r.values["best"]                                # 截至取消点的最优值，照常交付
```

保证：只在**迭代边界**停、返回历史最优序、照常写产物、不抛异常；不传 `stop_check` 时
行为与原来逐字节一致。`MCMCSampler.sample(stop_check=...)` 同样接受它；**唯一例外**是
一条还没产出任何样本就被取消的 MCMC 链——那会抛 `RuntimeError` 而不是返回空统计。
底层开关是 `Ranker.run(stop_check=...)` → `optimisation_locale(stop_check=...)`，
后者把 `stats_out["cancelled"]` 写回 `run_metadata`。详见 08 章 8.9。

### 4.4.7 压缩输入直接喂

```python
rank("species.tree.gz", "constraints.tsv.gz")                      # 读取层透明解压
rank("species.tree", "examples/reconciliations.tgz", from_tool="ale")
#  ↑ 归档被安全解压到进程临时目录，成员各算一份输入
r.run_metadata["archives_expanded"]    # {归档路径: [成员路径…]}，供追溯
```

多成员归档被当作**一个**文本输入读取时抛 `CompressedArchiveError`（`ValueError` 子类），
消息含 `tar -xzf` / `tar -xzOf` 命令。边界与实测见 05 章 5.7。

## 4.5 适配器底层 API（进阶）

```python
from maxtic_next.tree.tree import Tree
from maxtic_next.constraints.adapters import registry

species = Tree(); species.read_newick(open("species.tree").read())

# 显式工具
cset = registry.convert("ranger", species, ["FAM1.dtl"], min_family_size=0)
# 自动检测
cset = registry.convert_auto(species, ["out.recphyloxml"])
# 可用工具名
print(registry.REGISTRY.names())   # ['ale','alerax','artra','eccetera','ranger']
```

`convert` / `convert_auto` 返回的 `ConstraintSet` 带 `.diagnostics`（端点命中率、各类丢弃计数、
`donor_endpoint_convention`、`weight_semantics`）。`--from auto` 遇到**混合格式**时不再抛错
了事：按工具分组各自解析后 `registry.merge_constraint_sets(...)` 合并，同时把“供体端点层级
不同 / 权重口径不可通约”打到 stderr 并记入 `diagnostics`。

各工具便捷函数：

```python
from maxtic_next.constraints.adapters import (
    convert_from_ale, convert_from_ranger_dtl, convert_from_eccetera,
    convert_from_artra, convert_from_alerax,
)
```

> 上一章：[03 · 命令行参考](03_cli_reference.md) ｜ 下一章：[05 · 输入输出格式](05_io_formats.md)
