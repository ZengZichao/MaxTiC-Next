# 适配器语义说明（端点层级 / 权重口径 / 可观测性）

本文件是 `src/maxtic_next/constraints/adapters/` 内的**规范性语义注记**。
它回答三个问题：每个工具的 donor 到底取哪一层、权重到底是不是 [0,1] 支持度、
以及“约束变少”时用户如何看见。CLI/手册（`docs/adapters/*.md`）不重复这些内容时，
以本文件为准。

## 1. 供体端点相差一个物种层级

一条转移事件 `donor -> receptor` 可以写成两种 MaxTiC 约束，二者**不等价**：

| 路径 | donor 端点取什么 | 诊断值 `donor_endpoint_convention` | 上游证据 |
|---|---|---|---|
| ALE（`.uml_rec`） | **供体的父**：`parent(donnor)`（trf）或 `donnor_search()` 上溯到的谱系（rec） | `parent_of_donor` | 原版 `constraints_from_reconciliations.py:160-161`；`ale.py` 忠实复刻 |
| RANGER-DTLx / ARTra | **供体本身**：`Mapping --> <donor>` | `donor_itself` | `DTL-algorithm.h:1840` |
| ecceTERA | **供体本身**：父 clade 的 `<branchingOut speciesLocation="D">` | `donor_itself` | `DTLGraph.cpp:2738-2869` |
| AleRax | **供体本身**：`*_transfers.txt` 第一列 | `donor_itself` | `scripts/extract_families_transfer.py` |

原版 MaxTiC 的 ALE 集成同时把“受体的子”作为下界表述（ALE `maxtic/README.md:131`：
"the **father of the donor** branch should be older than the **child of the
receptor** branch"），而 `Mapping -->`/`branchingOut` 给的是事件所在谱系本身。
**本包不换算层级**：每个适配器保持其工具的原生约定（改了就不再等价于该工具的
官方产物），只把约定**写进诊断**。

`registry.merge_constraint_sets(sets)` 与 `registry.convert_auto(...,
allow_mixed=True)` 在合并两种约定时：

* stderr 打印“待合并的约束集来自**不同的供体端点层级约定**”警告；
* `ConstraintSet.diagnostics["convention_conflict"] = ["donor_itself",
  "parent_of_donor"]`（同样地，口径不可通约时记 `semantics_conflict`）。

因此 `--from-auto` 混用 ALE 与 RANGER 系输出不会再静默地把两代端点相加进同一条
edge——它会说出来。

## 2. 权重口径：块数不是采样分母

`diagnostics["weight_semantics"]` 取以下四值之一，随每条约束的
`metadata["weight_semantics"]` 一起返回：

| 口径 | 含义 | 何时使用 |
|---|---|---|
| `support_over_declared_samples` | `count / N`，N 为**工具自己声明**的样本数 | ALE 的 `<N> reconciled`；RANGER-DTLx/ARTra 的 `Total number of optimal solutions: N`（`DTL-algorithm.h:1785/:2224`、`output.txt:206`） |
| `support_over_replicate_blocks` | `count / 块数` | 无声明值，且各块叶子集合完全一致（确为同一家族的重复调和）；会附警告 |
| `integer_count` | **整数计数**，不是 [0,1] 支持度 | 无声明值且块数不可信（各块叶子集合不同 = 多棵不同基因树），或一个样本数都拿不到 |
| `posterior_frequency` | 上游已预聚合的后验频率 | AleRax `*_meanTransfers.txt`（不再除以任何分母） |

要点：

1. **真实 ARTra 输出**只打印 1 个最优解却声明 `Total number of optimal
   solutions: 24` ⇒ 分母是 24。旧实现按“块数 = 1”归一化，于是全部权重为 `1.0`，
   `--min-support 0.95` 一类阈值成了空操作。
2. 口径为 `integer_count` 时，`min_support` 仍在**同一数值**上比较（不会假装通过），
   但 `[0,1]` 语义不再成立：`diagnostics["min_support_is_fractional"] = False`，
   并在 stderr 明确警告。
3. 不同口径的权重**不可相加比较**：ALE 的 0.3 与 ARTra 的计数 3 不是同一回事。

## 3. 丢弃必须可见

每个适配器的 `convert()` 返回的 `ConstraintSet.diagnostics` 至少含：

```
tool, donor_endpoint_convention, weight_semantics, sample_denominator,
files_seen, blocks_seen, trees_seen, transfers_seen, transfers_resolved,
transfers_unpaired, numeric_id_resolutions, constraints_kept, edges_kept,
dropped_by_label_miss, dropped_by_support, dropped_by_family_size,
self_loop_transfers, family_size_probe_ok, endpoint_hit_rate,
unresolved_endpoint_samples, species_label_samples, warnings[], per_file[]
```

并打印一行 stderr 摘要：

```
[artra] donor 端点约定=donor_itself 权重口径=support_over_declared_samples 分母=24
端点命中率=100.0% | 文件=1 块=1 转移=13 保留=13 丢弃[标签不匹配=0, 支持度=0,
家族规模=0, 家族规模探测失败=0, 自环=0]
```

三条硬性规则：

* **家族规模探测失败 ⇒ 跳过过滤**：`Leaf Node` 行（RANGER/ARTra）
  或 `<leaf>` 元素（ecceTERA）探测不到时，`family_size` 会是 0，旧实现
  `if family_size <= min_family_size: transfers = {}` 于是**清空全部转移**
  （`--ale-min-family-size 0` 必然触发）。现在跳过过滤 + 告警。
* **端点命中率守卫**：ecceTERA 等适配器对内部物种节点写的是数字
  `node->getId()`（`DTLGraph.cpp:2818-2821`；编号规则见
  `MySpeciesTree.cpp:112-136` 的自底向上逐层广度优先）。`_species.numeric_species_id_map`
  据此反解；反解后仍**全部**不命中 ⇒ 抛 `LabelMismatchError`（错误信息同时给出
  未命中标本与物种树标签样本），部分命中低于阈值（默认 50%，可经
  `min_endpoint_hit_rate` 配置）⇒ 告警。绝不静默返回 `{}`。
* **XML 解析失败 ⇒ 抛错**：`UpstreamParseError` 携带文件、行列号与原始异常；
  兼容旧签名的纯函数 `parse_recphyloxml(..., strict=False)` 至少打印含文件与行列的
  响亮警告，并另支持 `strict=True`。损坏文件与“该工具没检出转移”从此可区分。

## 4. 距离与自环

* `_species.distance_from` 对**不在物种树中的标签返回 `None`**（旧实现返回 `0.0`，
  与“供体=受体的真距离 0”不可区分，且该假 0 会被 `--min-transfer-distance`
  当成合法距离把约束静默删除）。`None` 在
  `ConstraintSet.filter_by_distance` 中与“文本输入没有第 4 列”同义：**保留**。
* `filter_by_distance(d)` 现在返回并记录 `diagnostics["distance_filter"] =
  {threshold, before, kept, dropped, without_distance_column,
  self_loops_dropped, ignored}`，其中 `ignored=True` 表示“阈值 > 0 但一条都没删”。
  调用方（`ranking/ranker.py`、`dry_run.py`）必须据此报 **error 级**“该选项对本次
  输入被忽略”，而不是反过来宣称“会丢弃全部约束”。
* 自环（donor == receptor，“供体=受体谱系”）计入 `self_loop_transfers`，
  不再从统计里消失（；文本路径由原版 `ranker` 的 "to itself" 分支处理）。

## 5. 缓存键含拓扑指纹

`_species.species_cache_identity(tree)` = `v{CACHE_KEY_VERSION}|` +
标签集合摘要 + **拓扑指纹**（`sha256` over 排序后的 `parent_label->child_label`
对）。ALE / RANGER / ARTra / ecceTERA 的文件级 pickle 缓存键都使用它，且带版本号
前缀：任何只含标签摘要的条目**永不会被读到**。

因此 `--ale-cache-dir` 在“标签集合相同、父子关系不同”的两棵树之间复用目录时，
不会再命中旧拓扑的 donor（无缓存真值 `E`、命中缓存得 `F`）。

## 6. ALE 默认口径 = `trf`

`ALEAdapter` / `convert_from_ale` 的 `source` 默认为 **`"trf"`**：ALE 官方仓库自己的
MaxTiC 集成（ALE 源码包内 `maxtic/README.md:131-169`）规定输入
`constraints_from_transfers` 由**转移事件**产生（donor = 供体的父），实测该路径与
官方产物逐条一致（5,328 边键 / 27,576 行 / 权重相等），而 `"rec"` 只有 13,995 条、
总权重 5482.98 vs 6438.6 —— 默认值换了整个约束集。`"rec"` 仍完整保留、可选。
`convert()` 结束时在 stderr 打印一行：本次由哪种口径产出、边数与总权重。

> 上层接线（`api.py` 的 `ale_source` 形参默认值、`cli.py` 的 `--ale-source`
> `default=`）需同步为 `"trf"`，否则 CLI/API 仍会覆盖本包的默认值。
