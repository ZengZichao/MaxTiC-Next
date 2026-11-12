# RANGER-DTLx 适配器规格说明书

> 🌐 English version: [ranger_dtlx.en.md](ranger_dtlx.en.md)

- **工具**：RANGER-DTLx（DTL 调和 + AggregateRanger），来自 `RANGER-DTLx-代码`。
- **适配器**：`maxtic_next.constraints.adapters.ranger_dtl`（`RangerDTLAdapter` /
  `convert_from_ranger_dtl`），解析核心 `_recon_report`（与 ARTra 共用）。
- **工具名（注册表 / CLI）**：`ranger`（`--from ranger`）。

## 1. 真实输出格式（源码求证）

RANGER-DTLx 逐基因节点打印人类可读调和报告。块头由
`RANGER-DTLx-代码/Ranger-DTLx/DTL-algorithm.h:2170` 打印；转移事件行则由
**同一文件的三处分支之一**生成，父节点片段随映射节点类型而变——`:1840-1841`（映射到
`isEdge` 边节点）打 `Transfer, Mapping --> …, Edge, Parent = …, Recipient --> …`、
`:1864-1865`（映射到 `isDummy` dummy 节点）打
`Transfer, Mapping --> …, Parent = …, Recipient --> …`、`:1886`
（`isDummy` 判定的 else 支，普通物种节点即由此打印）两者皆无。下文示例统一按
`:1840-1841` 那一支（带 `Edge,`）书写：

```
 ------------ Reconciliation for Gene Tree 1 (rooted) -------------
g1_0: Leaf Node
m1 = LCA[g1, g2]: Speciation, Mapping --> 59, Edge, Parent = 69
m2 = LCA[g3, g4]: Duplication, Mapping --> 61, Edge, Parent = 65
m3 = LCA[g1, g4]: Transfer, Mapping --> 61, Edge, Parent = 65, Recipient --> 62
```

> 分块只匹配 `Reconciliation for Gene Tree` 这一**前缀**（`_recon_report._BLOCK_RE`），故
> 本仓库样本（`tests/data/ranger_dtl_example/FAM1.dtl`）里的简写块头
> `Reconciliation for Gene Tree 1:` 与上游实际的 `(rooted)` 块头**都能识别**。

一个文件可含多个 `Reconciliation for Gene Tree N` 块。
⚠️ **块数不是采样分母**：上游按**输入基因树**逐棵打印块头，块 = 另一棵树，
而非另一次采样。支持度分母现优先取报告**自己声明**的
`Total number of optimal solutions: N`（`weight_semantics = support_over_declared_samples`）；
报告里出现**多个**声明值时（RANGER-DTLx 为每棵基因树各打印一行）取**最大值**作分母并告警；
仅当块数 > 1 且各块叶子集合一致（确为同一家族的重复调和）才退回
`support_over_replicate_blocks`（ecceTERA 不校验各块叶子集合，见 [eccetera.md](eccetera.md) §2）；
两者都不成立时权重按**整数计数**输出（`integer_count`），并明确告知
`--ale-min-support` 阈值此时**不具 [0,1] 支持度语义**。

> **格式已求证**：本适配器解析的就是 RANGER-DTLx 实际写出的 `Mapping --> / Recipient -->`
> 记录，字段来源已对照其源码与样本文件逐项核实。

## 2. 字段提取与语义映射

- 转移行的判定：`parse_transfer_line` 要求同一行**同时**含 `Mapping` 与 `Recipient`
  （只有转移事件才有 `Recipient`，损失 / 复制 / 物种形成行没有）；
- `Mapping --> D`（取到逗号 / 行尾）→ **供体** D；`Recipient --> R`（首个 token）→ **受体** R；
- 产出约束 `(D, R)`，语义与 ALE `T@D->R` 一致；
- 叶子行 `X: Leaf Node` 用于家族规模估算；
- Speciation / Duplication / Leaf 行不含 `Recipient -->`，自动忽略。

## 3. 支持度 / 权重 / 距离

块内计数 → 互反减法 → `support = count / denominator`（分母 = 采样数，即工具声明值；分母裁定
见 §1）→ `> min_support` 保留；
`distance = distance_from(receptor, donor)`；**自环 `donor == receptor` 例外**：写
`distance = None`（无距离信息），不被 `--min-transfer-distance` 删除。
供体端点取 `Mapping -->` 的**供体本身**（`donor_endpoint_convention = donor_itself`），
与 ALE 的 `parent_of_donor` **相差一个物种层级**，不可与 ALE 约束直接混加。

## 4. 错误处理与验证

| 情形 | 处理 |
|------|------|
| donor / receptor 不在物种树标签集合 | 该转移丢弃（不产伪约束） |
| 端点丢弃的可观测性 | 每条被丢弃的转移**按原因计数**（`dropped_by_label_miss`；数字 ID 反解仅在 ecceTERA 路径执行，本适配器 `numeric_id_resolutions` 恒为 0），转换结束打印 stderr 诊断摘要并挂在 `ConstraintSet.diagnostics`；端点命中率为 **0** 抛 `LabelMismatchError`，低于 `min_endpoint_hit_rate`（默认 `0.5`）告警（`--quiet-adapters` 只关打印，不影响计数与抛错） |
| 家族规模（叶子行数）`<= min_family_size` | 跳过 |
| 家族规模**探测失败**（报告里没有 `X: Leaf Node` 行） | **跳过** `min_family_size` 过滤并告警（`resolve_family_size_filter`），不再清空全部转移 |
| 无块头但有转移 | 块数记为 1 |
| 并行单文件异常 / 缓存命中 | 顺序回退 / 复用 |

## 5. 已知限制

- 端点为物种标签（内部节点常为数字，如 `61`），须与物种树 bootstrap 标签一致；
- AggregateRanger 的多样本聚合建议合并到同一文件的多个块，或分文件后由适配器跨文件求和。
- 本节为摘要：端点层级、权重口径与 diagnostics 计数等**完整语义**以规范性语义注记
  [README_semantics.md](../../src/maxtic_next/constraints/adapters/README_semantics.md) /
  [README_semantics.en.md](../../src/maxtic_next/constraints/adapters/README_semantics.en.md)
  为准。

## 6. 用法

```bash
maxtic-next species.tree FAM1.dtl --from ranger --ale-min-family-size 5
```

```python
from maxtic_next.constraints.adapters import convert_from_ranger_dtl
cset = convert_from_ranger_dtl(species_tree, ["FAM1.dtl", "FAM2.dtl"], min_family_size=5)
```

## 7. 样本

`tests/data/ranger_dtl_example/FAM1.dtl` + `species.tree`。回归测试见
`tests/test_ranger_dtl_adapter.py`、`tests/test_recon_report.py`。
