# RANGER-DTLx 适配器规格说明书

- **工具**：RANGER-DTLx（DTL 调和 + AggregateRanger），来自 `RANGER-DTLx-代码`。
- **适配器**：`maxtic_next.constraints.adapters.ranger_dtl`（`RangerDTLAdapter` /
  `convert_from_ranger_dtl`），解析核心 `_recon_report`（与 ARTra 共用）。
- **工具名（注册表 / CLI）**：`ranger`（`--from ranger`）。

## 1. 真实输出格式（源码求证）

RANGER-DTLx 逐基因节点打印人类可读调和报告。转移事件行由
`RANGER-DTLx-代码/Ranger-DTLx/DTL-algorithm.h:1840` 生成：

```
Reconciliation for Gene Tree 1:
g1_0: Leaf Node
 = LCA[g1, g2]: Speciation, Mapping --> 59, Parent = 69
 = LCA[g3, g4]: Duplication, Mapping --> 61, Parent = 65
 = LCA[g1, g4]: Transfer, Mapping --> 61, Edge, Parent = 65, Recipient --> 62
```

一个文件可含多个 `Reconciliation for Gene Tree N` 块。
⚠️ **块数不是采样分母**：上游按**输入基因树**逐棵打印块头，块 = 另一棵树，
而非另一次采样。支持度分母现优先取报告**自己声明**的
`Total number of optimal solutions: N`（`weight_semantics = support_over_declared_samples`）；
仅当各块叶子集合一致（确为同一家族的重复调和）才退回 `support_over_replicate_blocks`；
两者都不成立时权重按**整数计数**输出（`integer_count`），并明确告知
`--ale-min-support` 阈值此时**不具 [0,1] 支持度语义**。

> **格式已求证**：本适配器解析的就是 RANGER-DTLx 实际写出的 `Mapping --> / Recipient -->`
> 记录，字段来源已对照其源码与样本文件逐项核实。

## 2. 字段提取与语义映射

- 只有**转移事件**才带 `Recipient -->`，据此触发（正则 `_recon_report`）；
- `Mapping --> D`（取到逗号 / 行尾）→ **供体** D；`Recipient --> R`（首个 token）→ **受体** R；
- 产出约束 `(D, R)`，语义与 ALE `T@D->R` 一致；
- 叶子行 `X: Leaf Node` 用于家族规模估算；
- Speciation / Duplication / Leaf 行不含 `Recipient -->`，自动忽略。

## 3. 支持度 / 权重 / 距离

块内计数 → 互反减法 → `support = 计数 / 声明样本数`（见 §2 的分母裁定）→ `> min_support` 保留；
`distance = distance_from(receptor, donor)`；**自环 `donor == receptor` 例外**：写
`distance = None`（无距离信息），不被 `--min-transfer-distance` 删除。
供体端点取 `Mapping -->` 的**供体本身**（`donor_endpoint_convention = donor_itself`），
与 ALE 的 `parent_of_donor` **相差一个物种层级**，不可与 ALE 约束直接混加。

## 4. 错误处理与验证

| 情形 | 处理 |
|------|------|
| donor / receptor 不在物种树标签集合 | 该转移丢弃（不产伪约束） |
| 家族规模（叶子行数）`<= min_family_size` | 跳过 |
| 无块头但有转移 | 块数记为 1 |
| 并行单文件异常 / 缓存命中 | 顺序回退 / 复用 |

## 5. 已知限制

- 端点为物种标签（内部节点常为数字，如 `61`），须与物种树 bootstrap 标签一致；
- AggregateRanger 的多样本聚合建议合并到同一文件的多个块，或分文件后由适配器跨文件求和。

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
