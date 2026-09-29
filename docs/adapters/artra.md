# ARTra 适配器规格说明书

- **工具**：ARTra（Additive and Replacing Transfer classifier），来自 `ARTra-代码` /
  `ARTra-程序`。在 RANGER-DTL 风格调和上，用规则启发式 + 机器学习把每个转移分类为
  **加性（Additive）**或**替换（Replacing）**转移。
- **适配器**：`maxtic_next.constraints.adapters.artra`（`ARTraAdapter` /
  `convert_from_artra`），解析核心 `_recon_report`（与 RANGER-DTLx 共用）。
- **工具名（注册表 / CLI）**：`artra`（`--from artra`）。

## 1. 真实输出格式（样本求证）

ARTra 最终输出（`ARTra-程序/output.txt` 实证）与 RANGER-DTLx 调和报告**同形**，仅在
事件名多出 `Replacing` / `Additive` 前缀：

```
 ------------ Reconciliation for Gene Tree 1 (rooted) -------------
Species Tree:
(...);
Gene Tree:
(...);
Reconciliation:
H117_0: Leaf Node
 = LCA[H117, H14]: Replacing Transfer, Mapping --> H15, Recipient --> H125
 = LCA[H159, H17]: Additive Transfer, Mapping --> H17, Recipient --> H159
```

## 2. 字段提取与语义映射

- 与 RANGER 一致：`Mapping --> D` = 供体、`Recipient --> R` = 受体，产出约束 `(D, R)`；
- `Replacing Transfer` / `Additive Transfer` / 纯 `Transfer` 均被识别（触发条件为含
  `Recipient -->`）；
- `Species Tree:` / `Gene Tree:` 等段落头与 newick 行不含 `Recipient -->`，自动忽略。

## 3. 转移类别过滤（ARTra 专属）

MaxTiC 的时间约束语义与转移是加性还是替换无关，故**默认统计全部转移**
（`transfer_kind="all"`）。如研究只关心某一类别：

- `transfer_kind="replacing"`：仅统计替换转移；
- `transfer_kind="additive"`：仅统计加性转移。

（类别过滤在主进程内按行筛选后再解析，不使用并行 / 缓存，以免污染共享缓存键。）

## 4. 支持度 / 权重 / 距离 / 错误处理

同 RANGER-DTLx：块内计数 → 互反减法 → `support = 计数 / 声明样本数` → `> min_support`；
`distance = distance_from(receptor, donor)`（**自环 `donor == receptor` 例外**：写 `distance = None`，
不被 `--min-transfer-distance` 删除）；非法端点 / 过小家族被丢弃并**按原因计数上报**。
⚠️ ARTra 会打印 `Total number of optimal solutions: 24` 却只输出**一个**最优解，
因此以块数为分母在该工具上必然退化；分母不可确立时权重为**整数计数**，
`min_support` 不具 [0,1] 语义。供体端点同样取 `donor_itself`（与 ALE 相差一代）。

## 5. 已知限制

- 端点为物种标签，须与物种树标签一致；
- 机器学习分类置信度未纳入权重（权重仍为采样支持度）；如需以分类概率加权，可在上游
  预处理阶段生成加权约束后走文本约束路径。

## 6. 用法

```bash
maxtic-next species.tree FAM1.txt --from artra --artra-transfer-kind all
maxtic-next species.tree FAM1.txt --from artra --artra-transfer-kind replacing
```

```python
from maxtic_next.constraints.adapters import convert_from_artra
cset = convert_from_artra(species_tree, ["FAM1.txt"], transfer_kind="all", min_family_size=5)
```

## 7. 样本

`tests/data/artra_example/FAM1.txt` + `species.tree`。回归测试见
`tests/test_artra_adapter.py`。
