# ARTra 适配器规格说明书

> 🌐 English version: [artra.en.md](artra.en.md)

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
m3 = LCA[H117,H14]: Replacing Transfer, Mapping --> H15, Recipient --> H125
m19 = LCA[H159,H17]: Replacing Transfer, Mapping --> H17, Recipient --> H159
m93 = LCA[H10,H90]: Additive Transfer, Mapping --> H90, Recipient --> H10
```

> **块头亦与 RANGER-DTLx 相同**（不是 ARTra 特有，但由 ARTra **自己的 Java 代码**打印）：
> `参考软件-ARTra-代码/JavaCode/dtl_algo/DTLReconciliation.java:619`（stdout 版
> `printOutputOfDTRLConsole`）与 `:1387`（`PrintWriter` 版 `printOutputOfDTRL`）打印的就是
> ` ------------ Reconciliation for Gene Tree 1 (rooted) -------------`，与上游
> `RANGER-DTLx-代码/Ranger-DTLx/DTL-algorithm.h:2170` 的输出**同形**。分块只匹配
> `Reconciliation for Gene Tree` 这一**前缀**，故上游完整块头与 RANGER-DTLx 样本里的简写形式
> `Reconciliation for Gene Tree 1:` 都能被识别（ARTra 自身的样本
> `examples/adapters/artra/FAM1.txt` 用的就是完整 `(rooted)` 块头）。事件名由 `-h` 选定的分类
> 步骤改写为 `Replacing Transfer` / `Additive Transfer`：`DTLReconciliation.java:501` 处
> `heuristic` 默认为 `"ml"`，故 `ARTra-程序/output.txt` 的标签出自机器学习分类器
> `modifySolutionByML`（`:764` / `:767`）；规则启发式的改写点在 `:608` / `:610`（lost-gene）、
> `:1465` / `:1467`（gene-frequency）、`:1509-1518`（mapping-count）。事件行由 stdout 版
> `DTLReconciliation.java:637-641` 与 `PrintWriter` 版 `:1405-1409`（写 `output.txt`
> 的即此路径）打印——ARTra **从不写** `Edge, Parent =`
> （该片段只出现在 RANGER-DTLx 的 C++ 打印器里）。事件行的 `LCA[.., ..]` 内 ARTra 不写空格
> （`Node.java:466`）、RANGER-DTLx 写空格，解析只取 `Mapping -->` / `Recipient -->` 之后的标签，
> 两种写法都兼容。

## 2. 字段提取与语义映射

- 与 RANGER-DTLx 一致：`Mapping --> D` = 供体、`Recipient --> R` = 受体，产出约束 `(D, R)`；
- `Replacing Transfer` / `Additive Transfer` / 纯 `Transfer` 均被识别（触发条件为同一行**同时**
  含 `Mapping` 与 `Recipient`，见 `parse_transfer_line`）；
- `Species Tree:` / `Gene Tree:` 等段落头与 newick 行不含 `Recipient -->`，自动忽略。

## 3. 转移类别过滤（ARTra 专属）

MaxTiC 的时间约束语义与转移是加性还是替换无关，故**默认统计全部转移**
（`transfer_kind="all"`）。如研究只关心某一类别：

- `transfer_kind="replacing"`：仅统计替换转移；
- `transfer_kind="additive"`：仅统计加性转移。

（类别过滤在主进程内按行筛选后再解析，不使用并行 / 缓存，以免污染共享缓存键。）

## 4. 支持度 / 权重 / 距离 / 错误处理

同 RANGER-DTLx：块内计数 → 互反减法 → `support = count / denominator`（分母 = 采样数，即工具
声明值；无声明值时，**只有块数 > 1 且各块叶子集合完全一致**才以块数作分母，单块报告不在此列——
ecceTERA 不校验叶子集合，见
[eccetera.md](eccetera.md) §2）→ `> min_support`；
`distance = distance_from(receptor, donor)`（**自环 `donor == receptor` 例外**：写 `distance = None`，
不被 `--min-transfer-distance` 删除）；非法端点 / 过小家族被丢弃并**按原因计数上报**。
⚠️ ARTra 会打印 `Total number of optimal solutions: 24` 却只输出**一个**最优解，
因此以块数为分母在该工具上必然退化；分母不可确立时权重为**整数计数**，
`min_support` 不具 [0,1] 语义。供体端点同样取 `donor_itself`（与 ALE 相差一个物种层级）。
丢弃的可观测性与其余适配器相同：每条被丢弃的转移都计入
`dropped_by_label_miss`（数字 ID 反解仅在 ecceTERA 路径执行，本适配器
`numeric_id_resolutions` 恒为 0），转换结束打印 stderr 诊断摘要并挂在
`ConstraintSet.diagnostics`；端点命中率为 **0** 抛
`LabelMismatchError`，低于 `min_endpoint_hit_rate`（默认 `0.5`）告警（`--quiet-adapters`
只关打印，不影响计数与抛错）。

## 5. 已知限制

- 端点为物种标签，须与物种树标签一致；
- 机器学习分类置信度未纳入权重（权重仍为采样支持度）；如需以分类概率加权，可在上游
  预处理阶段生成加权约束后走文本约束路径。
- 本节为摘要：端点层级、权重口径与 diagnostics 计数等**完整语义**以规范性语义注记
  [README_semantics.md](../../src/maxtic_next/constraints/adapters/README_semantics.md) /
  [README_semantics.en.md](../../src/maxtic_next/constraints/adapters/README_semantics.en.md)
  为准。

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
