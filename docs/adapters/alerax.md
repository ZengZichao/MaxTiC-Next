# AleRax 适配器规格说明书

> 🌐 English version: [alerax.en.md](alerax.en.md)

- **工具**：AleRax（贝叶斯基因树-物种树调和与打分），来自 `AleRax-代码`。
- **适配器**：`maxtic_next.constraints.adapters.alerax`（`AleRaxAdapter` /
  `convert_from_alerax`）。
- **工具名（注册表 / CLI）**：`alerax`（`--from alerax`）。

## 1. 真实输出格式（脚本求证）

AleRax 在输出目录内，为每个基因家族生成转移频率文件：早期版本为
`reconciliations/summaries/<family>_transfers.txt`（经
`AleRax-代码/scripts/extract_families_transfer.py` 求证）；现行版本的
`AleOptimizer.cpp` 写 summaries 处生成 `<family>_meanTransfers.txt`。
两种命名均被适配器接受（大小写不敏感）。每行为：

```
<donor_species> <receptor_species> <frequency>
```

其中 `frequency` 是该转移在后验样本中的**期望频率**（浮点，上游取值通常为 0–1，也可能是
期望计数）。与 ALE / RANGER-DTLx 的 `support = count / denominator`（分母 = 采样数，即工具
声明值）不同，AleRax 已在家族内**预聚合为频率**，因此适配器**直接以频率为权重**
（不再除以块数）。`weight_semantics = posterior_frequency`。
⚠️ 该频率**不归一化**（本包不再二次归一化）：同一条边跨家族求和后可以**超过 1**（本仓库
样本 `FAM1_transfers.txt` 第 1 行的 `0.9` 与 `FAM2_transfers.txt` 第 2 行的 `0.4` 相加即
`1.3`），此时 `min_support` 阈值**不具 [0,1] 支持度语义**，只按数值大小比较。

供体端点取表文件的 donor 列本身（`donor_endpoint_convention = donor_itself`），与 ALE 的
`parent_of_donor` **相差一个物种层级**；`--from-auto` 混合 ALE 与 AleRax 输入时会分组解析，并
**分别就端点层级、权重口径各告警一次**（`registry` 的两条独立提示：不同的供体端点层级约定 /
权重口径不一致），不再**静默地**把两个物种层级上的端点相加进同一条边——它会明确告警。

## 2. 输入形式与自动展开

`inputs` 支持三种，均由 `expand_alerax_inputs` 统一展开为转移频率文件列表：

1. AleRax 输出**根目录**（自动拼接 `reconciliations/summaries` 后 glob）；
2. `reconciliations/summaries` **目录**（直接 glob）；
3. 直接的若干转移频率**文件**。

自动检测（`detect_format`）：目录含 `reconciliations/summaries`，或目录内含
`*_transfers.txt` / `*_meanTransfers.txt` → `alerax`；单个文件按后缀判定，即以
`_transfers.txt` / `_meanTransfers.txt` 结尾 → `alerax`（**大小写不敏感**，且先剥掉 `.gz`
等压缩后缀，故 `FAM1_meanTransfers.txt.gz` 同样识别）。

## 3. 字段提取与语义映射

- 每行 `donor receptor frequency` → 约束 `(donor, receptor)`，权重 = frequency；
- 端点须在物种树标签集合内，否则丢弃；
- 家族名 = 文件名去转移频率文件后缀；
- 跨家族按边键**求和**频率（家族内先做互反减法）。

## 4. 支持度 / 权重 / 距离

- `weight = support = frequency`（无块数分母）；仅保留 `frequency > min_support`（默认 0.05）；
- `distance = distance_from(receptor, donor)`；**自环例外**：`donor == receptor` 时写
  `distance = None`（无距离信息），不被 `--min-transfer-distance` 删除。

## 5. 错误处理与验证

| 情形 | 处理 |
|------|------|
| 行字段 < 3 或频率非数字 | 跳过该行 |
| donor / receptor 不在物种树 | 该转移丢弃 |
| 端点丢弃的可观测性 | 每条被丢弃的转移**按原因计数**（`dropped_by_label_miss`；数字 ID 反解仅在 ecceTERA 路径执行，本适配器 `numeric_id_resolutions` 恒为 0），转换结束打印 stderr 诊断摘要并挂在 `ConstraintSet.diagnostics`；端点命中率为 **0** 抛 `LabelMismatchError`，低于 `min_endpoint_hit_rate`（默认 `0.5`）告警（`--quiet-adapters` 只关打印，不影响计数与抛错） |
| 目录无 summaries / 无转移频率文件 | **报错**（绝不静默返回空约束集） |

## 6. 已知限制

- 转移频率文件不含基因树叶子信息，故**不做 `min_family_size` 过滤**（AleRax 已家族内
  聚合）；`min_family_size` 参数仅为接口一致性占位。
- 转移频率文件与 MaxTiC 原生空格约束格式（`donor receptor weight`）结构相同——若不经
  `--from alerax`、也不带约定后缀，将被当作通用文本约束解析（语义一致）。
- 本节为摘要：端点层级、权重口径与 diagnostics 计数等**完整语义**以规范性语义注记
  [README_semantics.md](../../src/maxtic_next/constraints/adapters/README_semantics.md) /
  [README_semantics.en.md](../../src/maxtic_next/constraints/adapters/README_semantics.en.md)
  为准。

## 7. 用法

```bash
maxtic-next species.tree alerax_run/ --from alerax          # 根目录
maxtic-next species.tree run/reconciliations/summaries/ --from alerax
maxtic-next species.tree famA_transfers.txt famB_transfers.txt --from alerax
```

```python
from maxtic_next.constraints.adapters import convert_from_alerax
cset = convert_from_alerax(species_tree, ["alerax_run/"], min_support=0.05)
```

## 8. 样本

`tests/data/alerax_example/run/reconciliations/summaries/FAM{1,2}_transfers.txt` +
`species.tree`。回归测试见 `tests/test_alerax_adapter.py`。
