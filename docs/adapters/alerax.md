# AleRax 适配器规格说明书

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

其中 `frequency` 是该转移在后验样本中的**期望频率**（浮点）。与 ALE / RANGER 的
“计数 / 采样数”不同，AleRax 已在家族内**预聚合为频率**，因此适配器**直接以频率为权重**
（不再除以块数）。`weight_semantics = posterior_frequency`。

供体端点取表文件的 donor 列本身（`donor_endpoint_convention = donor_itself`），与 ALE 的
`parent_of_donor` **相差一个物种层级**；`--from-auto` 混合 ALE 与 AleRax 输入时会分组解析并
告警“端点层级与权重口径不可通约”，不再直接相加进同一条边。

## 2. 输入形式与自动展开

`inputs` 支持三种，均由 `expand_alerax_inputs` 统一展开为转移频率文件列表：

1. AleRax 输出**根目录**（自动拼接 `reconciliations/summaries` 后 glob）；
2. `reconciliations/summaries` **目录**（直接 glob）；
3. 直接的若干转移频率**文件**。

自动检测：目录含 `reconciliations/summaries` 或 `*_transfers.txt` → `alerax`。

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
| 目录无 summaries / 无转移频率文件 | **报错**（绝不静默返回空约束集） |

## 6. 已知限制

- 转移频率文件不含基因树叶子信息，故**不做 `min_family_size` 过滤**（AleRax 已家族内
  聚合）；`min_family_size` 参数仅为接口一致性占位。
- 转移频率文件与 MaxTiC 原生空格约束格式（`donor receptor weight`）结构相同——若不经
  `--from alerax`、也不带约定后缀，将被当作通用文本约束解析（语义一致）。

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
