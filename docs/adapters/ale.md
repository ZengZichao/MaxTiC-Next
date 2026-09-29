# ALE 适配器规格说明书

- **工具**：ALE_undated（AmalgamatedLikelihoodEstimation），上游独立项目。
- **适配器**：`maxtic_next.constraints.adapters.ale`（`ALEAdapter` / `convert_from_ale`）。
- **工具名（注册表 / CLI）**：`ale`（`--from ale`，向后兼容 `--from-ale`）。
- **来源**：忠实移植原版 MaxTiC 自带的 `constraints_from_reconciliations.py`
  （Python 2）的控制流，包括 **rec 约束仅由带 `T@` 标注的节点产生**（原版
  "2/ calculer la contrainte selon les reconciliations" 块位于 `if e[:2]=="T@"`
  门控体内，每个 T@ 事件执行一次）。

## 1. 真实输出格式

ALE_undated 对每个基因家族产出 `<family>.ale.uml_rec` 采样调和文件。结构：

```
<N> reconciled G-s:
<空行>
(gene_tree_1 with NHX-style event annotations);
(gene_tree_2 ...);
...（共 N 棵采样调和树）
```

内部节点事件标注带**前导点**（ALE_undated 恒输出 `.`+物种串），解析按
`.` 分隔事件，如 `.59.T@61->62`、`.65`。事件类型：

- `T@donor->receptor`：一次转移（donor 供体物种、receptor 受体物种）；
- `D@species`：重复；
- 其它：物种形成（映射物种）。

`N` 是采样调和树数量，用作**支持度分母**（support = 计数 / N）。

## 2. 字段提取与语义映射

- **转移事件（trf，默认）**：`T@donor->receptor` → 约束 `(donor, receptor)`，donor 取物种树中
  该节点的父（`donnor_search` 沿父链上溯首个非 T@/D@ 事件）。此即 **ALE 官方 MaxTiC 集成**
  （`maxtic/constraints_from_reconciliations.py` 的 `constraints_from_transfers`）规定的
  直接输入，故本包把默认值定为 `trf`。
- **调和事件（rec）**：仅对**标注中含 `T@` 事件的节点**执行（与原版门控一致）；
  受体取该节点自身映射物种（或经子节点方向匹配 + `receptor_search` 递归收集，`D@`
  节点递归合并两子树受体），donor 取 `donnor_search` 上溯结果。纯物种形成/重复节点
  不产生 rec 约束。
- **⚠️ `rec` 与 `trf` 不等价**：两者产出的约束集在**规模与权重量级上都不同**（官方示例规模
  下边数与总权重相差数倍），切换 `--ale-source` 等于**换掉整个输入约束集**，所有下游数字随之
  改变。需要复现以 `rec` 为默认的旧结果时，请显式 `--ale-source rec`。实际生效口径会打印到
  stderr，并记入 `Result.run_metadata["adapter_diagnostics"]`。
- **端点合法性（忠实原版的宽松语义）**：受体不得为现存物种（叶子）；donor 需可在
  物种树中解析出父节点。注意 ALE 与其余四个适配器（双端点均须 ∈ 物种树标签集）
  的校验强度不同——这是对原版逐行保真的结果；混用多工具输出时建议保证上游物种
  命名与用户物种树一致。
- **供体端点层级（与其余四个适配器相差一代）**：本适配器记
  `donor_endpoint_convention = "parent_of_donor"`（取 `parent(donor)`），而
  RANGER-DTLx / ARTra / ecceTERA / AleRax 记 `"donor_itself"`（取 `Mapping -->` 的供体本身 /
  `branchingOut` 的父节点等）。同一转移事件在两套口径下的端点**相差一个物种层级**，
  因此 `--from-auto` 遇到混合输入时按各自适配器分组解析、并告警“层级与权重口径不可通约”，
  不再把两者直接相加进同一条边。

## 3. 支持度 / 权重 / 距离

- 单家族内按边键计数，做**互反减法**（reciprocal subtraction）消除对称噪声；
- `support = 计数 / N`，仅保留 `support > min_support`（默认 0.05）；
- `distance = distance_from(receptor, donor)`：物种树拓扑距离（枝长视为 1）；**自环
  `donor == receptor` 例外**：写 `distance = None`（无距离信息），不被
  `--min-transfer-distance` 删除。这与 RANGER-DTLx / ecceTERA / ARTra / AleRax
  的写法一致——ALE 是该缺陷的**最后一个**对齐者，回归测试见
  `tests/test_ale_selfloop.py`，症状与实测见手册 03 章 3.2.3。

## 4. 错误处理与验证

| 情形 | 处理 |
|------|------|
| 家族规模 `叶子数 <= min_family_size`（默认 5） | 跳过该家族 |
| donor / receptor 不在物种树 | 该转移丢弃 |
| 单文件解析异常（并行） | 自动回退顺序重算该文件 |
| 缓存命中（`cache_dir`） | 直接复用（断点续传） |

## 5. 已知限制

- ALE 的 donor/receptor 为物种标签，须与物种树标签一致（本工具通常同源，故一致）。
- `.uml_rec` 事件标注格式对空白 / 分隔敏感，已按原版逐字节移植；非 ALE_undated 变体格式
  需另行适配。

## 6. 用法

```bash
maxtic-next species.tree fam1.ale.uml_rec fam2.ale.uml_rec --from ale \
    --ale-source rec --ale-min-support 0.05 --ale-min-family-size 5
```

```python
from maxtic_next.constraints.adapters import convert_from_ale
# source 的默认值是 "trf"（DEFAULT_SOURCE，ALE 官方 MaxTiC 集成口径）；
# 下面显式写 source="rec" 是为了演示"复现旧口径"这一用法，不是默认行为。
cset = convert_from_ale(species_tree, ["fam1.ale.uml_rec"], source="rec",
                        cache_dir="cache/", parallel="process")
print(cset.diagnostics["donor_endpoint_convention"])  # 'parent_of_donor'
```

## 7. 样本

`tests/data/ale_example/rec.uml_rec` + `species.tree`。回归测试见
`tests/test_ale_adapter.py` / `tests/test_ale_parallel.py`。
