# ALE 适配器规格说明书

> 🌐 English version: [ale.en.md](ale.en.md)

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
...
```

首行是**采样数声明行**：上游 ALE_undated 实际打印 `<N> reconciled G-s:`，精简样本
（如 `tests/data/ale_example/rec.uml_rec`）只有 `<N> reconciled`——适配器判定**只看第二个
token 是否恰为 `reconciled`**（`ale.py:261-262`），声明行其余部分不参与判定，故两种形态
都能识别。声明行与首棵调和树之间恒隔一行（该行不参与解析：上游输出里是空行，
本仓库样本中是家族名）。
内部节点事件标注带**前导点**（ALE_undated 恒输出 `.`+物种串），解析按 `.` 分隔事件，
如 `.59.T@61->62`、`.65`。事件类型：

- `T@donor->receptor`：一次转移（donor 供体物种、receptor 受体物种）；
- `D@species`：重复；
- 其它：物种形成（映射物种）。

`N` 是采样调和树数量，用作**支持度分母**（support = 计数 / N）。

## 2. 字段提取与语义映射

- **转移事件（trf，默认）**：`T@donor->receptor` → 约束 `(donor, receptor)`，donor 取物种树中
  该节点的**直接父**（`parent_map.get(donnor)`，即预先从物种树算好的 `label -> parent_label`
  查表；`donnor_search` 沿调和树父链上溯只用于下面的 rec 分支）。此即 **ALE 官方 MaxTiC 集成**
  （`maxtic/constraints_from_reconciliations.py` 的 `constraints_from_transfers`）规定的
  直接输入，故本包把默认值定为 `trf`。
- **调和事件（rec）**：仅对**标注中含 `T@` 事件的节点**执行（与原版门控一致）；
  受体取该节点自身映射物种（或经子节点方向匹配 + `receptor_search` 递归收集，`D@`
  节点递归合并两子树受体），donor 取 `donnor_search` 上溯结果。纯物种形成/重复节点
  不产生 rec 约束。
- **⚠️ `rec` 与 `trf` 不等价**：两者产出的约束集在**条数与权重上都不相同**（官方示例规模下
  `rec` 只有 13,999 条、总权重 5483.31，`trf` 总权重 6438.6），切换 `--ale-source`
  等于**换掉整个输入约束集**，所有下游数字随之改变。需要复现以 `rec` 为默认的旧结果时，
  请显式 `--ale-source rec`。实际生效口径会打印到 stderr，并记入
  `Result.run_metadata["adapter_diagnostics"]`。
- **端点合法性（忠实原版的宽松语义）**：受体不得为现存物种（叶子）；donor 需可在
  物种树中解析出父节点。注意 ALE 与其余四个适配器（双端点均须 ∈ 物种树标签集）
  的校验强度不同——这是对原版逐行保真的结果；混用多工具输出时建议保证上游物种
  命名与用户物种树一致。
- **无端点命中率守卫**：`convert_from_ale` **没有** `min_endpoint_hit_rate` 形参，内部以
  `finalize_adapter_diagnostics(..., 0.0, raise_on_zero_hit=False)` 收尾，因此端点命中率为
  **0** 时**不抛** `LabelMismatchError`，只向 stderr 打印一行错误级提示；其余四个适配器命中率为
  0 直接抛错。经 registry / CLI 传入的 `min_endpoint_hit_rate` 会被按签名丢弃，并记入
  `diagnostics["adapter_params_not_applicable"]`。
- **供体端点层级（与其余四个适配器相差一个物种层级）**：本适配器记
  `donor_endpoint_convention = "parent_of_donor"`（取 `parent(donor)`），而
  RANGER-DTLx / ARTra / ecceTERA / AleRax 记 `"donor_itself"`（取 `Mapping -->` 的供体本身 /
  父 clade 的 `branchingOut speciesLocation` 等）。同一转移事件在两套口径下的端点**相差一个物种层级**，
  因此 `--from-auto` 遇到混合输入时按各自适配器分组解析、并告警“层级与权重口径不可通约”，
  不再把两者直接相加进同一条边。

## 3. 支持度 / 权重 / 距离

- 单家族内按边键计数，做**互反减法**（reciprocal subtraction）消除对称噪声；
- **分母三级退化**（`_merge_payload`）：① 文件头声明的 `<N> reconciled` 采样数
  （`weight_semantics = support_over_declared_samples`）→ ② 无声明值时退为本文件**实际解析到
  的调和树数 `trees_seen`**（`support_over_replicate_blocks`，并告警）→ ③ 两者皆拿不到时为
  `integer_count`：权重按**整数计数**输出，此时 `min_support` **不具 [0,1] 支持度语义**
  （阈值退化为“至少出现一次”），`diagnostics["min_support_is_fractional"] = False`
  并在 stderr 告警；
- `support = 计数 / N`，仅保留 `support > min_support`（默认 0.05）；
- `distance = distance_from(receptor, donor)`：物种树拓扑距离（枝长视为 1）；**自环
  `donor == receptor` 例外**：写 `distance = None`（无距离信息），不被
  `--min-transfer-distance` 删除。这与 RANGER-DTLx / ecceTERA / ARTra / AleRax
  的写法一致——ALE 是该缺陷的**最后一个**对齐者，回归测试见
  `tests/test_ale_selfloop.py`，症状与实测见手册
  [03 章 3.2.3](../manual/zh/03_cli_reference.md)。

## 4. 错误处理与验证

| 情形 | 处理 |
|------|------|
| 家族规模 `叶子数 <= min_family_size`（默认 5） | 跳过该家族 |
| donor / receptor 不在物种树 | 该转移丢弃 |
| 端点命中率为 0 | 不抛错（ALE 无命中率守卫），stderr 打印错误级提示 + 记入诊断 |
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

- `tests/data/ale_example/rec_real.uml_rec` + `species.tree`：带前导点的**真实上游输出**
  （首行 `1 reconciled G-s:`，事件标注 `.59.T@61->62` / `.65`），与 §1 描述的格式一致；
- `tests/data/ale_example/rec.uml_rec` + `species.tree`：精简样本（首行 `1 reconciled`，
  事件标注 `59.T@61->62` / `X65`，无前导点）。

两种声明行形态（`rec_real.uml_rec` 的 `<N> reconciled G-s:` 与 `rec.uml_rec` 的
`<N> reconciled`）适配器都能解析；实测两个样本对同一 `species.tree` 产出的约束
`(donor, receptor, weight)` 多重集**逐条一致**——`--ale-source rec` 下均为一条
`65 -> 59, w = 1.0`，`--ale-source trf` 下均为一条 `65 -> 62, w = 1.0`（两样本还同时存在
§1 所述的前导点标注差异，该项同样不改变结果）；唯一区别是 `metadata["family"]`
取自文件名（`rec_real` / `rec`）。回归测试**分别**断言两个样本在两种 source 下的期望约束
对：`tests/test_ale_adapter.py::test_convert_endpoints_exist_in_species_tree` 与
`test_convert_transfer_source_produces_65_62` 用 `rec.uml_rec`，
`test_real_format_sample_end_to_end` 用 `rec_real.uml_rec`；测试并**没有**把两个样本互相
比对。`tests/test_ale_parallel.py` 比对的是**同一文件**顺序执行与并行执行的结果一致性。
