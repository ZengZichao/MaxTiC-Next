# MaxTiC-Next 多软件兼容适配框架

> 本目录是**统一约束适配框架**的规格与 API 文档。它把 5 个上游调和 / 转移推断工具的
> 输出统一转换为 MaxTiC-Next 的内部约束模型 `ConstraintSet`，并提供**插件式注册**与
> **输出格式自动检测**。每个工具的详细适配规格见对应 `docs/adapters/<tool>.md`。

## 1. 为什么需要适配框架

MaxTiC 核心只消费一种抽象：一组**加权有向时间约束** `donor → receptor`（“donor 不晚于
receptor”）。但上游工具输出格式各异——从 NHX 风格的 `.uml_rec`、人类可读的调和报告，
到标准 recPhyloXML、再到预聚合的转移频率表。适配框架的职责是：**在不丢失语义的前提下**，
把这些异构输出映射到同一 `Constraint(donor, receptor, weight, metadata)` 抽象。

关键设计原则：

- **语义映射而非正则套用**：每个适配器都基于对该工具**真实输出格式的源码 / 样本求证**
  实现，而不是猜测格式（历史上 RANGER/ecceTERA 适配器曾按臆造格式解析，已按源码重写）。
- **端点合法性校验**：所有 donor / receptor 必须存在于物种树标签集合（叶子名 + 内部
  bootstrap 标签）；不匹配的端点被丢弃，**绝不静默产出错误约束**。
- **丢弃必须可见**：每一条被丢弃的转移都**按原因计数**，转换结束后向
  **stderr** 打印一行诊断摘要，并把同一份字典挂在 `ConstraintSet.diagnostics`（Python API 侧即
  `result.run_metadata["adapter_diagnostics"]`）。端点命中率为 **0** 直接抛
  `LabelMismatchError`；低于 `min_endpoint_hit_rate`（默认 0.5）告警。`--quiet-adapters`
  **只关打印**，不影响计数与抛错。
- **供体端点层级不可混用**：ALE 取 `parent(donor)`
  （`donor_endpoint_convention = "parent_of_donor"`），RANGER-DTLx / ARTra / ecceTERA / AleRax
  取**供体本身**（`"donor_itself"`）——同一事件两套口径**相差一个物种层级**。`--from-auto`
  遇混合输入按各自适配器分组解析，并告警“层级与权重口径不可通约”，不再直接相加进同一条边。
- **权重口径不可通约**：`weight_semantics` 四取之一——
  `support_over_declared_samples`（分母 = 工具**自己声明**的样本数，如 ALE 的 `N reconciled`、
  RANGER/ARTra 的 `Total number of optimal solutions: N`）、`support_over_replicate_blocks`
  （仅当各块叶子集合一致、确为同一家族重复调和时才成立）、`integer_count`（分母不可确立 →
  权重是**整数计数**，此时 `min_support` 阈值**不具 [0,1] 支持度语义**）、
  `posterior_frequency`（AleRax 上游已预聚合，本包不再归一化）。
  ⚠️ **“调和块数”不是采样数**：RANGER-DTLx / ARTra 的块头按输入**基因树**逐棵打印，ARTra
  还会打印 `Total number of optimal solutions: 24` 而只输出 1 个解。因此**不能**把"块数"当作
  支持度分母：适配器优先采用工具声明的样本数，探测不到则如实标为 `integer_count`。
- **统一元数据**：每条约束携带 `metadata = {family, support, distance}`，其中 `distance`
  为物种树拓扑距离，供 `--min-transfer-distance` 过滤。**自环例外**：`donor == receptor`
  的事件在**全部五个**适配器（RANGER-DTLx / ecceTERA / ARTra / AleRax / ALE）中一律写
  `distance = None`（无距离信息），因此按“无距离列的约束一律保留”通过过滤，`to itself`
  统计不再被 `--d` 吃掉；文本输入路径本来就是这个行为。ALE 是最后一个对齐的，
  回归测试 `tests/test_ale_selfloop.py`，实测与症状说明见 [ale.md](ale.md) 与手册 03 章 3.2.3。

## 2. 支持的工具一览

| 工具 | 输入 | 真实格式（求证来源） | 转移语义 |
|------|------|----------------------|----------|
| **ale** | `*.uml_rec` | 采样调和，NHX `T@donor->receptor` / `D@`（`constraints_from_reconciliations.py`） | `T@` 事件 donor→receptor |
| **ranger** | 调和报告文本 | `m# = LCA[..]: Transfer, Mapping --> D, ..., Recipient --> R`（`DTL-algorithm.h:1840`） | Mapping=donor, Recipient=receptor |
| **eccetera** | recPhyloXML | `branchingOut speciesLocation=D` + 子 `transferBack destinationSpecies=R`（`DTLGraph.cpp`） | 父 branchingOut=donor, 子 transferBack=receptor |
| **artra** | `output.txt` | 与 RANGER 同形，`Replacing/Additive Transfer` 前缀（`ARTra-程序/output.txt`） | Mapping=donor, Recipient=receptor |
| **alerax** | `reconciliations/summaries/*_transfers.txt` 或输出目录 | `donor receptor frequency`（`extract_families_transfer.py`） | 每行 donor→receptor，频率为权重 |

> RANGER-DTLx 与 ARTra 输出**同形**，共用解析核心 `_recon_report`；因此 `--from ranger`
> 与 `--from artra` 对同一报告产出一致结果（ARTra 额外支持按转移类别过滤）。

## 3. 统一 API

### 3.1 注册表分发

```python
from maxtic_next.tree.tree import Tree
from maxtic_next.constraints.adapters import registry

species = Tree(); species.read_newick(open("species.tree").readline())

# 显式指定工具
cset = registry.convert("ranger", species, ["FAM1.dtl"], min_family_size=0)

# 自动检测（按内容 / 扩展名）
cset = registry.convert_auto(species, ["out.recphyloxml"])

print(registry.REGISTRY.names())      # ['ale', 'alerax', 'artra', 'eccetera', 'ranger']
print(registry.detect_format("run/")) # 'alerax'
```

### 3.2 便捷函数（每工具）

```python
from maxtic_next.constraints.adapters import (
    convert_from_ale, convert_from_ranger_dtl, convert_from_eccetera,
    convert_from_artra, convert_from_alerax,
)
cset = convert_from_eccetera(species, ["FAM1.recphyloxml"], min_family_size=5)
```

### 3.3 顶层 `api.rank`（端到端）

```python
from maxtic_next import api
# 直接对上游工具输出排序（一步）
api.rank("species.tree", ["FAM1.dtl"], from_tool="ranger")
# 两阶段：仅生成统一约束文件（供 Snakemake/Nextflow 复用），不排序
api.rank("species.tree", ["run/"], from_tool="alerax", constraints_out="constraints.tsv")
```

### 3.4 通用参数

| 参数 | 含义 | 默认 | 适用 |
|------|------|------|------|
| `min_support` | 单家族内最小支持度 / 频率阈值，保留 `> min_support` | 0.05 | 全部 |
| `min_family_size` | 最小基因家族规模（叶子数），`<=` 阈值的家族被跳过 | 5 | ale/ranger/eccetera/artra（alerax 不适用） |
| `cache_dir` | 文件级 pickle 缓存目录（断点续传） | None | ale/ranger/eccetera/artra |
| `parallel` | `"process"`（进程池，绕过 GIL）/ `"thread"` | `"process"` | 全部 |
| `max_workers` | 并行 worker 数 | 自动 | 全部 |
| `output_path` | 写出统一约束（`family,donor,receptor,weight,distance`） | None | 全部 |
| `source` | ALE 约束来源 `"rec"`/`"trf"` | **`"trf"`**（`DEFAULT_SOURCE`；ALE 官方 MaxTiC 集成口径。） | ale |
| `transfer_kind` | ARTra 转移类别 `"all"`/`"replacing"`/`"additive"` | `"all"` | artra |
| `min_endpoint_hit_rate` | 端点命中率告警阈值；命中率为 0 直接报错（`<= 0` 关闭告警） | `0.5` | 全部 |
| `quiet` | 不打印 stderr 诊断摘要（计数与抛错不受影响） | `False` | 全部 |

## 4. 命令行

```bash
# 显式工具
maxtic-next species.tree FAM1.dtl        --from ranger
maxtic-next species.tree out.recphyloxml --from eccetera
maxtic-next species.tree out.txt         --from artra --artra-transfer-kind replacing
maxtic-next species.tree alerax_run/     --from alerax
maxtic-next species.tree rec1.uml_rec rec2.uml_rec --from ale   # 等价 --from-ale

# 自动检测
maxtic-next species.tree some_output      --from-auto

# 两阶段：只产约束（供流水线）
maxtic-next species.tree FAM1.dtl --from ranger -o constraints.tsv
```

`--ale-min-support` / `--ale-min-family-size` / `--ale-cache-dir` / `--ale-parallel` 对
所有调和适配器通用（历史命名沿用 `--ale-` 前缀以向后兼容）。

## 5. 自动检测规则

`detect_format(path)` 顺序（由具体到宽松）：

1. **目录**：含 `reconciliations/summaries` 或 `*_transfers.txt` → `alerax`；
2. **扩展名**：`.uml_rec` → ale；`_transfers.txt` → alerax；`.recphyloxml` → eccetera；
3. **内容嗅探**（前 64 KB）：含 `<recPhylo`/`<recGeneTree` → eccetera；含 `T@` 且
   `reconciled`/`D@` → ale；含 `Recipient -->` 且 `Replacing/Additive Transfer` → artra；
   仅含 `Recipient -->` → ranger；
4. 均不匹配 → `None`（视为通用文本约束，走原生解析）。

`detect_formats(paths)`：多路径必须检测到**同一**工具，否则抛 `ValueError`（混合输入需
显式 `--from <tool>`）。

## 6. 扩展：注册自定义适配器

```python
from maxtic_next.constraints.adapters.registry import REGISTRY, AdapterEntry

def convert_mytool(species_tree, inputs, **kwargs):
    ...  # 返回 ConstraintSet
    return cset

REGISTRY.register(AdapterEntry(
    "mytool", convert_mytool, "我的调和工具", "*.myrec 文件"))
```

注：自动检测的内容嗅探规则内置于 `detect_format`；自定义工具如需自动检测，可在应用侧
先行判定后用 `registry.convert("mytool", ...)` 显式分发。

## 7. 已知限制（通用）

- **内部节点标签对齐**：约束端点必须与用户传入 MaxTiC 的物种树标签一致。ecceTERA 对内部
  物种节点使用其自身数字 ID（见 `docs/adapters/eccetera.md`），若与物种树 bootstrap 标签
  不一致会被 `valid_labels` 过滤丢弃——这是**保守正确**（不产伪约束），但可能少召回。
- **家族规模估算**：调和报告 / recPhyloXML 的家族规模由叶子事件数估算，非精确基因树叶子数。
- **缓存信任**：文件级缓存为本地 pickle（与 ALE 一致），仅读取本次运行在 `cache_dir` 内
  自建的可信文件；勿指向不可信来源。
