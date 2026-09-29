# 06 · 上游软件配套方案

本章是 MaxTiC-Next 与 **5 个上游调和 / 转移推断工具**端到端配套的操作手册。核心思路：

> 上游工具（推断 HGT 转移事件）→ **适配器**统一为 `Constraint(donor, receptor, weight, metadata)`
> → MaxTiC-Next 排序。你既可以“一步到位”（`--from <tool>` 直接排序），也可以“两阶段”
> （先 `-o` 生成统一约束文件，再排序）。

## 6.0 统一心智模型

所有适配器都产出同构的 `Constraint(donor, receptor, weight, metadata={family, support, distance})`：

- **donor → receptor**：语义恒为“donor 不晚于 receptor”；
- **端点合法性**：donor/receptor 必须存在于**你传给 MaxTiC 的物种树**标签集合（内部 bootstrap 标签 + 叶名），否则该约束被丢弃（**绝不静默产伪约束**）；
- **权重/支持度**：多为“采样计数 / 采样数”或后验频率；
- **距离**：`distance_from(receptor, donor)` 的物种树拓扑距离，供 `--d` 过滤。

通用参数（对所有适配器生效，历史沿用 `--ale-` 前缀）：

| 参数 | 含义 | 默认 |
|------|------|------|
| `--ale-min-support` | 保留 `support > 阈值` | `0.05` |
| `--ale-min-family-size` | 家族叶子数 `<=` 阈值则跳过（alerax 不适用） | `5` |
| `--ale-cache-dir` | 文件级 pickle 缓存（断点续传） | 无 |
| `--ale-parallel` | `process`（默认）/ `thread` | `process` |

> **小样本提示**：演示用的小家族容易被 `--ale-min-family-size 5` 滤掉，可临时设 `--ale-min-family-size 0`；真实数据用默认或按需调。

自动检测（`--from-auto`）规则概要：目录含 `reconciliations/summaries` 或 `*_transfers.txt` → alerax；
`.uml_rec` → ale；`.recphyloxml` / `<recPhylo>` → eccetera；含 `Recipient -->` 且
`Replacing/Additive Transfer` → artra；仅 `Recipient -->` → ranger。多文件必须同一工具，否则需显式 `--from`。

---

## 6.1 ALE（`ALEml_undated` / `ALEmcmc_undated`）

**工具定位**：Amalgamated Likelihood Estimation，对每个基因家族做采样调和，输出 `.uml_rec`。

**真实格式**：`<N> reconciled trees:` 头 + N 棵 NHX 风格采样调和树；事件标注在节点标签上，
`T@donor->receptor`（转移）、`D@species`（重复）。`N` 为支持度分母。

**上游命令（示意）**：

```bash
# 1. 观测 & 调和（每个基因家族）
ALEobserve gene_family_XXX.treelist
ALEml_undated species_tree.nwk gene_family_XXX.treelist.ale
# 得到 gene_family_XXX.treelist.ale.uml_rec

# 2. 收集所有 .uml_rec
ls *.uml_rec > all_rec_files    # 或直接把它们作为参数列出
```

**MaxTiC-Next 配套（一步到位）**：

```bash
maxtic-next species_tree.nwk fam1.ale.uml_rec fam2.ale.uml_rec ... \
    --from ale --ale-source rec \
    --ale-min-support 0.05 --ale-min-family-size 5 \
    --ale-cache-dir ale_cache/
```

**两阶段（推荐用于大批量 / 流水线）**：

```bash
# 阶段 1：仅生成统一约束文件（可复用、可缓存续传）
maxtic-next species_tree.nwk *.uml_rec --from ale \
    --ale-cache-dir ale_cache/ -o constraints.tsv
# 阶段 2：对约束排序
maxtic-next species_tree.nwk constraints.tsv --ls 180
```

**关键点**：
- `--ale-source` 默认 **`trf`**：每个 `T@donor->receptor` 事件产出
  `parent(donor) -> receptor`，即 **ALE 官方 MaxTiC 集成** `constraints_from_transfers`
  的口径。`rec` 取自 `.uml_rec` 的调和记录块（该块在上游脚本里位于 `T@` 门控体内，
  纯物种形成/重复节点不产生约束）。**两者不是等价写法**：约束集的规模与权重量级都不同，
  切换即换掉整个输入，需复现旧结果请显式 `--ale-source rec`。实际口径会打到 stderr。
- 供体端点约定为 **`parent_of_donor`**（供体的父节点）——与 RANGER-DTLx / ARTra /
  ecceTERA / AleRax 的 `donor_itself`（供体本身）**相差一个物种层级**，
  `--from-auto` 混合输入时会分组解析并告警，不再静默相加。
- 受体不得为叶子（现存物种）；donor/receptor 须在物种树标签集内。
- 单家族内做**互反减法**消除对称噪声；`support = 计数 / N`，N 取报告自己声明的
  `N reconciled`。支持度低于 `--ale-min-support`（默认 0.05）的约束被丢弃并计数。
- `--ale-cache-dir` 的缓存键含**物种树拓扑指纹**：标签集相同、拓扑不同的
  两棵树不会复用同一缓存条目。
- 详见 [../adapters/ale.md](../../adapters/ale.md)。

---

## 6.2 RANGER-DTLx

**工具定位**：DTL 调和 + AggregateRanger，输出**人类可读**调和报告。

**真实格式**（源码 `DTL-algorithm.h:1840`）：逐基因节点打印，转移行形如

```
 = LCA[g1, g4]: Transfer, Mapping --> 61, Edge, Parent = 65, Recipient --> 62
```

`Mapping --> D` = 供体、`Recipient --> R` = 受体 → 约束 `(D, R)`。一个文件可含多个
`Reconciliation for Gene Tree N` 块，块数为支持度分母。

**上游命令（示意）**：

```bash
Ranger-DTL -i gene_family_XXX.nwk -s species_tree.nwk -o FAM_XXX.dtl
# 多最优解 / 多样本可合并到同一文件的多个块，或分文件由适配器跨文件求和
```

**MaxTiC-Next 配套**：

```bash
# 一步到位
maxtic-next species_tree.nwk FAM1.dtl FAM2.dtl --from ranger \
    --ale-min-family-size 5

# 两阶段
maxtic-next species_tree.nwk *.dtl --from ranger -o constraints.tsv
maxtic-next species_tree.nwk constraints.tsv --ls 180
```

**关键点**：内部节点常为数字标签（如 `61`），须与物种树 bootstrap 标签一致。
详见 [../adapters/ranger_dtlx.md](../../adapters/ranger_dtlx.md)。

---

## 6.3 ecceTERA

**工具定位**：DTL 调和，输出标准 **recPhyloXML**（Duchemin et al. 2018）。

**真实格式**（源码 `DTLGraph.cpp`）：供体节点 `<branchingOut speciesLocation="D">`，
被转移到的受体子节点 `<transferBack destinationSpecies="R">` → 约束 `(D, R)`。

**上游命令（示意）**：

```bash
ecceTERA species.tree=species_tree.nwk gene.trees=gene_family_XXX.nwk \
    output.dir=out/ recPhyloXML.reconciliation=true
# 得到 out/..._recPhyloXML.xml（recPhyloXML）
```

**MaxTiC-Next 配套**：

```bash
maxtic-next species_tree.nwk fam1.recphyloxml fam2.recphyloxml --from eccetera
# 或自动检测（recPhyloXML 会被识别为 eccetera）
maxtic-next species_tree.nwk out.recphyloxml --from-auto
```

**关键限制（重要）**：ecceTERA 对**内部物种节点**用其自身**数字 ID**标注（叶子用名字），
这些 ID 未必等于你物种树的 bootstrap 标签——不一致的端点会被保守丢弃（不产伪约束，但可能少召回）。
建议使物种树内部标签与 ecceTERA 命名对齐。含 `<!DOCTYPE>`/`<!ENTITY>` 的 XML 会被**拒绝解析**（防 XXE）。
详见 [../adapters/eccetera.md](../../adapters/eccetera.md)。

---

## 6.4 ARTra（Additive / Replacing Transfer 分类）

**工具定位**：在 RANGER 风格调和上把每个转移分类为**加性（Additive）**或**替换（Replacing）**。

**真实格式**：与 RANGER 调和报告**同形**，仅事件名多 `Replacing`/`Additive` 前缀：

```
 = LCA[H117, H14]: Replacing Transfer, Mapping --> H15, Recipient --> H125
 = LCA[H159, H17]: Additive Transfer, Mapping --> H17, Recipient --> H159
```

**MaxTiC-Next 配套**：

```bash
# 统计全部转移（默认）
maxtic-next species_tree.nwk FAM1.txt --from artra --artra-transfer-kind all
# 仅替换转移
maxtic-next species_tree.nwk FAM1.txt --from artra --artra-transfer-kind replacing
# 仅加性转移
maxtic-next species_tree.nwk FAM1.txt --from artra --artra-transfer-kind additive
```

**关键点**：MaxTiC 的时间约束语义与转移是加性/替换无关，故默认统计全部；
类别过滤只是让你聚焦某一类。与 RANGER 共用解析核心，同一报告结果一致。
详见 [../adapters/artra.md](../../adapters/artra.md)。

---

## 6.5 AleRax

**工具定位**：贝叶斯基因树-物种树调和与打分，在家族内**预聚合**转移为期望频率。

**真实格式**（脚本 `extract_families_transfer.py`）：输出目录内每个家族一份
`reconciliations/summaries/<family>_transfers.txt`，每行 `donor receptor frequency`。
频率已是期望值，**直接作为权重**（不再除以块数）。

**MaxTiC-Next 配套**（输入支持三种形式，自动展开）：

```bash
# 1. AleRax 输出根目录（自动拼 reconciliations/summaries）
maxtic-next species_tree.nwk alerax_run/ --from alerax
# 2. summaries 目录
maxtic-next species_tree.nwk alerax_run/reconciliations/summaries/ --from alerax
# 3. 直接给若干 *_transfers.txt
maxtic-next species_tree.nwk famA_transfers.txt famB_transfers.txt --from alerax
```

**关键点**：`_transfers.txt` 不含叶子信息，故**不做 `min_family_size` 过滤**；
仅保留 `frequency > --ale-min-support`（默认 0.05）；跨家族按边键求和。
详见 [../adapters/alerax.md](../../adapters/alerax.md)。

---

## 6.6 五工具对照速查

| 工具 | 输入 | 供体/受体来源 | 权重 | 家族过滤 | CLI |
|------|------|---------------|------|----------|-----|
| ALE | `*.uml_rec` | `T@D->R` / 调和事件集合 | 计数/N | ✔ | `--from ale` |
| RANGER-DTLx | 调和报告文本 | `Mapping-->D` / `Recipient-->R` | 计数/块数 | ✔ | `--from ranger` |
| ecceTERA | recPhyloXML | `branchingOut` / `transferBack` | 计数/块数 | ✔ | `--from eccetera` |
| ARTra | `output.txt` | 同 RANGER | 计数/块数 | ✔ | `--from artra` |
| AleRax | `*_transfers.txt` | 每行 `donor receptor` | 频率 | ✘ | `--from alerax` |

## 6.7 端到端最佳实践

1. **先预检**：`--dry-run` 检查物种树/标签唯一性/端点合法性（尤其 ecceTERA 的 ID 对齐问题）。
2. **两阶段**：大批量家族先 `-o constraints.tsv` 生成统一约束（配 `--ale-cache-dir` 断点续传），
   再对约束反复调参排序（`--ls`/`--ts`/`--d`），避免每次重解析。
3. **自动检测**：不确定格式用 `--from-auto`；但**混合工具输入**必须显式 `--from`。
4. **距离过滤**：仅在约束确有 distance 列（ALE 5 列 / 空格 4 列）时用 `--d`。
5. **可复现**：固定 `--seed`，记录版本与参数（HTML 报告与 stdout 头部均含运行元数据）。

## 6.8 自定义适配器（扩展）

```python
from maxtic_next.constraints.adapters.registry import REGISTRY, AdapterEntry

def convert_mytool(species_tree, inputs, **kwargs):
    ...  # 解析你的工具输出，返回 ConstraintSet
    return cset

REGISTRY.register(AdapterEntry("mytool", convert_mytool, "我的调和工具", "*.myrec"))
```

> 上一章：[05 · 输入输出格式](05_io_formats.md) ｜ 下一章：[07 · 流程封装与部署](07_workflows_deployment.md)
