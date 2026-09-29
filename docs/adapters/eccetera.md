# ecceTERA 适配器规格说明书

- **工具**：ecceTERA（DTL 调和），来自 `ecceTERA-代码`。
- **适配器**：`maxtic_next.constraints.adapters.eccetera`（`EcceTERAAdapter` /
  `convert_from_eccetera`）。
- **工具名（注册表 / CLI）**：`eccetera`（`--from eccetera`）。

## 1. 真实输出格式（源码求证）：recPhyloXML

ecceTERA 的结构化调和输出为 **recPhyloXML** 标准（Duchemin et al. 2018），由
`ecceTERA-代码/src/DTLGraph.cpp:getRecPhyloXMLReconciliation` 生成。转移事件编码：

- **供体**节点：`<eventsRec><branchingOut speciesLocation="D"></branchingOut></eventsRec>`；
- 被转移到的**受体**子节点：`<eventsRec><transferBack destinationSpecies="R"></transferBack>
  <speciation.../leaf.../></eventsRec>`。

```xml
<recPhylo>
 <recGeneTree><phylogeny rooted="true">
  <clade><name>1</name>
   <eventsRec><branchingOut speciesLocation="61"></branchingOut></eventsRec>
   <clade><name>2</name>
    <eventsRec><leaf speciesLocation="61" geneName="g1"></leaf></eventsRec></clade>
   <clade><name>3</name>
    <eventsRec><transferBack destinationSpecies="62"></transferBack>
               <leaf speciesLocation="62" geneName="g2"></leaf></eventsRec></clade>
  </clade>
 </phylogeny></recGeneTree>
</recPhylo>
```
上例 → 约束 `61 -> 62`。

> **格式已求证**：本适配器按 ecceTERA 的真实输出（recPhyloXML）解析，字段来源已对照
> 其源码与样本文件逐项核实。

## 2. 字段提取与语义映射

- 递归遍历 `<clade>` 树并向下传递父节点的 `branchingOut speciesLocation`；
- 对含 `transferBack destinationSpecies="R"` 的子 clade，其**父** clade 的 branchingOut
  `speciesLocation="D"` 即供体 → 产出约束 `(D, R)`；
- 每个 `<recGeneTree>` 视为一个调和块；`<leaf geneName>` 计家族规模。
  ⚠️ 块数**不自动**充当支持度分母：只有各块叶子集合一致时才按
  `support_over_replicate_blocks` 归一化，否则权重为整数计数并如实标注口径。

## 3. 支持度 / 权重 / 距离

同 RANGER：块内计数 → 互反减法 → `support = 计数 / 声明样本数` → `> min_support`；
供体端点取 `branchingOut speciesLocation` 的**供体本身**（`donor_itself`），与 ALE 的
`parent_of_donor` 相差一个物种层级。
`distance = distance_from(receptor, donor)`；**自环 `donor == receptor` 例外**：写
`distance = None`（无距离信息），不被 `--min-transfer-distance` 删除。

## 4. 错误处理与验证

| 情形 | 处理 |
|------|------|
| 含 `<!DOCTYPE>` / `<!ENTITY>` | **拒绝解析**（防实体展开 / XXE，stdlib 无需 defusedxml） |
| XML 解析错误 | 返回空（不崩溃） |
| donor / receptor 不在物种树 | 该转移丢弃 |
| 家族规模 `<= min_family_size` | 跳过 |

## 5. 已知限制（重要）

- ecceTERA 对**内部物种节点**在 recPhyloXML 中标注为其**自身数字节点 ID**（叶子用名字）。
  这些 ID 与用户传入 MaxTiC 的物种树 bootstrap 标签**未必一致**——不一致的端点会被
  `valid_labels` 过滤丢弃（保守正确，不产伪约束，但可能少召回）。
- 建议：使物种树内部标签与 ecceTERA 物种命名对齐，或仅依赖叶子-叶子 / 标签一致的转移。
- ecceTERA 亦可输出 Sylvx / `.mr` 格式；本适配器只支持标准 recPhyloXML。

## 6. 用法

```bash
maxtic-next species.tree fam1.recphyloxml --from eccetera
```

```python
from maxtic_next.constraints.adapters import convert_from_eccetera
cset = convert_from_eccetera(species_tree, ["fam1.recphyloxml"], min_family_size=5)
```

## 7. 样本

`tests/data/eccetera_example/FAM1.recphyloxml` + `species.tree`。回归测试见
`tests/test_eccetera_adapter.py`。
