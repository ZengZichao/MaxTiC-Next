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
- **支持度分母裁定**（`_parse_core`，三级）：① 文件里工具**自己声明**的场景数
  （`Total number of optimal solutions: N`，出现多个时取**最大值**）→
  `support_over_declared_samples`；② 无声明值时，**只要块数 > 0 就用块数**作分母
  （`support_over_replicate_blocks`），同时告警“若 ecceTERA 实际枚举了更多最优调和，本次
  支持度会被高估”——⚠️ 与 RANGER-DTLx / ARTra 不同，这里**不校验各块叶子集合是否一致**；
  ③ 一个 `<recGeneTree>` 都没解析到（块数 = 0）时才退化为 `integer_count`（**整数计数**），
  此时 `min_support` 不具 [0,1] 支持度语义。

## 3. 支持度 / 权重 / 距离

同 RANGER-DTLx：块内计数 → 互反减法 → `support = 计数 / 分母`（分母裁定见 §2）→ `> min_support`；
供体端点取 `branchingOut speciesLocation` 的**供体本身**（`donor_itself`），与 ALE 的
`parent_of_donor` 相差一个物种层级。
`distance = distance_from(receptor, donor)`；**自环 `donor == receptor` 例外**：写
`distance = None`（无距离信息），不被 `--min-transfer-distance` 删除。

## 4. 错误处理与验证

| 情形 | 处理 |
|------|------|
| 含 `<!DOCTYPE>` / `<!ENTITY>` | **拒绝解析**（防实体展开 / XXE，stdlib 无需 defusedxml） |
| XML 解析错误（**适配器路径** `_process_eccetera_worker`） | **抛 `UpstreamParseError`**，异常信息含**文件名、行列号与原始异常**——“文件损坏”与“该工具没检出转移”从此可区分，不再静默返回空 |
| XML 解析错误（**兼容纯函数** `parse_recphyloxml`，默认 `strict=False`） | 保留旧行为：返回空三元组 `({}, 0, 0)`，但向 stderr 打出含文件与行列的**响亮警告**；`strict=True` 时同样抛 `UpstreamParseError` |
| donor / receptor 解析不进物种树（纯数字 ID 已先反解，见 §5） | 该转移丢弃并计入 `dropped_by_label_miss` |
| 家族规模 `<= min_family_size` | 跳过 |
| 家族规模**探测失败**（一个 `<leaf>` 都没探到） | **跳过** `min_family_size` 过滤并告警（不清空全部转移） |

## 5. 已知限制（重要）

- ecceTERA 对**内部物种节点**在 recPhyloXML 中标注为其**自身数字节点 ID**（叶子用名字，
  `DTLGraph.cpp:2818-2821`）。编号规则求证自 `MySpeciesTree.cpp:112-136`：排序后的叶子先取
  `0..L-1`，其后各代节点**自底向上、逐层广度优先**依次编号。这些 ID 与用户传入 MaxTiC 的
  物种树 bootstrap 标签**未必一致**。
- **数字 ID 反解**：端点先按物种树标签（叶子名 + 内部 bootstrap 标签）直接匹配；匹配不上且
  端点是**纯数字**时，适配器按上述编号规则经 `_species.numeric_species_id_map` 把数字 ID
  **反解**回物种树标签，计入 `diagnostics["numeric_id_resolutions"]` 并在 stderr 告警
  （可用 `resolve_numeric_ids=False` 关闭，即旧行为：数字端点一律丢弃）。**反解仍失败**的端点
  才被 `valid_labels` 过滤丢弃（保守正确，不产伪约束，但可能少召回）；若反解后端点**全部**
  不命中，则抛 `LabelMismatchError` 而不是静默返回空集。
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
