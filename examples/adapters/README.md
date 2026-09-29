# 多软件兼容适配 —— 各上游工具样例

本目录为 5 个上游调和 / 转移推断工具提供**真实格式**样例输入，演示如何用
`--from <tool>` / `--from-auto` 把它们的输出直接喂给 MaxTiC-Next 排序。
物种树统一用 `../minitree.tree`（内部节点标签为 bootstrap 数字，如 61/62/67）。

## 目录

| 子目录 / 文件 | 工具 | 格式 |
|---------------|------|------|
| `ranger/FAM1.dtl` | RANGER-DTLx | `Mapping --> / Recipient -->` 调和报告 |
| `eccetera/FAM1.recphyloxml` | ecceTERA | recPhyloXML（`branchingOut` / `transferBack`） |
| `artra/FAM1.txt` | ARTra | `Replacing/Additive Transfer` 报告 |
| `alerax/run/reconciliations/summaries/*_transfers.txt` | AleRax | `donor receptor freq` |
| `run_adapters.sh` | — | 一键演示全部 + 自动检测 + 两阶段 |

## 一键运行

```bash
micromamba run -n python-3.11 bash examples/adapters/run_adapters.sh
```

## 单独运行示例

```bash
# 显式工具
maxtic-next examples/minitree.tree examples/adapters/ranger/FAM1.dtl \
    --from ranger --ale-min-family-size 0

# recPhyloXML 自动检测
maxtic-next examples/minitree.tree examples/adapters/eccetera/FAM1.recphyloxml \
    --from-auto --ale-min-family-size 0

# ARTra 仅统计替换转移
maxtic-next examples/minitree.tree examples/adapters/artra/FAM1.txt \
    --from artra --artra-transfer-kind replacing --ale-min-family-size 0

# AleRax 目录（自动展开 *_transfers.txt）
maxtic-next examples/minitree.tree examples/adapters/alerax/run \
    --from alerax

# 两阶段：只产统一约束文件（供 Snakemake/Nextflow 复用）
maxtic-next examples/minitree.tree examples/adapters/ranger/FAM1.dtl \
    --from ranger -o constraints.tsv --ale-min-family-size 0
```

> 说明：这些为**小型演示样例**，故加 `--ale-min-family-size 0` 以免家族规模过滤把它们
> 滤掉；真实数据请用默认 `5` 或按需调整。各工具的适配规格与已知限制见
> `docs/adapters/<tool>.md`。
