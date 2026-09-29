# MaxTiC-Next 使用手册 · User Manual

> 多级、中英双语的完整操作手册。涵盖软件自身的安装、用法、输入输出、Python API，
> 以及与 5 个上游调和/转移推断工具（ALE / RANGER-DTLx / ecceTERA / ARTra / AleRax）
> 的**具体配套方案**、流程封装与部署、进阶功能与故障排查。
>
> A multi-level, fully bilingual (Chinese + English) manual. It covers installation,
> usage, I/O formats, and the Python API of the software itself, plus the **concrete
> integration recipes** for the five upstream reconciliation / transfer-inference tools
> (ALE / RANGER-DTLx / ecceTERA / ARTra / AleRax), workflow packaging & deployment,
> advanced features, and troubleshooting.

---

## 中文手册（Chinese）

| 章节 | 文件 | 内容 |
|------|------|------|
| 00 · 总览 | [zh/00_index.md](zh/00_index.md) | 手册导航、术语、软件定位 |
| 01 · 安装 | [zh/01_installation.md](zh/01_installation.md) | pip / 源码 / micromamba / 依赖 |
| 02 · 快速开始 | [zh/02_quickstart.md](zh/02_quickstart.md) | 5 分钟跑通第一个例子 |
| 03 · 命令行参考 | [zh/03_cli_reference.md](zh/03_cli_reference.md) | 全部 CLI 参数逐一说明 |
| 04 · Python API | [zh/04_python_api.md](zh/04_python_api.md) | `rank` / `build_constraints` / `Result` |
| 05 · 输入输出格式 | [zh/05_io_formats.md](zh/05_io_formats.md) | 物种树、约束双格式、输出文件（含 `--random-trees` 的第 4 个）、stdout、退出码 |
| 06 · 上游软件配套方案 | [zh/06_upstream_integration.md](zh/06_upstream_integration.md) | ALE/RANGER/ecceTERA/ARTra/AleRax 端到端 |
| 07 · 流程封装与部署 | [zh/07_workflows_deployment.md](zh/07_workflows_deployment.md) | Snakemake / Nextflow / Docker / Singularity |
| 08 · 进阶功能 | [zh/08_advanced_features.md](zh/08_advanced_features.md) | 稳健性/敏感性、MCMC、剪裁、预检、检查点 |
| 09 · 常见问题与排错 | [zh/09_faq_troubleshooting.md](zh/09_faq_troubleshooting.md) | FAQ、错误信息、复现性注意事项 |

## English Manual

| Chapter | File | Content |
|---------|------|---------|
| 00 · Overview | [en/00_index.md](en/00_index.md) | Navigation, glossary, scope |
| 01 · Installation | [en/01_installation.md](en/01_installation.md) | pip / source / micromamba / deps |
| 02 · Quickstart | [en/02_quickstart.md](en/02_quickstart.md) | Run the first example in 5 minutes |
| 03 · CLI reference | [en/03_cli_reference.md](en/03_cli_reference.md) | Every CLI flag explained |
| 04 · Python API | [en/04_python_api.md](en/04_python_api.md) | `rank` / `build_constraints` / `Result` |
| 05 · Input/Output | [en/05_io_formats.md](en/05_io_formats.md) | Trees, dual constraint formats, outputs |
| 06 · Upstream integration | [en/06_upstream_integration.md](en/06_upstream_integration.md) | End-to-end recipes for the 5 tools |
| 07 · Workflows & deployment | [en/07_workflows_deployment.md](en/07_workflows_deployment.md) | Snakemake / Nextflow / Docker / Singularity |
| 08 · Advanced features | [en/08_advanced_features.md](en/08_advanced_features.md) | Robustness, MCMC, pruning, dry-run, checkpoint |
| 09 · FAQ & troubleshooting | [en/09_faq_troubleshooting.md](en/09_faq_troubleshooting.md) | FAQ, error messages, reproducibility |

---

## 阅读建议 · Reading Path

- **新用户 / New users**：01 → 02 → 03 → 05
- **要接上游工具 / Connecting upstream tools**：05 → 06 → 07
- **要嵌入流水线 / Pipeline integration**：04 → 07
- **关心不确定性 / Uncertainty analysis**：08

## 相关文档 · Related Docs

- 适配框架规格：[../adapters/README.md](../adapters/README.md) 及 `../adapters/<tool>.md`
- 顶层说明：[../../README.md](../../README.md)

## 许可与引用 · License & Citation

本软件继承原版 **CeCILL 2.1** 许可，保留原作者 Eric Tannier 与原始引用
（MaxTiC, Biorxiv doi.org/10.1101/127548）。详见仓库根目录 `LICENSE`。

This software inherits the original **CeCILL 2.1** license and preserves the original
author Eric Tannier and citation (MaxTiC, Biorxiv doi.org/10.1101/127548). See `LICENSE`.
