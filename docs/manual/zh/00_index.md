# 00 · 总览

## 0.1 MaxTiC-Next 是什么

`MaxTiC-Next` 是 Python 2 时代系统发生学工具 **MaxTiC**（Eric Tannier, Inria）的
**Python 3 忠实重写与工程化扩展**。它解决一个明确的问题：

> 给定一棵**有根物种树**，以及一组由水平基因转移（HGT）推断得到的**加权时间约束**
> （形如 “节点 A 不晚于节点 B”），寻找该物种树内部节点的一个**排序（相对年代序）**，
> 使其**最大化满足的约束累计权重**（等价于反馈弧集 / Feedback Arc Set 最小化）。

一句话：**用 HGT 提供的“谁比谁老”的信息，给物种树的内部节点排出一个尽量自洽的时间先后顺序。**

## 0.2 核心概念与术语

| 术语 | 含义 |
|------|------|
| **物种树（species tree）** | 有根 Newick 树，**内部节点标签写在 bootstrap 字段**（如 `...)61:2,...`，`61` 是标签） |
| **约束（constraint）** | 有向加权边 `donor → receptor`，含义为“donor 不晚于 receptor” |
| **供体 / 受体（donor / receptor）** | 转移的来源分支 / 目标分支所对应的物种树节点 |
| **权重（weight）** | 约束的置信度分数（如采样支持度、后验频率） |
| **MTC（最大时间一致性）** | Maximum Time Consistency：最大化不冲突约束权重 |
| **FAS（反馈弧集）** | Feedback Arc Set，被违反的约束集合；MTC 等价于最小化 FAS 权重 |
| **信息性约束（informative）** | 既非树内祖先-后代关系、也非平凡环、端点均为内部节点的约束 |
| **偏序（partial order）** | 与最优排序一致的约束集合，构成节点间的相对年代偏序 |
| **Kendall 相似度** | 输出排序与输入树（若为超度量）排序的归一化一致度，∈[0,1] |

## 0.3 三种排序启发式

MaxTiC-Next 对同一问题同时运行三种启发式，取最优：

1. **贪婪（greedy，`order_from_graph`）**：按权重降序逐条加入不产生环的约束，再拓扑排序。
2. **混合（mixing，`opt`/`mix`）**：自底向上用动态规划合并左右子树的中序。
3. **局部搜索（local search，`optimisation_locale`）**：以 Metropolis 准则对当前最优序做
   区间旋转扰动，在给定时间预算内改进（`--local-search` 打开时）。

## 0.4 软件能力速览

- **算法等价**：核心算法逐函数忠实移植原版。**默认路径（`--ls 0`）固定 `--seed` 即逐字节
  可复现**；`--ls > 0` 以墙钟为停止条件、迭代次数依机器负载而变，需再加
  `--local-search-max-iters` 才跨机器可复现（详见 03 章 3.8）。
- **双输入格式**：空格 `donor receptor [weight] [distance]` 与逗号 `family,donor,receptor,[weight],[distance]`。
- **五大上游工具适配**：ALE、RANGER-DTLx、ecceTERA、ARTra、AleRax，`--from <tool>` 一键转换，或 `--from-auto` 自动检测。
- **交互式 HTML 报告**：默认生成，含权重分布、冲突比例、节点位置稳健性、排序树可视化。
- **稳健性/敏感性分析**：局部搜索收集近优解，统计节点位置分布与成对次序频率。
- **MCMC 采样器**（可选，**初步实现**）：在与拓扑相容的线性扩展上做可逆 MH 采样；
  收敛诊断未经验证，样本**不得**当作后验样本使用（见 08 章 8.2）。
- **工程化**：现代 CLI + Python API、Docker/Singularity、Snakemake/Nextflow、预检、剪裁、
  断点续传、增量打分、输出覆盖保护（`--force`）。

## 0.5 目录结构（源码）

```
MaxTiC-Next/                  # 仓库根（https://github.com/ZengZichao/MaxTiC-Next）
├── src/maxtic_next/           # 软件包（导入名 maxtic_next）
│   ├── cli.py  api.py         # 命令行 / Python API 入口
│   ├── tree/                  # 轻量 Tree 类
│   ├── constraints/           # 约束数据结构、双格式解析、5 适配器
│   │   └── adapters/          # ale/ranger_dtl/eccetera/artra/alerax + registry
│   ├── ranking/               # 三启发式 + 主类 Ranker
│   ├── robustness/            # 稳健性摘要 + MCMC
│   ├── io/  report/           # 三文件输出 / HTML 报告
│   └── benchmarks/            # 基准套件
├── docs/                      # 本手册与适配规格（桌面端说明在 Studio 仓库）
├── examples/                  # 示例数据（minitree + Cyano + 各工具样例）
├── workflows/                 # Snakemake / Nextflow
├── tests/                     # 测试套件
├── Dockerfile  Singularity.def
└── pyproject.toml
```

## 0.6 命名约定（三层）

| 层级 | 取值 |
|------|------|
| 产品 / 显示名 | `MaxTiC-Next` |
| GitHub 仓库 | [`ZengZichao/MaxTiC-Next`](https://github.com/ZengZichao/MaxTiC-Next) |
| Python 导入包名 | `maxtic_next` |
| 控制台命令 | `maxtic-next` |

> 下一章：[01 · 安装](01_installation.md)
