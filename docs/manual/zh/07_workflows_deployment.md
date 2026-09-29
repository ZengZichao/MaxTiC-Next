# 07 · 流程封装与部署

MaxTiC-Next 面向现代生信流水线，提供容器镜像与两大工作流引擎的封装。

> 本章的命令与仓库内 `Dockerfile`、`Singularity.def`、`workflows/` 逐字对齐。
> **注意**：本仓库开发环境未安装 `nextflow` / `snakemake` / `docker`，
> 因此引擎侧改动只做了语法与路径解析层面的静态核对，**未端到端实跑**；首次在真实集群
> 上运行时请按 7.8 的自检清单逐项确认。

## 7.1 两阶段流程模型

无论用哪种引擎，都遵循同一两阶段模型：

```
[Stage 1 · 可选] 上游工具输出 (.uml_rec / .dtl / recphyloxml / ...)
                 │  maxtic-next ... --from <tool> -o constraints.tsv
                 ▼
              统一约束文件 constraints.tsv (5 列 ALE 逗号格式)
                 │  maxtic-next species.tree constraints.tsv --ls ...
                 ▼
[Stage 2]     三输出文件（+ `--random-trees>0` 时的第 4 个）+ HTML 报告
```

- 若已有文本约束，直接从 Stage 2 开始。
- Stage 1 产物可复用、可缓存（`--ale-cache-dir` 断点续传；缓存键含物种树拓扑指纹）。
- `-o` **只在适配器模式下有意义**：文本约束模式给 `-o` 会直接报错，而不是静默忽略。
- Stage 1 与 Stage 2 是**两次进程调用**：`-o` 写出约束后即返回，不进入排序。

## 7.2 Docker

镜像基于 `python:3.11-slim`，入口为 `maxtic-next`（CLI 退出码即容器退出码）。

```bash
# 构建（在仓库根目录；名字/标签与 nextflow.config、Singularity、01 章保持一致）
docker build -t maxtic-next:0.1.0 .

# 基本运行（挂载数据目录）
docker run --rm -v "$PWD/examples":/data maxtic-next:0.1.0 \
    /data/minitree.tree /data/Cyano_CUTConstraints.tsv --seed 42

# 开启 MCMC（初步实现，收敛诊断未经验证；样本不得当作后验样本）
docker run --rm -v "$PWD/examples":/data maxtic-next:0.1.0 \
    /data/minitree.tree /data/Cyano_CUTConstraints.tsv --seed 42 \
    --mcmc --mcmc-iters 2000
```

要点：

- **镜像名必须全小写**：`docker build -t MaxTiC-Next .` 会直接报
  `invalid reference format`。统一用 `maxtic-next:0.1.0`。
- 镜像内已安装 `jinja2` + `plotly`，HTML 报告开箱可用。
- **`.dockerignore`**（新增）：`Dockerfile` 用 `COPY . /app` 拷入构建上下文，
  `.dockerignore` 已排除 `.git`、`tests`、`__pycache__`、`.pytest_cache`、`.DS_Store`、
  `benchmarks/datasets`、`docs`。保留的是 `pip install .` 与镜像内冒烟测试真正需要的东西：
  `pyproject.toml`、`src/`、`examples/`、`requirements*.txt`、`LICENSE`、`README.md`。
  因此上面的构建/运行命令不变，只是镜像更小、不带仓库元数据。
- **挂载目录写产物**：若上次运行已在挂载目录留下产物，CLI 会以**退出码 3** 拒绝覆盖
  。请加 `-f/--force`，或用 `-p /data/run2` 换前缀。
- 冒烟测试：`docker run --rm maxtic-next:0.1.0 --version` → `MaxTiC-Next 0.1.0`。

## 7.3 Singularity / Apptainer

用于 HPC 集群（无 root）。构建上下文 = `Singularity.def` 所在目录，即**仓库根**：

```bash
# 在仓库根目录构建
sudo singularity build maxtic-next_0.1.0.sif Singularity.def
# 或 apptainer build maxtic-next_0.1.0.sif Singularity.def

# 镜像内自带 examples/，可直接冒烟
singularity run maxtic-next_0.1.0.sif \
    examples/minitree.tree examples/Cyano_CUTConstraints.tsv --seed 42 --no-html

# 真实数据：挂载宿主目录
singularity run -B "$PWD/data":/data maxtic-next_0.1.0.sif \
    /data/minitree.tree /data/Cyano_CUTConstraints.tsv --seed 42
```

`%files` 是逐条白名单（与 Docker 的 `.dockerignore` 无关），必须包含
`pyproject.toml`、`README.md`（`pyproject` 的 `readme` 字段会读它）、`src/`、
`examples/`，否则 `pip install .` 或冒烟测试失败。

## 7.4 Snakemake

封装文件：`workflows/Snakefile` + `workflows/config.yaml`
（conda 环境定义：`workflows/envs/maxtic_next.yaml`）。

`config.yaml` 关键字段（与 CLI 选项同名）：

```yaml
species_tree: "examples/minitree.tree"
constraints: "examples/Cyano_CUTConstraints.tsv"
seed: 42
from_ale: false          # true 时启用 Stage 1（ALE .uml_rec -> 约束）
local_search: 0
min_transfer_distance: 0
threshold_constraints: 0
random_trees: 0          # >0 时产出第 4 个文件 *.mt.random_dist.tsv
output_style: "short"    # CLI 默认；legacy = 原版长名
use_conda: false         # true 时把 envs/maxtic_next.yaml 挂到规则的 conda: 上
outdir: "maxtic_next_out"
ale:
  rec_files:
    - "examples/adapters/ale/rec.uml_rec"
  ale_cache_dir: "ale_cache"
  ale_source: "trf"      # trf=默认（ALE 官方 MaxTiC 集成口径）/ rec=调和事件
```

运行（**两种 CWD 都已支持**）：

```bash
# 从仓库根目录
snakemake -s workflows/Snakefile --configfile workflows/config.yaml --cores 4
# 或在 workflows/ 目录内
cd workflows && snakemake --configfile config.yaml --cores 4
```

- **路径以仓库根为基准解析**：若以当前工作目录为基准，“在 `workflows/` 内运行”会去找
  `workflows/examples/...` 而失败。因此 `config.yaml` 里的相对路径一律相对**仓库根**
  求值，绝对路径原样使用。
- `from_ale: false`：只跑 Stage 2（对文本约束排序）。
- `from_ale: true`：先跑 `generate_constraints`（ALE → 约束），再排序。
- 产物（三文件 + HTML，`random_trees > 0` 时再加一个）写入 `outdir`，`rule all` 的目标
  集合与 `--output-style` 自动同步。
- 两条规则都带 `--force`：`outdir` 归流程管理器所有，重跑必须能覆盖上次产物
  （否则 CLI 的覆盖保护会以退出码 3 让任务失败）。
- `use_conda: true` 时 Snakemake 为每条规则创建 `envs/maxtic_next.yaml` 环境。该环境
  **现已包含 `jinja2` / `plotly`**——过去只装解释器，于是 HTML 报告静默降级、
  “流程跑通了但没有报告”。

> 注：仓库不随附 snakemake 引擎，需自行 `pip install snakemake`。

## 7.5 Nextflow

封装文件：`workflows/main.nf` + `workflows/nextflow.config`。

```bash
# 在仓库根目录
nextflow -C workflows/nextflow.config run workflows/main.nf \
    --species_tree examples/minitree.tree \
    --constraints examples/Cyano_CUTConstraints.tsv \
    --seed 42 --outdir results/
```

ALE 两阶段模式（物种树须与 ALE 的物种命名一致）：

```bash
nextflow -C workflows/nextflow.config run workflows/main.nf \
    --species_tree examples/adapters/ale/species.tree --from_ale \
    --ale_rec_files 'examples/adapters/ale/rec.uml_rec'
```

- 物种树与约束文件经 channel 输入自动暂存到任务目录；参数名为
  `--species_tree`（对应 `main.nf` 的 `params.species_tree`），与 CLI 一一对应透传给
  `maxtic-next`。适合云 / 集群批量。
- **`--output-style` 默认 `short`（与 CLI 一致）**，且 `RANK` 的 output 声明按
  `params.output_style` **推导**，不固定长名。`main.nf` 的 params、
  `workflows/config.yaml` 与 output 声明三处口径必须一致，否则按 CLI 默认运行时会报
  "output file not found"。改跑 `--output_style legacy` 会得到原版长名产物。
- `RANK` 的产物声明含**可选的第 4 个文件**（`--random-trees > 0` 时才有），用
  `optional: true` 表达，因此 `--random_trees 0` 不会让流程失败。
- 也带 `--force`（理由同上）。
- 用 Docker 跑：先 `docker build -t maxtic-next:0.1.0 .`，再在 `nextflow.config` 里把
  `docker.enabled` 置 true 并取消 `image = 'maxtic-next:0.1.0'` 的注释。

DSL2 结构要点：

- 两个 `process` 都声明在**顶层**，`workflow` 里只用 `if (params.from_ale)` 决定是否调用
  Stage 1。`process` 不能写在 `if` 块里，否则 DSL2 拒绝解析。
- `GENERATE_CONSTRAINTS` 的输入是 `tuple path(species), path(rec)`，调用处
  `GENERATE_CONSTRAINTS(species_ch, rec_ch)` **同序传参**：声明为
  `path species; each path(rec)`，实参一旦对调，任务起步即失败。
- 物种树被两个 Stage 共用，故用 `Channel.value(...)`（可重复订阅并自动广播）而不是
  只能订阅一次的 `Channel.fromPath` 热通道。

## 7.6 环境定义（conda）

`workflows/envs/maxtic_next.yaml` 提供 conda 环境定义（名字 `MaxTiC-Next`），供 Snakemake
的 `conda:` 指令（`use_conda: true`）或手动创建。本包尚未发布到 PyPI，创建环境后仍需安装
本包（核心零第三方依赖）：

```bash
micromamba env create -f workflows/envs/maxtic_next.yaml
micromamba run -n MaxTiC-Next pip install -e <仓库根目录>
```

环境内已声明 `jinja2 >=3.0` 与 `plotly >=5.0`（等价于 `pip install "MaxTiC-Next[report]"`）。

## 7.7 大批量部署建议

| 场景 | 建议 |
|------|------|
| 数千基因家族 → 约束 | Stage 1 开 `--ale-cache-dir`（断点续传）、`--ale-parallel process`（多核） |
| 反复调排序参数 | 固化 Stage 1 产物 `constraints.tsv`，只重跑 Stage 2 |
| 长局部搜索 | 增量打分已默认开启（加速 3.1–8.2×）+ `--checkpoint` 防中断丢失（见 08 章） |
| 需要跨机器**完全一致**的搜索结果 | `--local-search-max-iters N`（墙钟时长本身不可移植） |
| 集群可复现 | 固定 `--seed`，锁 Python/依赖版本，容器化并**固定镜像标签**（`:0.1.0` 而非 `:latest`） |
| 上游输出可疑时 | 先看 `--dry-run`（适配器模式会**真正解析**上游输出，零约束/命中率/解析失败都报 error 并退出 1），再看 stderr 的 `[tool] …` 诊断行 |
| 上游示例是压缩包 | `.gz` / gzip 流 / 单成员归档**透明解压**；上游官方 `.tgz`（如 ALE 的 `reconciliations.tgz`）可整包喂给 `CONSTRAINTS`，会自动展开为成员输入并记入 `run_metadata["archives_expanded"]`（见 05 章 5.7） |
| 需要更宽的近优邻域 | `--near-optimal-top-k K`（默认 50）放大稳健性/敏感性摘要的支持集；内存与 summary 计算随 K·n² 增长（见 03 章 3.6b） |
| 退出码 | `0` 成功 / `1` 预检 error（含 `--from` 模式下适配器解析失败）/ `2` 参数非法 / `3` 拒绝覆盖（缺 `--force`）/ `4` 前置条件不满足。非 0 可被 Snakemake/Nextflow 正确捕获重试 |

> 桌面端 `maxtic-studio`（独立仓库 `MaxTiC-Next-Studio`）复用同一张退出码表，
> 其中 **`3`** 另有一义：缺 PySide6 时打印安装指引并以 3 退出（见 03 章 3.10）。
> 工作流引擎只调 `maxtic-next`，不受该入口点影响。

## 7.8 改动后的自检清单

在没有引擎的机器上改完封装后，建议按此清单复验：

```bash
# 1) 镜像可构建、入口可用、退出码传播
docker build -t maxtic-next:0.1.0 .
docker run --rm maxtic-next:0.1.0 --version
docker run --rm maxtic-next:0.1.0 examples/minitree.tree \
    examples/Cyano_CUTConstraints.tsv --seed 42 --no-html; echo "exit=$?"

# 2) 镜像里确实没有 .git / tests / docs
docker run --rm --entrypoint sh maxtic-next:0.1.0 -c 'ls /app'

# 3) Nextflow 语法（不真跑也要过 lint）
nextflow lint -C workflows/nextflow.config workflows/main.nf
nextflow -C workflows/nextflow.config run workflows/main.nf --outdir /tmp/nf_check

# 4) Snakemake dry-run（两种 CWD 都要通）
snakemake -s workflows/Snakefile --configfile workflows/config.yaml -n -p
cd workflows && snakemake --configfile config.yaml -n -p
```

> 上一章：[06 · 上游软件配套方案](06_upstream_integration.md) ｜ 下一章：[08 · 进阶功能](08_advanced_features.md)
