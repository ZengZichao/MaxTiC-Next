# MaxTiC-Next

> 🌐 English version: [README.en.md](README.en.md) ｜ 📖 完整多级手册（中英双语）：[docs/manual/README.md](docs/manual/README.md)

`MaxTiC-Next` 是 Python 2 年代发生学工具 **MaxTiC**（Eric Tannier, Inria）的 Python 3 重写与工程化扩展。
它在受控条件下（固定种子、固定 Python 版本、固定解析器）**功能等价**移植原版排序算法的基础上，
提供现代 CLI/Python API、面向 **5 个上游工具**（ALE / RANGER-DTLx / ecceTERA / ARTra / AleRax）
的插件式适配框架、稳健性/敏感性摘要与可选的 MCMC 采样器（初步实现）、交互式 HTML 报告，
以及完整的打包与流程集成（Docker/Singularity、Snakemake/Nextflow）。

> **文档导航**：安装/用法/输入输出/Python API/上游软件配套/流程部署/进阶功能/FAQ 的
> 完整中英双语手册见 [`docs/manual/`](docs/manual/README.md)（共 10 章 × 2 语言）。

## 算法等价性

一个核心设计目标是在受控条件下（固定种子、固定 Python 版本、固定解析器）**功能等价**地移植原版算法：

- `path()` / `value()` / `edgeweights()` / `mix()` / `opt()` / `order_from_graph()` /
  `optimisation_locale()` 等核心函数均按原版 `MaxTiC.py` 逐行忠实移植；
- 所有 `dict.keys()/values()/items()` 视图均用 `list()` 包裹后再 `.sort()` / `del`；
  整数除法 `/` → `//`；`print` 改为函数；`cmp` → `functools.cmp_to_key`；
- 所有随机性统一由单个 `random.Random(seed)` 实例（`maxtic_next.random_.RandomWrapper`）驱动，
  默认 `--seed 42`；
- `--min-transfer-distance`（原 `d`）**按 phylogenetic distance 列过滤**（保留 `distance > 阈值`），
  **绝不变更为按权重过滤**；无距离列的约束一律保留（此时该选项被忽略，程序会明确报出）；
- stdout 摘要与输出文件内容逐行复刻原版 `MaxTiC.py` 的顺序与措辞。

### 等价性的**准确边界**（不要引用成“逐字节等价”就完事）

1. **口径修正**：本版有意修正原版若干统计/口径缺陷（见下方“与原版的有意偏离”），
   因此 stdout 与产物**不等于**原版输出。任何“逐字节等价”的表述都必须限定为
   **“除已裁定偏离外逐字节等价”**。
2. **可复现性**：**默认路径（`--local-search 0`）固定 `--seed` 即逐字节可复现**。
   `--ls > 0` 以**墙钟时长**为停止条件，迭代次数依机器负载而变，故同种子**不保证跨机器复现**
   （程序自身会在 stdout 末尾打 WARNING 说明这点）；需要确定性请同时设置
   `--local-search-max-iters`。此外贪婪启发式在**同权重并列**时依赖 `edge` 字典插入顺序，
   所以“可复现”始终限定为“固定种子 + 同解析器 + 同 Python 版本”。
3. **验证手段分两类**（别混淆）：
   - `tests/baselines/`（`stdout.txt` + 三个产物文件快照）是**本移植版本自身在 `seed=42`
     下的确定性输出快照**，作用是**防回归**——同一份代码必须产出同一份快照。
     它**不是**原版输出，也**不**证明跨实现等价；
   - **真正的跨实现等价门禁**是 `tests/test_reference_equivalence.py`：按原版
     `MaxTiC.py` 的控制流与公式**独立重写**一份 Python 3 参考实现（不复用被测代码的
     `ReachabilityMatrix` / `EdgeBuilder` / `order_from_graph` / `mix` / `opt` / `value` /
     `RandomWrapper`），在真实蓝细菌示例上**逐字段差分比对**。核心算法或统计口径一漂移它即失败。
     它同时**显式断言**已裁定的偏离（ uninformative 分母、 偏序哨兵、 阈值后总权重、
      哨兵边不可删；并记录  与原版一致、故不算偏离）。
     Newick 解析层不在其覆盖范围内（该文件复用本项目的 `Tree`，文件内已注明）。
   - `tests/test_equivalence_stats.py` 里那条“跑原版做对比”的用例需要环境中有 `python2`
     与原版 `MaxTiC.py`；本环境没有则跳过，上面的差分测试照常运行作为兜底。

## 安装

```bash
pip install -e .            # 基本安装（核心功能，仅依赖标准库）
pip install -e ".[report]"  # 完整安装（含交互式 HTML 报告所需的 jinja2/plotly）
```

> 控制台命令名为小写的 `maxtic-next`（见 pyproject `[project.scripts]`）。
> 本仓库**只提供命令行与 Python API**，不含任何 GUI 代码与 GUI 依赖。

## 图形前端在哪里

桌面端 **MaxTiC-Next Studio**（PySide6，中英双语 + 亮暗双主题）在独立仓库
**[`MaxTiC-Next-Studio`](https://github.com/ZengZichao/MaxTiC-Next-Studio)**，独立安装、独立打包、独立发版：

```bash
pip install .                # 本仓库（核心算法）
pip install "git+https://github.com/ZengZichao/MaxTiC-Next-Studio.git"  # Studio，提供 maxtic-studio 命令
maxtic-studio
```

依赖方向是单向的：Studio 的 `pyproject.toml` 声明 `MaxTiC-Next>=0.1.0`，本仓库不引用
Studio 的任何东西，因此两个仓库放在哪里、本地目录叫什么名字都无关紧要。

免安装的 macOS 应用 `MaxTiC-Next-Studio.app`（PyInstaller 冻结，双击即可运行）由 Studio
仓库的 `packaging/build_studio_app.sh` 产出到其 `release/` 目录；使用说明见该仓库的
`docs/studio.md`。两个仓库共享同一份算法：Studio 同进程调用本仓库的 `maxtic_next.api.rank`，
不重写、不 fork 任何计算逻辑。

离线 / 不安装时，从仓库根目录设置：

```bash
export PYTHONPATH=src
python -m maxtic_next --help
python -m maxtic_next --version   # MaxTiC-Next 0.1.0
```

## 用法

```bash
# 标准库运行（无需联网装包）
PYTHONPATH=src python3 -m maxtic_next examples/minitree.tree examples/Cyano_CUTConstraints.tsv --seed 42

# 或安装后
maxtic-next examples/minitree.tree examples/Cyano_CUTConstraints.tsv --seed 42
```

### 参数（核心）

| 参数 | 说明 | 默认 |
|------|------|------|
| `species_tree` | 物种树（Newick，内部节点标签写在 bootstrap 字段；须**二叉**、内部标签与叶名各自唯一；透明支持 `.gz` / gzip 流与归档） | 必填 |
| `constraints` | 一个或多个约束文件（空格 `donor receptor [weight] [distance]` 或逗号 `family,donor,receptor,[weight],[distance]`；同样透明支持 `.gz` / 归档，上游官方 `.tgz` 示例可整包喂） | 必填 |
| `--seed` | 随机种子（驱动 mix 平局与局部搜索） | `42` |
| `--local-search` / `--ls` | 局部搜索**墙钟时长**（秒），0 表示关闭 | `0` |
| `--local-search-max-iters` | 局部搜索迭代上限（0 = 只受时长约束）；设正数以获得跨机器确定性 | `0` |
| `--temperature` / `--t` | Metropolis 温度 | `0.001` |
| `--random-type` / `--r` | 随机化类型，只接受 0/1/2 | `0` |
| `--min-transfer-distance` / `--d` | 最小转移距离阈值（按 phylogenetic distance 列过滤） | `0` |
| `--threshold-constraints` / `--ts` | 约束权重阈值比例，取值域 [0, 1] | `0.0` |
| `--random-trees` / `--rd` | 随机树采样数量（p 值用 (k+1)/(n+1) 校正） | `0` |
| `--near-optimal-top-k` | 近优解收集容量：稳健性/敏感性摘要最多考察多少个**去重**近优排序（API `top_k`；须 ≥ 1，仅 `--ls > 0` 时有意义） | `50` |
| `-f` / `--force` | 允许覆盖既有产物（默认拒绝并以退出码 3 中止） | 关 |

> 完整参数（上游适配 `--from` 与 `--min-endpoint-hit-rate`、类群剪裁、`--dry-run`、`--mcmc*`
> （含 `--mcmc-burn-in` / `--mcmc-thin` / `auto` 温度）、`--incremental`、检查点、
> `--output-style`、`-p/--output-prefix` 等）见
> [CLI 参考手册](docs/manual/zh/03_cli_reference.md)。

### Python API

```python
from maxtic_next import rank, build_constraints

result = rank("examples/minitree.tree", "examples/Cyano_CUTConstraints.tsv", seed=42)
print(result.best_source, result.similarity_to_input)   # greedy heuristic 1.0
```

## 上游工具集成

MaxTiC-Next 为 5 个上游调和/转移推断工具——**ALE / RANGER-DTLx / ecceTERA / ARTra / AleRax**——
提供适配器，经插件注册表 + 自动检测统一为加权约束：

```bash
maxtic-next species.tree FAM1.dtl        --from ranger
maxtic-next species.tree out.recphyloxml --from-auto
maxtic-next species.tree alerax_run/     --from alerax
```

- **ALE 的 `--ale-source` 默认是 `trf`**（转移事件，`parent(donor) -> receptor`），即 ALE 官方
  MaxTiC 集成 `constraints_from_transfers` 规定的口径；`rec`（调和事件）仍可选，但两者产出的
  约束集规模与量级**不同**，切换等于换掉整个输入。
- 每次转换在 **stderr** 打印一行诊断（端点约定、权重口径、采样分母、端点命中率与各类**丢弃计数**），
  并写入 `result.run_metadata["adapter_diagnostics"]`；端点命中率为 0 直接报错，低于
  `--min-endpoint-hit-rate`（默认 0.5）告警。`--quiet-adapters` 只关闭打印。

详见[上游软件配套方案](docs/manual/zh/06_upstream_integration.md)与 [`docs/adapters/`](docs/adapters/README.md) 适配规格。

## 输出

产物前缀默认为**第一个约束文件**的路径（多文件时会在 stderr 声明；可用 `-p` 覆盖）。
默认 **short** 命名（简洁、带 `.tsv` 扩展名）：

- `<constraints>.mt.informative.tsv` —— 过滤后的加权信息性约束
- `<constraints>.mt.conflicts.tsv` —— 与最优序冲突的约束
- `<constraints>.mt.partial_order.tsv` —— 偏序（含 black/green 标记）
- `<constraints>.mt.random_dist.tsv` —— 随机序的 `value similarity` 分布，
  **仅 `--random-trees > 0` 时产出**（第 4 个数据产物）

用 `--output-style legacy` 可切换回原版 MaxTiC 的长名（便于按原版文件名比对）：

- `<constraints>_MT_output_filtered_list_of_weighted_informative_constraints`
- `<constraints>_MT_output_list_of_constraints_conflicting_with_best_order`
- `<constraints>_MT_output_partial_order`
- `<constraints>_distribution_random`

> 两种风格下**内容逐字节一致**，仅文件名不同；stdout 摘要不受命名风格影响。
> 产物为**原子写入**（临时文件 + `os.replace`），且**默认拒绝覆盖**已有文件（退出码 3，
> 需 `-f/--force`）。

默认还生成交互式 HTML 报告（`--no-html` 关闭；缺 `jinja2`/`plotly` 时降级为基础模板而非报错）；
stdout 摘要包含内部节点数、约束总权重、uninformative 分项与占比、trivially conflicting、
各启发式 value、best order、排序树 Newick、与输入树的 Kendall 相似度，并在末尾追加运行期
`NOTE:` / `WARNING:` / `ERROR:` 口径声明。

### 退出码

| 码 | 含义 |
|----|------|
| `0` | 成功（含 `--dry-run` 通过） |
| `1` | `--dry-run` 存在 error 级问题（`--from` 模式下“上游输出根本解析不了”也算） |
| `2` | 参数取值域错误（argparse；含 `--near-optimal-top-k < 1`） |
| `3` | 产物已存在且未加 `--force` |
| `4` | 输入不满足算法前置条件（非二叉、标签重复、端点缺失、权重非法、适配器 0 命中、多成员归档被当作一个输入读等） |

## 与原版的有意偏离（Intentional deviations）

以下偏离**都是刻意的**，每条对应原版的一处缺陷或本版的口径裁定。它们也逐项写在
`tests/test_reference_equivalence.py` 与 `tests/test_equivalence.py` 的声明里，
由测试显式断言而非被“整体相等”掩盖。

### 修正原版缺陷（优于原版）

| 项 | 原版行为 | 本版行为 |
|----|----------|----------|
| 自环约束 `X X w` | `MaxTiC.py:201-204` **死循环** | `greedy.py` 的 `v != current` 守卫 + 归入 `to_itself` 统计；实测不再挂 |
| 单内部节点树的 Kendall | `ZeroDivisionError`（`MaxTiC.py:138`） | 返回 `1.0` 并注释说明 |
| 标签唯一 / 端点存在 | 裸 `KeyError` | 带上下文（文件名、边键、修复建议）的 `ValueError` |
| `d=MIN_TRANSFER_DIST` | `int(words[1])` **截断**距离（`MaxTiC.py:45`） | `float` |
| 阈值把约束删空 | `IndexError`（`MaxTiC.py:483`） | 守卫 + 显式 `WARNING`：0.0 的 value 只表示“无数据可检验”，不是“完美一致” |
| 2/3 列逗号误用 | `IndexError` / 错位读入 | `parsers.py` 明确拒绝并给出可操作提示 |
| 随机源 | 全局 `random`（未播种） | 单一 `RandomWrapper(seed)` 驱动全部随机性 |
| 跨来源自动检测 | 无此功能 | 新增；混合输入按工具分组解析并告警，不再静默相加 |
| `uninformative` 百分比分母 | `total_transfers + uninformative`（**双重计入**无信息权重）→ 同例报 `12%` | `total_weight` → 同例报 `13.9%` |
| 偏序产物的哨兵判定 | 写死 `edge[e] < 100000`（其 `MAX_NUMBER` 为 `1e10`）→ 静默删掉权重 ≥ 1e5 的真实信息性约束 | 以 `MAX_NUMBER` 判定 |
| `--threshold-constraints` 的总权重口径 | 在删边**之前**求和，分母含被删权重 | 在删边**之后**统计 |
| 阈值删除集合 | 可把 `MAX_NUMBER` 谱系硬约束一并删掉 → 产出**不受物种树约束**的“合法”排名 | 哨兵边排除在可删除集合之外 |
| 局部搜索后的 stdout | 打印的是**搜索前**的值与来源，与交付树矛盾 | 回写 `values` / `best_source` 并补打 `after local search … rejected` / `best found solution …` 两行 |
| `--random-trees` 的 p 值 | `k/n`，50 次采样的最小非零 p 被压成 `0.0`；且检验的是**搜索前**的序 | `(k+1)/(n+1)` 校正（Phipson & Smyth 2010），且检验**最终交付序** |
| 未播种/未定义的取值域 | `--r 3`、`--ts 5`、`--ls -5` 等静默通过 | argparse 层强制取值域，非法即退出码 2 |
| 折行 Newick、BOM、行内 `#` 注释 | 只取首行 / BOM 污染 / `float('#')` 崩 | 整份文件拼接、`utf-8-sig`、支持行内注释；错误含文件名与行号 |
| 非法权重 | `nan`/`inf`/负值一路带进目标函数 | 解析层即拒绝 |
| 非二叉物种树 | 无检查（多歧分支被静默丢弃） | 上升为**前置条件**，报错退出码 4 |
| 重复叶子名 | 无检查（节点标识塌缩） | 报错并指出重复名与父节点 |
| ALE 约束源默认值 | ALE 官方 MaxTiC 集成规定为 `constraints_from_transfers`（trf） | 默认 **`trf`**（选错口径会静默换掉整个约束集） |
| ALE 文件级缓存键 | —— | 键含**物种树拓扑指纹**（标签集相同、拓扑不同的树不再互相污染） |
| 适配器静默丢弃 | —— | 全部丢弃路径**计数上报**（stderr 一行 + `run_metadata`），端点 0 命中直接报错 |
| 输出覆盖与原子性 | `open(...,"w")` 直接覆盖、可留半截文件 | **默认拒绝覆盖**（`--force` 显式放行）+ 临时文件 + `fsync` + `os.replace` |
| 自环事件（donor == receptor）按距离被删 | 拓扑距离恒为 0 的自环会被 `d` 阈值吃掉，`to itself` 静默少算 | **五个适配器（RANGER-DTLx / ecceTERA / ARTra / AleRax / ALE）对自环一律写 `distance = None`**（无距离信息），按“无距离列一律保留”通过过滤，`to itself` 不再被 `--d` 吃掉；文本路径本来就是这个行为。ALE 是这条缺陷的最后一个对齐者，现已闭合，回归测试 `tests/test_ale_selfloop.py`，实测见手册 03 章 3.2.3 |
| 压缩输入 | 只读未压缩纯文本 | `.gz` / gzip 流按魔数、单成员 tar 归档按 `ustar` 透明解压；上游官方 `.tgz` 示例包可整包喂（自动展开为成员输入）。**多成员归档永不拼接为“一个”输入**，而是给出含 `tar -xzf` / `tar -xzOf` 的报错 |
| `--dry-run` 对上游输出 | ——（原版无预检） | 预检**真正解析**上游输出：损坏/截断的 XML、零约束、端点 0 命中一律 error 级，CLI 因此**退出 1**（旧实现在这条路径上判“通过”并退出 0） |
| 图形前端取消 | ——（原版无 GUI） | 「取消」按钮经 `stop_check` 在**迭代边界**真正打断，并交付截至当下的最优序；不再是“等跑完再丢弃结果” |

### 术语纪律（不是偏离，但必须知道）

- 局部搜索输出**只**称“基于局部搜索访问解的稳健性/敏感性摘要”，不使用“置信区间 / 后验概率”。
- 该摘要是**局部搜索访问集合的去重视图**（`is_resampling_robustness` 恒为 `False`），其支持集
  由 `--near-optimal-top-k`（API `top_k`，默认 **50**）封顶：摘要最多考察目标值最小的 50 个去重
  近优序，而链实际走过的不同序往往远多于此。调大 K 只是拓宽支持集，**不**让它变成后验；
  报告这些频率时**必须**一并给出 K 的生效值（`run_metadata["near_optimal_top_k"]`）。
- `--mcmc` 是**初步实现**：收敛诊断未经验证（实测 ESS 远小于样本量、相邻样本高重复），
  其样本**不得**当作后验样本使用。程序自身每行输出都带这句限定。

## 许可证

继承原版 **CeCILL 2.1**，保留原作者 Eric Tannier 与引用（Biorxiv doi.org/10.1101/127548）。
详见 `LICENSE`。

## 当前能力

- **核心重写**：忠实移植排序算法（贪婪 / 混合 / 局部搜索），口径缺陷已修正并有差分测试兜底。
- **接口**：`maxtic-next` 命令行与 `maxtic_next` Python API（桌面端在独立仓库
  [`MaxTiC-Next-Studio`](https://github.com/ZengZichao/MaxTiC-Next-Studio)，复用同一套 API）。
- **上游适配器**：ALE / RANGER-DTLx / ecceTERA / ARTra / AleRax，经插件注册表 + 自动检测接入，
  带端点命中率与丢弃计数诊断。
- **不确定性**：局部搜索稳健性/敏感性摘要（支持集由 `--near-optimal-top-k` 封顶，默认 50）
  + 可选 MH 采样器（初步实现，术语严格区分）。
- **报告**：交互式、无后端依赖的 HTML 报告（缺依赖时降级而非失败）。
- **输入**：约束与物种树透明支持 `.gz` / gzip 流 / tar 归档，上游官方 `.tgz` 示例包可整包喂。
- **长跑任务控制**：协作式取消（`api.rank(stop_check=...)`，Studio 的「取消」按钮即绑此钩子）、
  局部搜索迭代上限、检查点续跑。
- **分析辅助**：类群剪裁、`--dry-run` 预检、随机树 p 值（(k+1)/(n+1) 校正）、增量打分、检查点。
- **部署**：Docker/Singularity 镜像与 Snakemake/Nextflow 封装，支持并行约束生成。

## 性能说明（实际接线状态）

- **无条件启用的加速**：`ranking/reachability.py` 的 `ReachabilityMatrix`（增量传递闭包，
  O(1) 查表替代 O(V) 遍历）在**每次排序中都生效**，且已按 400/400 例证明与原版 `path()` 贪婪一致。
  它是**默认开启**的，不是可选开关。
- **默认已开启的增量打分**：`--incremental`（局部搜索每次区间旋转只算增量 O(b−a)，而非全量
  O(|E|)）现在是**默认路径**。依据：`--incremental` 与全量重算在 4 组输入 × 5 个种子 ×
  {stdout 摘要、ranked newick、best_order、values、三个 TSV} 上**逐字节一致**，实测加速
  3.1–8.2 倍；`--no-incremental` 保留原版全量口径作为逃生门，等价性由
  `tests/test_output_equivalence.py` 长期看守。
- **默认路径提速**（不改任何数值顺序）：`mix` 动态规划里的前缀成员判定由列表切片改为增量维护的
  集合，n=3199 / |E|=12796 的合成基准点单次排序中位耗时从 21.7–32.9 s 降到 **3.7–7.6 s**
  （同一输入两轮各取中位数，本机轮间噪声约 ±2 倍；口径为确定性合成树，不是真实生物数据集）。
- **可选开关**：`--checkpoint`（仅覆盖 `--local-search` 阶段的续跑）、
  ALE 的 `--ale-cache-dir`（解析阶段按家族缓存）。
- `tree/cache.py` 的 `PathCache` **不在排序热路径上**（保留导出是为了让公开 import 面
  保持稳定）。它**不能**经 `--incremental` / `--checkpoint` “启用”——那两个开关
  指向的是上面两条真实路径。相关说明已在其模块 docstring 中写明。

## 未来工作

- 在严格的等价保持回归控制下，评估把 `PathCache` 一类的树拓扑缓存真正接入热路径的收益。
- 更广的上游工具真实样本覆盖；跨版本（Python 2 原版 vs 本版）容器化等价验证进入 CI
  （当前由 `tests/test_reference_equivalence.py` 的独立参考实现差分作为门禁）。
- MCMC 的收敛诊断与自适应温度（在其通过诊断之前，`--mcmc` 保持“初步实现”的口径）。

## 已知限制（Known Limitations）

- **`--ls > 0` 不可跨机器复现**（除非同时设 `--local-search-max-iters`）：停止条件是墙钟时长。
- **`--checkpoint` 仅覆盖局部搜索阶段**：检查点/续跑作用于 `--local-search` 的 Metropolis 搜索；
  ALE 解析阶段另有按基因家族的文件级缓存（断点续传），不受此限制。
- **物种树必须二叉、标签唯一**：本版不推广 `opt`/`maximum_distance`/`random_order` 到多叉，
  而是把二叉假设作为前置条件强制（见“与原版的有意偏离”）。
- **不同适配器的权重口径不可通约**：ALE 是 [0,1] 支持度，RANGER-DTLx / ARTra 在报告未声明
  样本数时是**整数计数**，AleRax 是上游预聚合的后验频率。混用时 `--ale-min-support` 的语义
  随之改变，诊断行会明确告诉你当前是哪一种。
