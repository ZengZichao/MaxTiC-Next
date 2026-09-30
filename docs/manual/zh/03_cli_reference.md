# 03 · 命令行参考

```
maxtic-next SPECIES_TREE CONSTRAINTS [CONSTRAINTS ...] [选项]
```

- `SPECIES_TREE`：物种树文件路径（Newick，内部节点标签在 bootstrap 字段）。**必填**。
- `CONSTRAINTS`：一个或多个约束文件路径。文本模式为约束文件；`--from <tool>` 模式下为对应
  上游工具的输出文件/目录。**必填**（至少一个）。

> 两个位置参数都**透明接受** `.gz` / gzip 流与 tar 归档（`.tar.gz` / `.tgz` / `.tar`），
> 不需要先手工解压；支持范围与边界见 [05 章 5.7](05_io_formats.md)。

> 核心参数等价原版 `MaxTiC.py` 的裸 `sys.argv` 语义；长短选项等价（如 `--local-search` = `--ls`）。
> 下表逐条对照 `src/maxtic_next/cli.py` 的 `build_parser()`（默认值取自
> `src/maxtic_next/config.py`），并与 `maxtic-next --help` 一致。

## 3.0 元信息

| 参数 | 说明 |
|------|------|
| `--version` | 打印 `MaxTiC-Next <版本>` 后退出 0（当前 `MaxTiC-Next 0.1.1`） |
| `-h`, `--help` | 打印完整参数列表 |

## 3.1 核心排序参数（等价原版）

| 参数 | 别名 | 类型 / 取值域 | 默认 | 说明 |
|------|------|--------------|------|------|
| `--seed` | — | int ≥ 0 | `42` | 随机种子，驱动 mix 平局与局部搜索 |
| `--local-search` | `--ls` | float ≥ 0 | `0.0` | 局部搜索**墙钟时长**（秒），0 = 关闭；打开后触发稳健性/敏感性摘要 |
| `--temperature` | `--t` | float > 0 | `0.001` | 局部搜索的 Metropolis 温度（接受劣解的尺度） |
| `--random-type` | `--r` | int ∈ {0,1,2} | `0` | 随机化类型（`0` 保持数据原样；`2` 完全随机化节点） |
| `--min-transfer-distance` | `--d` | float ≥ 0 | `0` | 最小转移距离阈值，**按 phylogenetic distance 列过滤**，保留 `distance > 阈值` |
| `--threshold-constraints` | `--ts` | float ∈ **[0, 1]** | `0.0` | 约束权重阈值比例，按权重升序剔除累计占比最低的一部分边 |
| `--random-trees` | `--rd` | int ≥ 0 | `0` | 随机树采样数量，>0 时计算 value / similarity 的经验 p 值并写分布文件 |

> **取值域在 argparse 层强制**：`--r 3`、`--ts 5`、`--ls -5`、`--t -1`、
> `--rd -1` 之类的旧写法过去被静默接受（等价于悄悄改了口径），现在一律**退出码 2** 并报出
> 可操作的提示。`--ts` 的语义是“删除累计权重占比最低的约束”，故它是**比例**而非绝对值：
> `--ts 5` 不再合法，`--ts 0.5` 表示删除占比 50% 的低权重约束。

> ⚠️ **`--d` 的真实语义**：它**只**按系统发育距离过滤，
> **绝不**按权重过滤；而“没有距离列”的约束**一律保留**，不是“全部丢弃”。
> 因此对 3 列空格 / 4 列逗号这类**无距离列**的输入，`--d 5` 一条也不会删。
> 这种情况现在是 **error 级**（不再是轻描淡写的 warning），stdout 末尾打印：
>
> ```
> ERROR: --min-transfer-distance=5.0 被忽略（123/123 条约束不含 phylogenetic distance 列，……）
> ```
>
> 并退出 0（它是口径声明，不改变计算）；`--dry-run` 同口径报 error 并退出 1。
> 需要按距离过滤请提供 4 列（空格格式）或 5 列（ALE 逗号格式）输入；不需要请去掉该选项。

## 3.2 上游工具适配参数

| 参数 | 取值 | 默认 | 说明 |
|------|------|------|------|
| `--from` | `ale`/`ranger`/`eccetera`/`artra`/`alerax`/`auto` | 无 | 把 CONSTRAINTS 视为该工具输出，经注册表转为统一约束 |
| `--from-auto` | — | 关 | 等价 `--from auto`：按内容/扩展名自动检测上游格式 |
| `--from-ale` | — | 关 | 向后兼容开关，等价 `--from ale`（CONSTRAINTS 为 `.uml_rec` 列表） |
| `--ale-min-support` | float ∈ [0,1] | `0.05` | 单家族内最小支持度/频率阈值，保留 `> 阈值`（所有调和适配器通用） |
| `--ale-min-family-size` | int ≥ 0 | `5` | 最小基因家族规模（调和树叶子数），`<=` 阈值的家族被跳过（alerax 不适用） |
| `--ale-cache-dir` | 路径 | 无 | 基因家族文件级 pickle 缓存目录（断点续传） |
| `--ale-source` | `rec`/`trf` | **`trf`** | ALE 约束来源，见下方专门说明 |
| `--ale-parallel` | `process`/`thread` | `process` | 并行解析模式：进程池（默认，绕 GIL）/ 线程池 |
| `--artra-transfer-kind` | `all`/`replacing`/`additive` | `all` | ARTra 统计的转移类别 |
| `--min-endpoint-hit-rate` | float ∈ [0,1] | `0.5` | 适配器**端点命中率**阈值：低于该比例向 stderr 警告；命中率为 0 直接报错。设为 `0` 关闭告警（0 命中仍报错） |
| `--quiet-adapters` | — | 关 | 不打印适配器在 stderr 上的例行口径声明与诊断摘要；**只影响打印**，诊断量仍写入 `Result.run_metadata["adapter_diagnostics"]`，错误照旧抛出 |

> 历史命名沿用 `--ale-` 前缀，但 `--ale-min-support` / `--ale-min-family-size` /
> `--ale-cache-dir` / `--min-endpoint-hit-rate` / `--quiet-adapters` 对
> `ranger/eccetera/artra/alerax` 同样生效。

### 3.2.1 `--ale-source`：默认已是 `trf`

| 取值 | 生成方式 | 语义 |
|------|----------|------|
| **`trf`（默认）** | 按**转移事件**生成，等价 ALE 官方 MaxTiC 集成脚本的 `constraints_from_transfers` | “供体分支的**父节点**须早于受体分支的**子节点**” |
| `rec` | 按**调和事件**（`.uml_rec` 的 reconciliation 记录）生成 | 约束集规模与权重量级与 `trf` **不同** |

- 默认取 `trf`，因为 **ALE 自带**的 MaxTiC 集成文档规定的直接输入就是
  `constraints_from_transfers`（trf）。这不是“两种等价写法”：**切换即换掉整个输入约束集**
  （官方示例规模下边数与总权重差异达数倍），所有下游数字随之改变。
- 实际生效口径会打印到 stderr，例如：
  `[ale] 约束口径：source=trf（……ALE 官方 MaxTiC 集成口径），供体端点约定=parent_of_donor，产出边数=…，总权重=…`
- 需要 `rec`（调和事件）口径时显式加 `--ale-source rec`。

### 3.2.2 适配器诊断行（stderr）

每次适配器转换结束打印一行，字段固定：

```
[tool] donor 端点约定=<parent_of_donor|donor_itself> 权重口径=<support_over_declared_samples|support_over_replicate_blocks|integer_count|posterior_frequency> 分母=<N|不可确立> 端点命中率=NN.N% | 文件=N 块=N 树=N 转移=N 保留=N 丢弃[标签不匹配=N, 支持度=N, 家族规模=N, 家族规模探测失败=N, 自环=N]
```

- `donor 端点约定`：ALE 取 `parent_of_donor`（供体的父节点），RANGER-DTLx/ARTra/ecceTERA/
  AleRax 取 `donor_itself`（供体本身）——**相差一个物种层级**，`--from-auto` 混合输入时
  会另行告警，不再静默相加。
- `权重口径`：ALE 是 [0,1] 支持度；RANGER/ARTra 在报告未声明样本数时是**整数计数**，
  此时 `--ale-min-support` 阈值不具 [0,1] 语义。不同口径的权重**不可直接比较**。
- `自环=N`：`donor == receptor`（供体与受体是同一谱系）的事件条数。这个数字是**计数**，
  不是“丢了 N 条”：自环约束**仍然保留**，并把自己的权重计入 stdout 摘要里原版口径的
  **`to itself`** 分项（，见 3.2.3）。五个适配器（含 ALE）一律如此。
- 同样的字典在 Python 侧可取：`result.run_metadata["adapter_diagnostics"]`。

### 3.2.3 自环转移与 `--d`

供体与受体同为一个谱系时，物种树上的拓扑距离恒为 0。若适配器把这个 0 写进 distance 列，
`--min-transfer-distance` 就会以“距离太近”为由把它删掉，于是这类事件从 `to itself`
统计里**静默消失**。现在 **五个适配器 RANGER-DTLx / ecceTERA / ARTra / AleRax / ALE**
对自环一律写 **`distance = None`**（= 无距离信息），按 `--d` 的既有口径“无距离列的约束
一律保留”，它们因此**不再被 `--d` 吃掉**。ALE 是最后一个对齐的，
回归测试为 `tests/test_ale_selfloop.py`。

实测一（RANGER-DTLx 输入：用 `examples/adapters/ale/species.tree` + 一份把 `` 的
Recipient 改写成其自身 Mapping（62 → 62）的 RANGER-DTLx 报告）：

```
$ maxtic-next species.tree selfloop.dtl --from ranger --d 1
[ranger] donor 端点约定=donor_itself … 端点命中率=100.0% | 文件=1 块=1 转移=2 保留=2
         丢弃[标签不匹配=0, 支持度=0, 家族规模=0, 家族规模探测失败=0, 自环=1]
1.0 uninformative (  0.0 to a descendant, 0.0 to a leaf 0.0 to an ancestor 1.0 to itself 50.0%)
```

`--d 1` 与 `--d 0` 两次运行的 `to itself` 完全相同（都是 `1.0`），说明自环穿过了距离过滤。

实测二（同一份物种树，输入换成 **`--from ale`**：一份含 `T@62->65` 的 `.uml_rec`——默认
`trf` 口径下 donor = `parent(62)` = `65` = 受体，正是一次“供体 = 受体谱系”的自环。文件全文
即 `1 reconciled` / `FAM_SELFLOOP` /
`(((CYAP8:1.0,CYAP0:1.0)59:1.0,(CYAA5:1.0,CYAP2:1.0)61.T@61->59:1.0)65:1.0,(NOSP7:1.0,ANAVT:1.0)62.T@62->65:1.0)69:1.0;`）：

```
$ maxtic-next species.tree selfloop.uml_rec --from ale
[ale] 约束口径：source=trf（trf=转移事件 parent(donor)->receptor，ALE 官方 MaxTiC 集成口径），供体端点约定=parent_of_donor，产出边数=2，总权重=2.0
[ale] donor 端点约定=parent_of_donor 权重口径=support_over_declared_samples 分母=1 端点命中率=100.0% | 文件=1 块=1 转移=2 保留=2 丢弃[标签不匹配=0, 支持度=0, 家族规模=0, 家族规模探测失败=0, 自环=1]
selfloop.uml_rec
tree with  5 internal nodes
2.0 total weight of constraints from transfers
1.0 uninformative (  0.0 to a descendant, 0.0 to a leaf 0.0 to an ancestor 1.0 to itself 50.0%)
```

即 ALE 与其余四个适配器结论一致：`自环=1` 只是**计数**，那条约束照样保留并贡献了
`to itself` 的 `1.0`（默认 `--d 0` 下没有丢弃，因此也没有下面那种 NOTE）。
两阶段产物（`-o`）同样如实反映：自环行的第 5 列**留空**，而不是写 `0.0`。

```
$ maxtic-next species.tree selfloop.uml_rec --from ale -o unified.tsv
[ale] …（两行 stderr 诊断同上）…
[约束生成] ale 约束已写出至：unified.tsv
$ cat unified.tsv
#family,donor,receptor,weight,distance
selfloop,65,59,1.0,2.0
selfloop,65,65,1.0,
```

**文本输入路径一直是这个行为**：3 列空格 / 4 列逗号输入没有距离列，自环自然全部保留并计入
`to itself`。只有当你自己在文本输入里**显式写出** 0 距离时才会被删，此时程序如实报出：

```
NOTE: --min-transfer-distance=1.0 丢弃 1/2 条约束（保留 1 条）。
NOTE: 上述丢弃中有 1 条是**自环约束**（donor == receptor，其拓扑距离恒为 0）；
      它们不会进入原版口径的 "to itself" 统计。
```

即：`to itself` 是否非零，取决于输入有没有“距离列 + 距离 0”，而不是取决于事件是不是自环。

## 3.3 两阶段流程参数

| 参数 | 别名 | 说明 |
|------|------|------|
| `-o` | `--constraints-out` | **仅 `--from <tool>` / `--from-ale` 模式**：把转换得到的统一约束集写出到该路径（5 列 ALE 逗号格式，可被文本解析器读回），随后**提前返回不排序**，供 Snakemake/Nextflow 复用。与 `--dry-run` **同时给出时也会写出**该文件（`-o` 不再被静默忽略）。文本约束模式下给出 `-o` 会**直接报错**（输入已是统一格式，静默忽略会误导） |

## 3.4 剪裁与预检

| 参数 | 说明 |
|------|------|
| `--target-clade <标签>` | 目标类群根节点标签；默认仅保留两端点均在类群内的约束 |
| `--target-clade-ancestor-map` | 开启“外部祖先映射到根”的保守剪裁规则（需配合 `--target-clade`，默认关闭以避免虚假偏序） |
| `--dry-run` | 静态预检：树格式/二叉性/内部标签唯一/叶子名唯一/约束格式/权重合法性/端点合法性/`d` 列语义；`--from` 模式下**真正解析上游输出**并报告零约束与命中率。**任何一次适配器解析失败都会让预检判为失败**并以**退出码 1** 中止（，实例见 08 章 8.4b）。不运行排序 |

## 3.5 报告与输出

| 参数 | 取值 | 默认 | 说明 |
|------|------|------|------|
| `--no-html` | — | 关（即默认**生成** HTML） | 关闭交互式 HTML 报告 |
| `--output-style` | `short`/`legacy` | `short` | 输出文件命名风格；`legacy` 为原版长名。**两种风格内容逐字节一致** |
| `-p` | `--output-prefix` | 无 | 输出前缀；默认取**第一个**约束文件路径（见下） |
| `-f` | `--force` | 关 | 允许覆盖已存在的输出文件。默认**拒绝覆盖**并以**退出码 3** 中止 |

> **前缀与多文件**：产物前缀取自 `CONSTRAINTS` 列表的**第一个**文件，
> 其余文件参与计算但不出现在文件名里。多文件且未指定 `-p` 时会向 **stderr** 声明一次：
>
> ```
> NOTE: 输入了 2 个约束文件，三输出文件与 HTML 报告均以第一个文件为前缀：<file1>
> （其余 1 个文件已参与计算：['<file2>']）。如需为不同来源分别留存产物，
> 请对每次运行显式指定 -p/--output-prefix。
> ```
>
> 同时写入 `result.run_metadata["multi_file_prefix_note"]`。
>
> **写入语义**：产物先写同目录临时文件、`fsync` 后 `os.replace` 原子落盘，中途崩溃不会
> 留下“半截”文件；缺失的父目录自动创建（创建失败给出带路径的错误）。

## 3.6 稳健性 / MCMC（进阶）

| 参数 | 类型 / 取值域 | 默认 | 说明 |
|------|--------------|------|------|
| `--mcmc` | — | 关 | 在拓扑相容的线性扩展上做可逆 Metropolis–Hastings 采样。**注意：这是一条尚未通过收敛诊断验证的链，其样本不得当作后验样本使用** |
| `--mcmc-iters` | int > 0 | `1000` | MH 链步数 |
| `--mcmc-temperature` | `auto` 或 float > 0 | **`auto`** | Boltzmann 温度。`auto` = `max(总权重, 1) / 100`，有效值会打印并记入 `run_metadata["params"]`；也可给正数以固定温度运行 |
| `--mcmc-burn-in` | int ≥ 0 | `iters` 的 **50%** | burn-in 步数：该段只推进链、不进入统计 |
| `--mcmc-thin` | int > 0 | `max(1, iters // 100)` | 抽稀间隔，用于降低相邻样本自相关 |
| `--near-optimal-top-k` | int ≥ 1 | `50` | **近优解收集容量**：稳健性/敏感性摘要最多可以考察多少个**去重**近优排序（API `top_k`，见 3.6b）。仅在 `--ls > 0` 时有意义 |

`--mcmc` 会在 stdout 追加口径声明、有效参数、样本统计与**收敛诊断**（实测 1000 步、
`--seed 42`、蓝细菌样例）：

```
[MCMC] Metropolis–Hastings over linear extensions (preliminary; convergence diagnostics not validated)
[MCMC] 口径声明：这是一条**初步实现**的 Metropolis–Hastings 链，收敛诊断未经验证；链上样本彼此自相关（非独立观测），未通过下列诊断前不得用作不确定性度量
[MCMC] 温度 T=22.184（auto：max(总权重,1)/100）；iters=1000、burn-in=500、thin=10
[MCMC] 保留 100 个样本，其中唯一排序 51 个；能量 mean=835.7672, min=795.6300（唯一样本数远低于样本数说明链未混合，均值不可用作不确定性度量）
[MCMC] 收敛诊断（convergence diagnostics，按抽稀前的记录链计）：唯一样本占比=0.074、相邻重复率=0.907、最长同态连续段占比=0.056、有效样本量 ESS=7.6（0.008 倍样本量）、状态覆盖=n/a（状态空间过大，未精确枚举）；能量 sd=21.2451
```

> 温度默认不再是旧文档写的固定 `0.01`：在总权重 2218.4 的样例上，`auto` 给出 T=22.184，
> 与 0.01 相差 3 个数量级（旧默认值下链几乎冻结，只覆盖极少数状态）。
> `--mcmc-temperature 0.01` 仍可显式固定。

### 3.6b `--near-optimal-top-k`：稳健性摘要的支持集上限

`--near-optimal-top-k`（API `api.rank(top_k=...)`、`Ranker.run(top_k=...)`）**封顶**
“基于局部搜索访问解的稳健性/敏感性摘要”最多能考察多少个**互不相同**（去重）的近优排序：

- 默认 **50**（`config.DEFAULT_NEAR_OPTIMAL_TOP_K`），与历史行为一致；过去这个数字写死在
  收集器里，摘要因此永远不可能看过第 51 个序；
- 取值由 `ranking/ranker.py:check_near_optimal_top_k()` 校验：**必须是 ≥ 1 的整数**。
  CLI 越界（如 `--near-optimal-top-k 0`）在 argparse 层即**退出码 2**
  （`argument --near-optimal-top-k: 近优解收集容量 必须 > 0，收到 0`）；
  Python 侧传 `0` / `-1` / `2.5` / `'50'` / `True` / `inf` 一律 `ValueError`，
  整数值浮点 `100.0` 归一化为 `100`；
- **生效值**记入 `Result.run_metadata["near_optimal_top_k"]`，同时回显在
  `run_metadata["params"]` 字符串（`--near-optimal-top-k 12`）和 HTML 报告
  （`top_k：12（本次报告使用 12）`）里，摘要 dict 里另有 `top_k` / `top_k_reported`。

**为什么这个上限要紧**：`--ls > 0` 时程序给出的**不是**对输入数据重抽样的稳健度
（摘要 dict 的 `is_resampling_robustness` 恒为 `False`），而是**局部搜索访问过的解集合**
的去重视图。所以 K 直接决定了这些频率统计的**支持集**有多宽：K 以下的序根本不参与计数。
实测（`examples/minitree.tree` + `examples/Cyano_CUTConstraints.tsv`，`--ls 1`、
`--local-search-max-iters 60000`、`--t 60`、`--seed 42`）：

| `--near-optimal-top-k` | 摘要保留并统计的去重序数 | 局部搜索累计见过的不同序数 |
|------------------------|----------------------------|----------------------------|
| `3` | 3 | 9141 |
| `10` | 10 | 9139 |
| `50`（默认） | 50 | 8982 |
| `400` | 400 | 7940 |

也就是说：链走出的不同序远多于 K，摘要只看其中目标值最小的前 K 个。调大 K 可让摘要
考察更宽的近优邻域，代价是内存与 `summary()` 计算量按 **K·n²** 增长（成对次序频率要
在 K 个长度 n 的序上逐对累加）。

> ⚠️ **术语规则不变**：调大 K **不会**把它变成后验。本模块仍只称“稳健性/敏感性摘要”，
> **严禁**使用“后验（posterior）/ 置信（confidence）/ 置信区间”等表述（见 08 章 8.1）。
> 报告它时必须同时给出 `--near-optimal-top-k` 的生效值，否则“节点 61 排在第 3 位的频率”
> 这类数字失去了可解释的分母。

## 3.7 性能与续跑（进阶）

| 参数 | 类型 / 取值域 | 默认 | 说明 |
|------|--------------|------|------|
| `--incremental` / `--no-incremental` | — | **开** | 增量 scoring：局部搜索每次区间旋转仅算增量 O(b−a)，远优于全量 O(\|E\|)。默认开启（与全量逐字节等价、实测 3.1–8.2×）；`--no-incremental` 退回原版全量口径 |
| `--checkpoint <路径>` | str | 无 | 检查点文件；搜索中断后可续跑（仅 `--local-search>0` 生效） |
| `--checkpoint-interval` | float > 0 | `60.0` | 检查点自动保存间隔（秒） |
| `--local-search-max-iters` | int ≥ 0 | `0` | **局部搜索迭代次数上限**（0 = 只受 `--ls` 时长约束）。设为正数可让固定 `--seed` 的局部搜索路径**跨机器可复现** |

## 3.8 可复现性的准确边界

| 路径 | 固定 `--seed` 是否可复现 | 说明 |
|------|--------------------------|------|
| 默认（`--ls 0`） | **是**，逐字节可复现 | 三启发式与输出全部由 `seed` 决定的随机流驱动；约束行乱序、多文件顺序均不影响产物 |
| `--ls > 0` 且未设迭代上限 | **不保证跨机器复现** | 停止条件是**墙钟时长**，迭代次数依机器负载而变。程序会主动在 stdout 末尾打出一条 WARNING 说明这点 |
| `--ls > 0` + `--local-search-max-iters N` | **是**（同一 Python/平台下） | 迭代次数确定，随机流消耗确定；这是需要确定性搜索时的正确做法 |
| `--r 1` / `--r 2` | 由 `--seed` 决定 | 随机化会改变数据本身，因此“与默认路径一致”无从谈起，但同种子仍自洽可复现 |

> 换言之：手册与 README 过去笼统承诺的“固定 `--seed` 可复现”，
> 准确表述是**“默认路径（`--ls 0`）可复现；`--ls > 0` 需要
> `--local-search-max-iters` 才跨机器可复现”**。

## 3.9 常用命令速查

```bash
# 基本
maxtic-next species.tree constraints.tsv --seed 42

# 局部搜索 180 秒 + 关 HTML（要跨机器确定性就加 --local-search-max-iters）
maxtic-next species.tree constraints.tsv --ls 180 --local-search-max-iters 2000000 --no-html

# 按距离过滤（需 distance 列，否则该选项被忽略并报 ERROR）
maxtic-next species.tree ale_constraints.tsv --d 3

# 权重阈值过滤 + 随机树 p 值
maxtic-next species.tree constraints.tsv --ts 0.1 --rd 1000

# 上游工具（显式 / 自动）
maxtic-next species.tree FAM1.dtl --from ranger
maxtic-next species.tree out.recphyloxml --from-auto

# ALE：显式回到旧口径 rec（默认已是 trf）
maxtic-next species.tree fam.uml_rec --from ale --ale-source rec

# 目标类群剪裁
maxtic-next species.tree constraints.tsv --target-clade 61

# 预检
maxtic-next species.tree constraints.tsv --dry-run

# MCMC（自动温度 + 默认 burn-in/thin；要固定温度就 --mcmc-temperature 0.01）
maxtic-next species.tree constraints.tsv --mcmc --mcmc-iters 2000

# 让稳健性/敏感性摘要考察更多近优拓扑（默认 K=50）
maxtic-next species.tree constraints.tsv --ls 180 --near-optimal-top-k 400

# 压缩输入直接喂（无需先解压，见 05 章 5.7）
maxtic-next species.tree.gz constraints.tsv.gz --seed 42
maxtic-next species.tree examples/reconciliations.tgz --from ale

# 上游输出可疑时：预检会真正解析它，解析失败即退出码 1
maxtic-next species.tree out.recphyloxml --from eccetera --dry-run; echo "exit=$?"

# 长搜索 + 增量 + 续跑
maxtic-next species.tree constraints.tsv --ls 600 --incremental \
    --checkpoint run.ckpt --checkpoint-interval 30

# 重跑同一目录、明确允许覆盖上次产物
maxtic-next species.tree constraints.tsv --seed 42 --force
```

## 3.10 退出码

| 退出码 | 含义 |
|--------|------|
| `0` | 成功（含 `--dry-run` 预检通过） |
| `1` | `--dry-run` 预检存在 **error 级**问题（依据结构化 `run_metadata["dry_run_ok"]`，不是子串匹配）。`--from` 模式下“上游输出根本解析不了”也走这一条 |
| `2` | argparse 取值域错误（非法 `--r`、越界 `--ts`、负 `--ls`、未知 `--ale-source`、`--near-optimal-top-k < 1` 等），由 argparse 自身返回 |
| `3` | 输出文件已存在且未加 `--force`（覆盖保护） |
| `4` | 运行时异常：输入不满足算法前置条件（非二叉树、重复叶子名、内部标签不唯一、约束端点不在树中、权重非法/为负、适配器 0 命中等） |

> 容器中该退出码直接传播（`ENTRYPOINT ["maxtic-next"]`）。

> **图形前端的退出码**：桌面端 `maxtic-studio` 在独立仓库 `MaxTiC-Next-Studio`（见 README 的「图形前端在哪里」），
> 它复用同一张表的语义（缺 PySide6 时打印安装指引并以 **`3`** 退出）。本仓库不再注册该入口点。

> 上一章：[02 · 快速开始](02_quickstart.md) ｜ 下一章：[04 · Python API](04_python_api.md)
