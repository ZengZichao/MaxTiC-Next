# 05 · 输入输出格式

## 5.1 输入 1：物种树（Newick）

- **有根**的 Newick 树。
- **所有内部节点标签必须写在 bootstrap 字段**（`)标签:枝长`），这是 MaxTiC 的核心约定；
  约束里的 `donor`/`receptor` 就用这些标签引用内部节点。
- 叶子用叶名引用。
- **必须是二叉树**：每个内部节点恰有 2 个子节点。原版 `opt`/`mix` 强依赖二叉假设，
  本版把它上升为**前置条件**——多歧分支直接报错（退出码 4），而不是静默丢节点
  。`--dry-run` 也会以 error 级报出（见 08 章 8.4）。
- **内部节点标签必须唯一**，**叶子名也必须唯一**：重复叶名会让节点标识塌缩。
- 超度量（ultrametric）**非必需**，但若提供可用于与输出排序做 Kendall 相似度比较。
- 文件可以**折行**书写：解析前把整份文件的行拼接为一条序列（原版只看首行，这是已记录的
  偏离）。以 `utf-8-sig` 读取，带 BOM 的文件不会污染首列。
- 可以是 `.gz` / gzip 流 / **单成员** `.tar.gz`·`.tgz`·`.tar` 归档，无需先解压（，
  见 5.7）。

示例（内部标签 `42/41/59/61/...`；下面这种三行缩进写法可直接解析）：

```
(((((CYAP8:1,CYAP0:1)42:7,(CYAA5:7,UCYNA:7)41:1)59:2,
  (MICAN:9,(CYAP2:3,CYAP7:3)43:6)45:1)61:2,(SYNY3:11,SYNP2:11)56:1)65:1,
 (TRIEI:6,(NOSA0:5,(NOSP7:4,(NOSS1:2,ANAVT:2)40:2)46:1)62:1)67:7)69;
```

## 5.2 输入 2：约束文件（双格式）

一行一条约束，两种格式**自动识别**（含逗号即视为逗号格式）：

### 5.2.1 空格格式

```
donor receptor [weight] [distance]
```

- `weight` 可选，缺失默认 `1.0`；
- 第 4 列 `distance` 可选，为**系统发育距离**，供 `--min-transfer-distance` 过滤。
  **缺失距离列的约束一律保留**——即 `--d > 0` 对纯 3 列输入完全无操作，
  并且会打出一条 ERROR 级口径声明告诉你“该选项被忽略”（，见 03 章 3.1）。

```
61 62 15.04
61 67 15.48 3.0     # 第 4 列 3.0 为 distance，行内 # 之后是注释
```

### 5.2.2 逗号格式

```
gene_family,donor,receptor,[weight],[distance]
```

- **首列 `gene_family` 被丢弃**（仅用于追溯）；
- 4 列 = `family,donor,receptor,weight`；
- 5 列 ALE = `family,donor,receptor,weight,distance`（第 5 列为距离）。

```
fam001,147,149,0.09
fam002,197,187,0.12,2.5
```

### 5.2.3 ⚠️ 严禁的三列逗号误用

```
donor,receptor,weight     # 错误！会被拒绝并报错
```

因为逗号格式**首列一律丢弃**，三列逗号会把 `donor` 当 family 丢掉，导致 donor/receptor
错位、权重丢失。解析器检测到三列逗号会抛 `ValueError`（退出码 4）。若你本意是空格格式，
请改用空格分隔。

### 5.2.4 注释、特殊行与取值校验

- 以 `#` 开头的行被跳过；**行内 `#` 注释也被支持**（`61 67 15.48 # 说明`，
  `#` 前须为空白或位于行首，以免切断含 `#` 的类群名）；
- 含 `FRQ` 的行被跳过（兼容 ALE 汇总行）；
- **权重与距离必须为有限实数、权重必须 ≥ 0**：`nan` / `inf` / 负权重在**解析层**即被
  拒绝，错误信息含**文件名与行号**，非法值不会一路带进目标函数。
  数值列只接受十进制实数，不接受 `-` / `NA` / 空串。
- 同样透明支持 `.gz` / gzip 流与归档输入（，见 5.7）。

### 5.2.5 自环约束与 distance 列

`donor == receptor`（供体与受体同一谱系）的约束拓扑距离恒为 0，会被 `--min-transfer-distance`
以“距离太近”删掉、从而从 `to itself` 统计里消失。**五个适配器 RANGER-DTLx / ecceTERA / ARTra /
AleRax / ALE 都不会给出这样的 0**：自环一律写 `distance = None`（无距离信息），按“`--d` 保留无距离
列的约束”的既有口径通过过滤（；ALE 是最后一个对齐的，回归测试 `tests/test_ale_selfloop.py`，
见 03 章 3.2.3）。`--from ale -o` 写出的 5 列文件因此把自环行的第 5 列**留空**
（`selfloop,65,65,1.0,`），读回后端点与统计一致——解析器把**空的距离列**视同“无距离列”，
这是 5.2.4“数值列不接受空串”的唯一例外，且仅适用于距离列。
文本输入一直是这个行为——3 列空格 / 4 列逗号没有距离列，自环自然保留；只有你**自己**在第 4 列
（或逗号第 5 列）显式写 `0.0` 时它才会被删，此时 stdout 会打一条
`NOTE: 上述丢弃中有 N 条是**自环约束**…` 说明 `to itself` 少算了什么。

## 5.3 输出：文件产物

前缀默认为**第一个约束文件的路径**（CLI `-p/--output-prefix`，API `output_prefix`）。
输入多个约束文件时，其余文件参与计算但不出现在文件名里，且会在 **stderr** 显式声明一次
（，见 03 章 3.5）。命名风格由 `--output-style` 决定：

| 内容 | short（默认） | legacy（原版长名） |
|------|---------------|--------------------|
| 过滤后的加权信息性约束 | `<前缀>.mt.informative.tsv` | `<前缀>_MT_output_filtered_list_of_weighted_informative_constraints` |
| 与最优序冲突的约束 | `<前缀>.mt.conflicts.tsv` | `<前缀>_MT_output_list_of_constraints_conflicting_with_best_order` |
| 偏序 | `<前缀>.mt.partial_order.tsv` | `<前缀>_MT_output_partial_order` |
| 随机序分布（**仅 `--random-trees > 0`**） | `<前缀>.mt.random_dist.tsv` | `<前缀>_distribution_random` |

> - 前四行在两种风格下**内容逐字节一致**，只是文件名不同。
> - **覆盖保护**：目标文件已存在且不是本次运行写出的，CLI 以**退出码 3** 中止；
>   显式 `-f/--force` 才覆盖。
> - **原子写入**：先写同目录临时文件、`fsync` 后 `os.replace`；崩溃不留半截产物。
>   缺失的父目录会自动创建。

### 各文件内容

- **informative**：`donor,receptor weight` 每行一条，过滤掉无信息约束后的加权信息性约束。
  与谱系同向而被置为 `MAX_NUMBER`（1e10）的硬约束、以及 0 权重边不出现在其中。
- **conflicts**：`donor,receptor weight`，与最优排序**冲突**（被违反）的约束。
- **partial_order**：`a b weight color`，与最优排序**一致**的约束（构成偏序）；
  - `black`：与输入树排序方向一致；
  - `green`：与输入树排序方向相反（即输出相对输入发生了翻转）。
  - 判定“哨兵边”的阈值是 `MAX_NUMBER`（1e10）而不是原版写死的 `100000`，
    因此权重 ≥ 1e5 的真实信息性约束不再被静默删除。
- **random_dist**：每行 `value similarity`，共 `--random-trees` 行，是 p 值的经验零分布。

## 5.4 输出：stdout 摘要

逐行复刻原版 `MaxTiC.py` 的 print 顺序与措辞（完整实例见
[02 章 2.3](02_quickstart.md)）：内部节点数、约束总权重、uninformative 分项与占比、
trivially conflicting 占比、input/greedy/mixing 三种排序的 value 与百分比、best order 标识、
排序树 Newick、与输入树的 Kendall 相似度。

随选项追加的行：

| 触发 | 追加内容 |
|------|----------|
| `--ls > 0` | `attempting a local search …` 提示行；摘要中 `best order is the …` 之后补 `after local search X rejected` 与 `best found solution Y (%)` 两行（，使 stdout 不与交付树矛盾） |
| `--rd N > 0` | 4 行：随机序 value 区间、value 的 p 值、随机序 similarity 区间、similarity 的 p 值；p 值用 **(k+1)/(n+1)** 校正且检验**最终交付序** |
| `--mcmc` | 5 行 `[MCMC]`：状态声明、口径声明、有效温度与 burn-in/thin、样本统计、收敛诊断 |
| 总是（按需） | 摘要**末尾**追加 `NOTE:` / `WARNING:` / `ERROR:` 运行期口径声明与告警（谱系硬约束计入、空约束集、`--d` 被忽略、`--ls` 不可跨机器复现、自环被距离过滤删掉等） |
| 被**取消**（API/GUI `stop_check`） | `NOTE: 本次运行按请求**取消**：局部搜索在第 N 次迭代边界停止（计划的搜索时长为 T 秒，未跑满）……`，交付的是截至取消点的历史最优序（见 08 章 8.9） |

百分比口径（本版与原版已裁定的偏离）：

- `uninformative` 的百分比 = `uninformative ÷ total_weight`。原版分母
  `total_transfers + uninformative` 把无信息权重量**双重计入**，同一实例上报 `12%`，
  本版报 `13.9%`。
- `--threshold-constraints` 生效时，总权重在**删边之后**统计（原版在删边之前，
  分母里仍含被删权重）。

## 5.5 输出：交互式 HTML 报告（默认）

- 纯静态单文件 `<前缀>.html`，无需后端，可离线打开。
- 内容：运行元数据（种子/参数/版本/适配器诊断）、约束权重分布、冲突比例、
  节点位置稳健性摘要（局部搜索开启时）、排序树可视化。
- 关闭：`--no-html`。
- 未安装 `jinja2` / `plotly`（`pip install "MaxTiC-Next[report]"`）时**不报错**，
  降级为基础模板：无交互图、且降级模板**不做 HTML 转义**。要完整报告请装 extras。

## 5.6 输出：稳健性/敏感性摘要（局部搜索开启时）

在 `Result.sensitivity_summary`（Python）与 HTML 报告中呈现：top-K 近优序、
每节点位置分布、成对次序频率。**仅称“稳健性/敏感性摘要”**，不使用“置信区间/后验概率”
（那是 MCMC 的专属语义，而 MCMC 本身也尚未通过收敛诊断，见 [08 章](08_advanced_features.md)）。

它是**局部搜索访问集合的去重视图**（`is_resampling_robustness` 恒为 `False`），
统计的支持集由 `--near-optimal-top-k`（API `top_k`）封顶：默认 **50**，即摘要最多看目标值
最小的 50 个去重近优序，而局部搜索实际走过的不同序数往往远多于此
（`n_unique_orders_total`）。生效值记入 `run_metadata["near_optimal_top_k"]`，
摘要 dict 里也有 `top_k` / `top_k_reported`；报告这些频率时必须同时给出 K，
详见 03 章 3.6b 与 08 章 8.1。

## 5.7 压缩输入：`.gz` / gzip 流与 tar 归档（透明解压）

上游工具的示例与结果**大多以压缩形式分发**（MaxTiC 的 `examples/reconciliations.tgz`、
ALE/Ranger/AleRax 输出常被打包成 `.tar.gz`，单文件示例常见 `xxx.gz`）。
`io/compression.py` 是全项目**唯一**的文本输入读取底层，因此**不需要先手工解压**。
它同时按**内容魔数**和**扩展名**判定，避免“只看扩展名”和“只看魔数”各自的失效场景：

| 输入 | 判定依据 | 行为 |
|------|----------|------|
| `x.gz` / 任何首两字节为 `1f 8b` 的文件 | gzip 魔数（RFC 1952） | 透明解压为文本。**无扩展名也能读** |
| 扩展名是 `.gz` 但内容**不是** gzip | 魔数不匹配 | 按**纯文本**读取（误命名的未压缩文件仍然可用，不再抛 `BadGzipFile`） |
| 恰好**一个**文件成员的 `.tar.gz` / `.tgz` / `.tar`（含被改名的） | `ustar` 魔数（偏移 257）+ `tarfile` | 透明取出该成员内容，等同普通文本文件 |
| **多个**文件成员的 tar 归档 | 同上 | **不**自动拼接：抛 `CompressedArchiveError`（`ValueError` 子类，CLI 退出码 4），消息里给出成员数、前若干成员名与**可直接执行**的 `tar -xzf 归档 -C 目录` / `tar -xzOf 归档 成员 > 文件` 命令 |

### 5.7.1 哪些入口真的吃得到

- **物种树（5.1）**：`.gz` / gzip 流 / **单成员**归档都可用（走 `io.parsing.read_newick_file`）。
  多成员归档在这里**就是报错**——一棵树只能有一个成员。
- **约束文件（5.2）**：同上（`io.parsing.read_constraints_file`、便捷封装
  `constraints.parsers.parse_constraints_file`）。
- **`--dry-run` 预检**：同一读取层，故压缩输入在预检里同样透明；多成员归档给出的是
  **error 级报告条目**而非裸回溯。
- **上游适配器**：ALE / RANGER-DTLx / ecceTERA / ARTra / AleRax 的逐文件解析都走本层，
  所以 `FAM1.dtl.gz`、`fam.uml_rec.gz` 可直接作为 `--from <tool>` 的输入；
  `--from-auto` 的格式嗅探也读的是**解压后**的文本（旧实现会把压缩字节当乱码判为“未识别”）。
  家族名取自**剥掉压缩后缀之后**的词干（`gene_1.modif.ale.uml_rec.gz` 与
  `gene_1.modif.ale.uml_rec` 得到**同一个** `metadata["family"]`），故 `.gz` 不改变按家族
  统计的口径。
- **归档输入（`CONSTRAINTS` 位置参数 / `api.rank` 的约束列表）**：`api.rank` 会先把
  `.tgz` / `.tar.gz` / `.tar` **安全解压**到进程临时目录，再把**成员当作输入**，
  所以多成员归档（例如 ALE 官方 `examples/reconciliations.tgz`，内含上千个 `*.uml_rec`）
  在约束位置参数上**可以直接喂**；此时 stderr 打一条 `NOTE: 输入 … 是 tar 归档，已自动展开为
  N 项输入（解包于 …）`，映射关系记入 `run_metadata["archives_expanded"]`
  （成员活在临时目录里，必须可追溯）。AleRax 的输出树归档按目录约定解出**根目录**。
  解压有安全边界：只解常规文件，绝对路径 / 含 `..` / 链接与设备成员一律跳过。

> **不要高估它**：拼接多成员归档为**一个**文本输入永远不做——归档内每个成员是一份独立的
> 上游输出，拼在一起会被当成一个基因家族，从而按 `--ale-min-family-size` /
> `--ale-min-support` 之类的**按家族**口径整体误算。多成员归档只有在“每个成员各是一份
> 输入”的位置（约束/上游输出列表）才会被自动展开。

### 5.7.2 实测

```bash
# .gz 约束输入 —— 产物与未压缩输入逐字节相同
maxtic-next minitree.tree constraints.tsv.gz --seed 42          # exit 0
diff constraints.mt.informative.tsv constraints_gz.mt.informative.tsv   # 无差异

# 单成员 .tar.gz（约束）与 .tgz（物种树）同样可用
maxtic-next minitree.tree single.tar.gz --seed 42               # exit 0
maxtic-next single.tgz constraints.tsv --seed 42                # exit 0

# 官方多成员 .tgz 直接作为上游输入
maxtic-next species.tree reconciliations.tgz --from ale         # 自动展开为 N 个 .uml_rec

# 多成员归档被当作"一个"输入读（此处是物种树位置）→ 退出码 4 + tar 命令
maxtic-next multi.tar.gz constraints.tsv --seed 42; echo "exit=$?"
错误：输入文件是一个含 2 个文件成员的 tar 归档（gzip）：/tmp/multi.tar.gz
  成员：Cyano_CUTConstraints.tsv, minitree.tree
  说明：MaxTiC-Next 透明解压**单文件** gzip（.gz / 无扩展名的 gzip 流）与**恰好一个成员**的
  tar 归档；多成员归档不会被自动拼接成一个输入，……
  请任选其一，把成员变成真正的输入：
    1) 解包整个归档，再把解出的文件（或目录，AleRax 支持目录）作为输入：
       tar -xzf /tmp/multi.tar.gz -C ./unpacked && ls ./unpacked
       maxtic-next species_tree ./unpacked/* [--from auto]
    2) 只取其中一个成员：
       tar -xzOf /tmp/multi.tar.gz Cyano_CUTConstraints.tsv > Cyano_CUTConstraints.tsv
  （Windows 10+ 自带 tar；--dry-run 会原样复述本提示）
exit=4
```

本模块只依赖标准库（`gzip` / `tarfile`），不引入任何第三方解压依赖。

> 上一章：[04 · Python API](04_python_api.md) ｜ 下一章：[06 · 上游软件配套方案](06_upstream_integration.md)
