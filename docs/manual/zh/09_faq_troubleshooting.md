# 09 · 常见问题与排错

## 9.1 FAQ

**Q1：MaxTiC-Next 和原版 MaxTiC 有什么关系？结果一样吗？**
A：MaxTiC-Next 是原版（Python 2）的 Python 3 忠实重写，核心算法逐函数移植；可达性矩阵加速
与原版 `path()` 贪婪经 400/400 例逐例核对一致。但**不是一句“逐字节等价”就能概括**：
本版有意修正了原版的若干口径缺陷（uninformative 百分比分母、偏序哨兵阈值、阈值过滤后的
总权重、局部搜索后的 value 回写、随机树 p 值校正、ALE 约束源默认值等），并加了几道原版没有
的前置条件校验。完整清单见 README 的“与原版的有意偏离”一节。
`--output-style legacy` 只改变**文件名**（长名与原版一致），两风格内容逐字节相同。

**Q2：约束里的 donor/receptor 用什么标识？**
A：内部节点用**物种树 bootstrap 字段里的标签**（如 `61`），叶子用叶名。约束端点必须在物种树标签集合内。

**Q3：`--seed` 能保证跨环境完全一致吗？**
A：分三种情况，别再说成一句“固定 seed 可复现”：
- **默认路径（`--ls 0`）**：可复现，逐字节一致。约束行乱序、多文件顺序都不影响产物。
- **`--ls > 0`**：停止条件是**墙钟时长**，迭代次数依机器负载而变 → 同种子在不同机器上
  **可能不同**（程序会在 stdout 末尾主动打这条 WARNING）。
- **`--ls > 0` + `--local-search-max-iters N`**：迭代次数确定 → 恢复确定性。

另外，贪婪启发式在**同权重并列**时依赖 `edge` 字典的插入顺序，Newick 解析顺序也参与其中，
所以“可复现”始终限定为“固定种子 + 同解析器 + 同 Python 版本”。建议锁版本、容器化。

**Q4：为什么我的 `--d` 没有生效？**
A：`--d` **只按 phylogenetic distance 列过滤**。空格格式需第 4 列 distance，逗号需 5 列 ALE 格式。
若约束**没有**距离列，则 `--d>0` **一条也不会删**（无距离列的约束一律保留，与原版一致）——
这是“被忽略”，不是“全部丢弃”。现在这一事实以 **error 级**告诉你：

```
ERROR: --min-transfer-distance=5.0 被忽略（123/123 条约束不含 phylogenetic distance 列，……）
```

排序本身仍退出 0；`--dry-run` 则报 error 并退出 1。要真按距离过滤请补距离列，不需要就去掉该选项。

**Q5：局部搜索的输出能当置信区间/后验概率吗？**
A：**不能**。局部搜索是 Metropolis 邻域搜索，只提供“稳健性/敏感性摘要”。
`--mcmc` 是唯一在术语上允许谈分布的模块，**但它本身是尚未通过收敛诊断的初步实现**——
实测 ESS 远小于样本量、相邻样本重复率高，链未混合。在其诊断通过之前，两边的数字都
**不应**作为不确定性度量发表。

**Q6：上游工具的小家族样例被过滤没了怎么办？**
A：演示小样本加 `--ale-min-family-size 0`；真实数据用默认 `5` 或按需调。

**Q7：ecceTERA 转换后约束很少 / 为空？**
A：ecceTERA 对内部物种节点用其自身数字 ID（叶子用名字），未必等于你物种树的 bootstrap 标签，
不一致端点会被丢弃并**计数上报**。请看 stderr 的诊断行（`端点命中率=…`、
`丢弃[标签不匹配=N,…]`），或取 `result.run_metadata["adapter_diagnostics"]`。
命中率低于 `--min-endpoint-hit-rate`（默认 0.5）会告警，**命中率为 0 直接报错**。
解决办法是对齐物种树内部标签与 ecceTERA 命名（或反向，把 bootstrap 换成数字 ID）。

**Q8：HTML 报告没生成？**
A：完整报告需 `jinja2 + plotly`（`pip install "MaxTiC-Next[report]"`）。缺依赖时**不会报错**，
而是降级为基础模板（无交互图，且降级模板不做 HTML 转义）；生成异常只在 stdout 打一行
`[警告] HTML 报告生成失败，已跳过：…`。排序与文件输出不受影响。也可 `--no-html` 主动关闭。

**Q9：可以只生成约束、不排序吗？**
A：可以。`--from <tool> -o constraints.tsv` 生成统一约束后提前返回（两阶段流程的 Stage 1）。
`-o` **只在适配器模式有效**：文本约束模式给 `-o` 会直接报错而不是静默忽略；
`--dry-run -o …` 现在也照常写出该文件。

**Q10：物种树是多歧树（polytomy）能跑吗？**
A：**不能**。原版 `opt`/`mix` 强假设二叉，本版把该假设上升为**前置条件**：多歧分支直接
`ValueError`（退出码 4），不再静默丢节点导致目标函数失真。请先解消多歧或重新定根
（未定根 Newick 的根天然三出）。`--dry-run` 同样报 error。

**Q11：为什么第二次运行同样的命令报了“输出文件已存在”？**
A：这是**覆盖保护**：默认拒绝覆盖上次产物，退出码 `3`，以免静默破坏已有结果。
确认要覆盖加 `-f/--force`，或用 `-p` 换前缀。同时产物现在是**原子写入**（临时文件 +
`os.replace`），崩溃不会留下半截文件。

**Q12：多个约束文件时产物只体现第一个文件名？**
A：是。前缀取自 `CONSTRAINTS` 列表的**第一个**文件，其余文件参与计算但不进文件名。
现在这种情况会向 **stderr** 打一条 `NOTE: 输入了 N 个约束文件…` 显式声明，并写入
`run_metadata["multi_file_prefix_note"]`。要按来源分别留存产物，请对每次运行显式 `-p`。

**Q13：`--ale-source` 到底该用 `rec` 还是 `trf`？**
A：默认是 **`trf`**（= ALE 官方 MaxTiC 集成规定的 `constraints_from_transfers`）。
两者产出的约束集**规模与权重量级都不同**，切换等于换掉整个输入。只有在复现以 `rec`
为默认的旧结果时才显式 `--ale-source rec`。生效口径会打到 stderr。

**Q14：上游示例是 `.tgz`，怎么直接喂给 MaxTiC-Next？**
A：**不用先解压**。`io/compression.py` 按 gzip 魔数（`1f 8b`）和 `ustar` 魔数判定输入，
所以 `.gz` / 无扩展名的 gzip 流 / **恰好一个成员**的 `.tar.gz`·`.tgz`·`.tar` 在**任何**入口
（物种树、文本约束、`--dry-run`、五个适配器、`--from-auto` 的格式嗅探）都是透明解压的；
家族名取自剥掉压缩后缀后的词干，所以 `.gz` 不改变按家族统计的口径。

上游工具的**官方示例包**现在也能直接喂：例如 ALE 官方
`examples/reconciliations.tgz`（里面是上千个 `*.uml_rec`）可以**整包**作为约束/上游输入 ——
`api.rank` / CLI 会先把归档安全解到进程临时目录，再把**每个成员**当作一份输入：

```bash
maxtic-next species.tree examples/reconciliations.tgz --from ale
# stderr: NOTE: 输入 … 是 tar 归档，已自动展开为 N 项输入（解包于 …）
# 映射记入 run_metadata["archives_expanded"]；产物前缀取第一项，必要时显式 -p
```

那“**含多个成员的归档被当作一个输入读**”时的报错是什么意思？——那是**故意不做拼接**：
归档里每个成员是一份独立的上游输出（各自的基因家族/样本），拼在一起会被当成**一个**家族，
从而按 `--ale-min-family-size` / `--ale-min-support` 这类**按家族**的口径整体误算。
所以这种输入（典型是把 `.tgz` 放在**物种树**位置，或绕过 CLI 直接调
`read_constraints_file`）抛的是 `CompressedArchiveError`（CLI 退出码 4、`--dry-run` 记为
error 级条目），消息里带成员清单与**可直接执行**的解包命令：
`tar -xzf 归档 -C ./unpacked`（然后喂 `./unpacked/*`），或
`tar -xzOf 归档 成员 > 文件`（只取一个成员）。Windows 10+ 自带 `tar`。详见 05 章 5.7。

**Q15：`--near-optimal-top-k` 是干什么的？要调吗？**
A：它是**稳健性/敏感性摘要的支持集上限**——局部搜索期间为摘要保留多少个**去重**近优排序，
默认 `50`（API `top_k`，`Ranker.run(top_k=...)`）。摘要不是重抽样，而是“访问过的解集合”的
去重视图，所以 K 之外的序压根不参与计数；实测链往往走过近万个不同序，摘要却只能看 50 个。
需要更宽的近优邻域就调大（内存/`summary()` 随 K·n² 增长）。生效值在
`run_metadata["near_optimal_top_k"]` 与 HTML 报告里，**报告数字时必须一并给出 K**。
取值必须是 ≥ 1 的整数：`--near-optimal-top-k 0` 退出码 2，Python 侧 `ValueError`。
调大 K **不会**让它变成后验（详见 03 章 3.6b、08 章 8.1b）。

**Q16：跑了一半想停，能保住已算出的结果吗？**
A：能。`api.rank(..., stop_check=...)` / `Ranker.run(..., stop_check=...)` 接受一个
`()`→`bool` 谓词，它在局部搜索主循环与 MCMC 链步之间被检查，于是取消发生在**迭代边界**：
程序**返回**截至当下的历史最优序、照常写产物，不抛异常也不毁掉解；是否真的打断看
`run_metadata["cancelled"]`（局部搜索 `stats_out["cancelled"]` 或 `mcmc_stopped_early`）。
Studio 的「取消 ■」按钮已经绑在这条路径上（Studio 仓库 `maxtic_studio/cancellation.py` 的 `StopToken`）——
过去那个按钮只是“等跑完再丢弃结果”。**唯一例外**：一条一个样本都没产出的 MCMC 链被取消时
会 `RuntimeError`，因为没有样本就没有可交付的量。详见 08 章 8.9。

**Q17：适配器诊断行里的 `自环=N` 是不是“丢了 N 条”？**
A：**不是**，那只是**计数**。`donor == receptor`（供体与受体同一谱系）的事件拓扑距离恒为 0，
若适配器真把这个 0 写进距离列，`--min-transfer-distance` 就会把它们删掉、让它们从原版口径的
`to itself` 统计里静默消失。现在**五个适配器 RANGER-DTLx / ecceTERA / ARTra /
AleRax / ALE** 对自环一律写 `distance = None`（= 无距离信息），按“`--d` 保留无距离列的约束”的
既有口径通过过滤，所以 `to itself` 不再被 `--d` 吃掉——这与**文本输入路径一直的行为一致**
（3 列空格/4 列逗号本就没有距离列）。ALE 是最后一个对齐的：`--from ale` 的自环现在同样保留，
实测 `to itself` 非零（`1.0 uninformative ( … 1.0 to itself 50.0%)`，见 03 章 3.2.3 的实测二，
回归测试 `tests/test_ale_selfloop.py`）。只有你在文本输入里**自己**写了距离 0（或人为把
`None` 改成 `0.0`）时它才会被删，此时 stdout 会打 `NOTE: --min-transfer-distance=… 丢弃 …` 与
`NOTE: 上述丢弃中有 N 条是**自环约束**…` 如实说明。见 03 章 3.2.3、05 章 5.2.5。

## 9.2 典型错误信息（逐字取自实现）

| 错误 / 现象 | 退出码 | 原因 | 解决 |
|-------------|--------|------|------|
| `错误：约束文件 <f> 第 N 行的逗号格式列数不足：'61,62,15.04'` | 4 | 三列逗号 `donor,receptor,weight` 误用 | 改空格格式，或补 family 首列成 4 列 |
| `错误：约束端点不在物种树中：donor='999', receptor='888'（约束文件/来源：<f>）` | 4 | donor/receptor 标签与物种树不符 | 核对标签；ecceTERA 需对齐内部数字 ID |
| `错误：物种树不是二叉树：内部节点 '1' 有 3 个子节点。` | 4 | 多歧分支 | 二叉解析或重新定根 |
| `错误：物种树存在重复的叶子名（共 1 个名字重复）：'A' 出现在父节点 ['9', '8']` | 4 | 叶名不唯一，节点标识塌缩 | 检查是否误把基因树/带编号拷贝当物种树 |
| `错误：物种树内部节点标签不唯一（重复标签：[…]）` | 4 | bootstrap 标签重复 | 修正物种树使内部标签唯一 |
| `约束文件 <f> 第 N 行的约束权重为负：'-3.0'` | 4 | 负权重使目标函数与百分比失去意义 | 校正权重（权重须为有限实数且 ≥ 0） |
| `约束文件 <f> 第 N 行的约束权重不是有限实数：'nan'` | 4 | `nan`/`inf` 会污染统计 | 用真实数值；缺失行请删除 |
| `输出文件已存在：<path>` + “请改用 -p… 或加 --force” | 3 | 覆盖保护（默认拒绝） | `-f/--force` 或换 `-p` 前缀 |
| `argument --threshold-constraints/--ts: 约束权重阈值比例 是比例，取值域 [0.0, 1.0]，收到 5.0` | 2 | `--ts` 越界（必须是 [0,1] 比例） | 用 0–1 之间的比例 |
| `argument --random-type/--r: invalid choice: '3' (choose from '0', '1', '2')` | 2 | 非法随机化类型 | 用 0/1/2 |
| `ERROR: --min-transfer-distance=… 被忽略（N/N 条约束不含 phylogenetic distance 列…）` | 0（`--dry-run` 为 1） | 无距离列却设了 `--d>0` | 去掉 `--d` 或提供 4/5 列输入 |
| `[tool] 错误：…` / `LabelMismatchError`（适配器端点 0 命中） | 4 | 上游命名与物种树不一致 | 对齐命名，或降低 `--min-endpoint-hit-rate`（0 命中仍报错） |
| `[error] 上游输出（格式 ranger）未产出任何约束：…` | 1（dry-run） | 上游格式误标 / 家族阈值过高 / 命名不一致 | 先看 stderr 诊断行；`--dry-run` 定位 |
| `[error] 上游输出（格式 eccetera）解析失败：UpstreamParseError: …recPhyloXML 解析失败（文件损坏 / 被截断 / 非 XML）…` | 1（`--dry-run`） | 上游输出本身读不出（截断 / 非 XML / 编码错） | 预检**必定**判失败并退出 1；重新导出或解压上游输出 |
| `错误：输入文件是一个含 N 个文件成员的 tar 归档（gzip）：<path>` + `tar -xzf` / `tar -xzOf` 命令 | 4（`--dry-run` 记 error 级条目） | 多成员归档被当作**一个**输入读（典型：放在物种树位置） | 按提示解包后喂 `./unpacked/*`；或 `tar -xzOf 归档 成员 > 文件` 只取一个。作为约束/上游输入时归档会自动展开（05 章 5.7） |
| `argument --near-optimal-top-k: 近优解收集容量 必须 > 0，收到 0` | 2 | 摘要支持集上限越界（须为 ≥ 1 的整数） | 用 ≥ 1 的整数；Python 侧同判据为 `ValueError` |
| `MaxTiC-Next Studio 需要 GUI 依赖 PySide6（以及 matplotlib）。请在 Studio 仓库根目录执行：pip install -e .` | 3 | 未装 GUI 依赖就运行 `maxtic-studio` | 在 Studio 仓库执行 `pip install -e .`；或改用 `maxtic-next`（CLI 不受影响） |
| `WARNING: 输入约束集为空（0 条时间约束）。…` | 0 | 数据或阈值导致无约束 | **别把 0.0 的 value 读成“完美一致”** |
| `WARNING: --local-search > 0 以墙钟时长为停止条件…` | 0 | 不可跨机器复现 | 加 `--local-search-max-iters` |
| `--target-clade 指定的节点 … 不在物种树中` | 4（dry-run 报 error） | 类群标签写错 | 用物种树内部 bootstrap 标签 |
| `maxtic-next: command not found` | — | 未 pip 安装 | `pip install -e .`，或用 `PYTHONPATH=src python3 -m maxtic_next` |
| HTML 报告“跑通了但很朴素” | 0 | 未装 jinja2/plotly → 静默降级 | `pip install "MaxTiC-Next[report]"` |

## 9.3 复现性核对清单

- [ ] 固定 `--seed`（默认 42）
- [ ] 走默认路径 `--ls 0`；若必须 `--ls > 0`，同时固定 `--local-search-max-iters`
- [ ] 锁定 Python 版本（推荐 3.11）与依赖版本
- [ ] 固定 `--ale-source`（默认 `trf`；切 `rec` 即换掉整个约束集）
- [ ] 使用同一 Newick 解析路径（本包自带 Tree）
- [ ] 记录完整命令与参数（HTML 报告与 `run_metadata["params"]` 均含运行元数据）
- [ ] 避免 `--random-type 1/2`（会改变数据本身）
- [ ] 容器化并**固定镜像标签**（`maxtic-next:0.1.1`，不要 `:latest`）
- [ ] 报告 p 值时写明 `--random-trees` 的 N 与“(k+1)/(n+1) 校正”
- [ ] 报告稳健性/敏感性摘要时写明 `--near-optimal-top-k` 的生效值（默认 50）——
      它是这些频率的支持集上限，见 03 章 3.6b
- [ ] 若上游输入是压缩包：`.gz` 与归档可直接喂，但**产物前缀**与临时解包位置会变
      （见 `run_metadata["archives_expanded"]`）；要按成员分别留存产物请显式 `-p`

## 9.4 性能建议

| 目标 | 手段 |
|------|------|
| 加速局部搜索 | 增量打分已默认开启（`--no-incremental` 可退回全量口径，两者结果一致） |
| 防长搜索中断丢失 | `--checkpoint` + `--checkpoint-interval` |
| 加速上游批量解析 | `--ale-parallel process` + `--ale-cache-dir`（断点续传，键含拓扑指纹） |
| 避免重复解析 | 两阶段：固化 `constraints.tsv` 后只重跑排序 |
| 控制随机检验成本 | `--random-trees` 按需取值（每次都要算一遍 value + similarity） |

## 9.5 获取帮助

- 手册总览：[../README.md](../README.md)
- 适配框架：[../../adapters/README.md](../../adapters/README.md)
- 原始工具作者：Eric Tannier（eric.tannier@inria.fr）；引用 MaxTiC, Biorxiv doi.org/10.1101/127548

> 上一章：[08 · 进阶功能](08_advanced_features.md) ｜ 返回：[手册总览](../README.md)
