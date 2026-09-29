# 08 · 进阶功能

## 8.1 稳健性 / 敏感性摘要（局部搜索）

打开局部搜索（`--local-search T > 0`）后，MaxTiC-Next 会在搜索过程中收集**近优解集合**，
并汇总：

- **节点位置分布**（`node_position_distribution`）：每个内部节点在近优序中出现在各排位的次数；
- **成对次序频率**（`pairwise_order_frequency`）：有序对 (a, b) 中 a 排在 b 之前的次数。

```bash
maxtic-next species.tree constraints.tsv --ls 30
```

```python
r = rank("species.tree", "constraints.tsv", local_search=30, print_summary=False)
summ = r.sensitivity_summary
print(summ["title"])                 # 基于局部搜索访问解的稳健性/敏感性摘要
print(summ["n_orders_collected"])    # 收集到的去重近优序数
print(summ["node_position_distribution"]["61"])   # 节点 61 的位置分布
```

> **术语铁律**：局部搜索是 Metropolis 邻域搜索，其访问解**不是**统计后验。因此本模块
> **只**称“稳健性/敏感性摘要”，**严禁**使用“置信区间 / 后验概率 / 后验分布”。
> 这条铁律**不因**下面调大 K 而有丝毫松动。

### 8.1b `--near-optimal-top-k`：摘要支持集的硬上限

摘要统计的不是“重抽样”，而是**局部搜索访问过的解集合**的去重视图
（`is_resampling_robustness` 恒为 `False`）。收集器按目标值升序只保留 **K 个**去重近优序，
因此 K 就是这些频率的**支持集上限**——K 以外的序根本不参与计数。它过去写死为 50，
现在由 CLI `--near-optimal-top-k` / API `api.rank(top_k=...)` / `Ranker.run(top_k=...)`
三层覆盖，取值由 `ranking/ranker.py:check_near_optimal_top_k` 校验（**≥ 1 的整数**；
CLI 越界退出码 2，Python 侧 `ValueError`），生效值记入
`Result.run_metadata["near_optimal_top_k"]`、回显在 `run_metadata["params"]` 与 HTML 报告，
摘要 dict 里另有 `top_k` / `top_k_reported`。

实测（`examples/minitree.tree` + `examples/Cyano_CUTConstraints.tsv`，`--ls 1`、
`--local-search-max-iters 60000`、`--t 60`、`--seed 42`）：

| `--near-optimal-top-k` | 摘要保留并统计的去重序数 | 局部搜索累计见过的不同序数 |
|------------------------|----------------------------|----------------------------|
| `3` | 3 | 9141 |
| `10` | 10 | 9139 |
| `50`（默认） | 50 | 8982 |
| `400` | 400 | 7940 |

```bash
maxtic-next species.tree constraints.tsv --ls 180 --near-optimal-top-k 400
```

```python
r = rank("species.tree", "constraints.tsv", local_search=180, top_k=400)
r.run_metadata["near_optimal_top_k"]        # 400 —— stdout 之外的机器可读入口
r.sensitivity_summary["n_orders_retained"]  # 参与统计的去重序数
```

- **调大 K** 让摘要考察更宽的近优邻域；代价是内存与 `summary()` 计算量按 **K·n²** 增长
  （成对次序频率要在 K 个长度 n 的序上逐对累加）。
- 调大 K **不**改变“这不是后验”这一事实，也**不**改善链的混合；它只改变统计的分母来源。
- 报告这些频率时必须同时给出 K 的生效值，否则“节点 61 出现在第 3 位的频率”没有可解释的分母。

> **可复现性**：`--ls > 0` 的停止条件是**墙钟时长**，迭代次数依机器负载而变，因此固定
> `--seed` **不保证跨机器复现**（程序自己也会打出这条 WARNING）。需要确定性请同时设置
> `--local-search-max-iters N`（见 8.5b）。默认路径 `--ls 0` 始终逐字节可复现。

## 8.2 MCMC 采样器（可选，初步实现）

`--mcmc` 开启在**与物种树拓扑相容的线性扩展**状态空间上的可逆 Metropolis–Hastings 采样器
（对称提议 + `min(1, exp((E_old−E_new)/T))`，细致平衡成立）。

```bash
maxtic-next species.tree constraints.tsv --mcmc --mcmc-iters 2000
```

| 选项 | 默认 | 说明 |
|------|------|------|
| `--mcmc-iters` | `1000` | 链步数（必须 > 0） |
| `--mcmc-temperature` | **`auto`** | Boltzmann 温度。`auto` = `max(总权重, 1) / 100`；也可给正数固定 |
| `--mcmc-burn-in` | `iters` 的 **50%** | burn-in：只推进链、不记录样本 |
| `--mcmc-thin` | `max(1, iters // 100)` | 抽稀间隔，降低相邻样本自相关 |

- **温度不再是旧文档写的固定 `0.01`**。旧默认值在真实总权重（数百到数千量级）下会让链
  几乎冻结、只覆盖极少数状态；现在默认 `auto` 按能量尺度自适应，并在 stdout 打印实际生效值：

  ```
  [MCMC] 温度 T=22.184（auto：max(总权重,1)/100）；iters=1000、burn-in=500、thin=10
  ```

- 结果在 `Result.mcmc_samples`（线性扩展样本列表），并打印样本统计与**收敛诊断**：

  ```
  [MCMC] 保留 100 个样本，其中唯一排序 51 个；能量 mean=835.7672, min=795.6300（唯一样本数远低于样本数说明链未混合，均值不可用作不确定性度量）
  [MCMC] 收敛诊断（convergence diagnostics，按抽稀前的记录链计）：唯一样本占比=0.074、相邻重复率=0.907、最长同态连续段占比=0.056、有效样本量 ESS=7.6（0.008 倍样本量）、状态覆盖=n/a（状态空间过大，未精确枚举）；能量 sd=21.2451
  ```

> ⚠️ **诚实的能力边界**：这是一条**尚未通过收敛诊断验证**的初步实现。程序自己
> 打出的第一行就是 `preliminary; convergence diagnostics not validated`，且实测显示链
> **未混合**（ESS 远小于样本量、相邻重复率高）。因此**在诊断通过之前，其样本不得用作
> 不确定性度量**。MCMC 只是唯一在术语上“允许”谈分布的模块，不代表它已经配得上这个资格。
> **默认关闭**；不开启时与 P0/P1 路径零漂移。
> 采样器同样**接受协作式取消**（`MCMCSampler.sample(stop_check=...)`，由 `api.rank` /
> `Ranker.run` 透传）：在链步之间检查、停止后返回已记录的样本，`run_metadata["mcmc_stopped_early"]`
> 回传是否真的被打断；若一条链**一个样本都没产出**就被取消，则 `RuntimeError` 而非静默返回空
> （见 8.9）。

## 8.3 目标类群剪裁

只关心物种树的某个子类群时，用 `--target-clade <标签>` 把排序聚焦到该类群：

```bash
# 默认：仅保留两端点都在类群 61 内的约束
maxtic-next species.tree constraints.tsv --target-clade 61

# 开启保守的"外部祖先映射到根"规则（谨慎使用）
maxtic-next species.tree constraints.tsv --target-clade 61 --target-clade-ancestor-map
```

- 默认策略：丢弃跨类群约束，**绝不静默引入虚假偏序**。
- `--target-clade-ancestor-map`：仅当外部端点可判定为类群祖先并落在 donor 侧时，才映射到根；
  旁支一律丢弃。默认关闭。
- 本步骤**不修改树拓扑 / 不重根**，只过滤 / 重映射约束端点。
- `--target-clade` 标签不存在时，`--dry-run` 与正式运行都会**报错**（不再等到深处 `KeyError`）。
- 提醒：被映射到根的约束在任何合法序中恒成立（根位次为 0），只会稀释权重分母、不提供信息。

## 8.4 预检（dry-run）

在正式排序前静态检查输入，**不消耗随机数、不修改数据**，可安全置于管线最前端：

```bash
maxtic-next species.tree constraints.tsv --dry-run
maxtic-next species.tree fam1.uml_rec --from ale --dry-run   # 真正解析上游输出
```

检查项（逐条对照 `src/maxtic_next/dry_run.py`）：

| 检查 | 严重度 |
|------|--------|
| 物种树可解析、有根 | error |
| 内部节点标签唯一 | error |
| **二叉树假设**（每个内部节点恰 2 个子节点） | error |
| **叶子名唯一** | error |
| `--target-clade` 标签存在于物种树 | error |
| 约束文件存在、可读、双格式可解析（拒三列逗号误用） | error |
| 约束权重/距离为有限实数且权重 ≥ 0（含文件名与行号） | error |
| 约束端点均在物种树中 | error |
| `--d > 0` 但约束**无** distance 列（选项将被忽略） | **error** |
| 物种树没有内部节点 | warning |
| 适配器模式：上游输出解析失败 / 零约束 / 端点全不命中 | error |
| 适配器模式：端点命中率低于 `--min-endpoint-hit-rate` | warning |

输出结构化报告并以 `结论：通过 / 存在致命问题，已中止 （error N, warning N）` 收尾。
合法输入（`examples/minitree.tree` + `examples/Cyano_CUTConstraints.tsv`）实测：

```
=== 预检（dry-run）报告 ===
未发现任何问题。
结论：通过 （error 0, warning 0）
```

> **退出码取自结构化字段** `run_metadata["dry_run_ok"]`，不再靠对报告文本做
> 子串匹配。预检通过退出 `0`，存在 error 退出 `1`。
> 适配器模式的预检**真正解析上游输出**：过去两行垃圾文本冒充 `.dtl` 会被判
> “通过”、去掉 `--dry-run` 后又在 0 总权重上产出一个“看似正常”的完整摘要；现在零约束与
> 端点全不命中都报 error。`--dry-run` 与 `-o` 同时给出时，`-o` 也照常写出约束文件。

### 8.4b 解析失败必定让预检判失败

预检在适配器模式下会真的跑一遍上游转换。**任何一次适配器解析失败**——损坏/截断的 XML、
读不了的文件、多成员归档——都转成 **error 级**报告项，并重算 `n_errors` / `ok`，
于是 CLI 以**退出码 1** 中止（旧实现在这条路径上追加了 error 却不重算 `ok`，
出现“一边打印 `[error]` 一边判定通过、退出码 0”的回归）。

实例：把 `examples/adapters/eccetera/FAM1.recphyloxml` 截断到 400 字节（`head -c 400`），
XML 因此未闭合：

```bash
$ head -c 400 examples/adapters/eccetera/FAM1.recphyloxml > truncated.recphyloxml
$ maxtic-next examples/adapters/ale/species.tree truncated.recphyloxml \
      --from eccetera --dry-run
=== 预检（dry-run）报告 ===
[error] 上游输出（格式 eccetera）解析失败：UpstreamParseError: /tmp/truncated.recphyloxml：
        recPhyloXML 解析失败（文件损坏 / 被截断 / 非 XML）（文件 /tmp/truncated.recphyloxml，
        第 13 行，第 11 列）；原始异常：ParseError: unclosed token: line 13, column 11
结论：存在致命问题，已中止 （error 1, warning 0）
$ echo $?
1
```

（错误全文见一次 `PYTHONPATH=src python3 -m maxtic_next … --from eccetera --dry-run`；
上面为了排版折了行。）对照实验：**未**截断的同一份文件预检退出 `0`
（`未发现任何问题。/ 结论：通过 （error 0, warning 0）`）。
Python 侧判定用 `result.run_metadata["dry_run_ok"]` / `dry_run_n_errors`，不要匹配文本。

## 8.5 增量打分（性能）

局部搜索的邻域移动是**区间旋转**。`--incremental` 只计算受影响区间的增量（O(b−a)），
而不是每次移动后全量重算 `value()`（O(|E|)），在大 |E| / 长搜索时显著加速，**结果不变**：

```bash
maxtic-next species.tree constraints.tsv --ls 600 --incremental
```

> **默认已开启**：增量与全量在 4 组输入 × 5 个种子 × 7 个产物字段上逐字节一致
> （由 `tests/test_output_equivalence.py` 看守），实测提速 3.1–8.2 倍。它只影响**性能**，
> 不改变搜索轨迹与结果；需要原版全量口径时用 `--no-incremental`。

## 8.5b 让局部搜索跨机器可复现

```bash
maxtic-next species.tree constraints.tsv --ls 600 --local-search-max-iters 2000000
```

- `--local-search-max-iters 0`（默认）= 只受 `--ls` 墙钟时长约束 → **迭代次数随机器负载而变**，
  同种子在不同机器上可能得到不同结果（程序会在 stdout 末尾打 WARNING 说明这点）。
- 设为正整数 N → 至多 N 次迭代后停止。随机流消耗确定，因此**固定 `--seed` 即确定**。
  给出 N 后 `--ls` 退化为“最多 also 这么久”的上界。
- 默认路径（`--ls 0`）不涉及此问题，始终可复现。

## 8.6 检查点 / 断点续跑

长局部搜索可能中途中断。`--checkpoint` 定期保存搜索状态，中断后再次运行同命令自动续跑：

```bash
maxtic-next species.tree constraints.tsv --ls 3600 \
    --checkpoint run.ckpt --checkpoint-interval 30
# 若中断，重复上面命令即从 run.ckpt 续跑
```

- 仅 `--local-search > 0` 时生效。
- `--checkpoint-interval`：自动保存间隔（秒，默认 60）。
- 检查点让**长搜索可续**，但它本身**不提供跨机器确定性**——那要靠
  `--local-search-max-iters`（8.5b）。二者互补：续跑解决“跑不完”，上限解决“不同机器
  跑出不同结果”。

> 注：ALE **解析阶段**另有按基因家族的文件级缓存（`--ale-cache-dir`），是独立的断点续传机制；
> 其缓存键含**物种树拓扑指纹**，因此标签集相同、拓扑不同的两棵树
> 不会再互相污染。

## 8.7 随机树显著性（p 值）

`--random-trees N` 生成 N 个**合法拓扑序**的随机排序，计算最终交付序的 value 与 similarity
的经验 p 值，并写出分布文件：

```bash
maxtic-next species.tree constraints.tsv --rd 1000
```

`--rd 50` 在蓝细菌样例上的实测输出：

```
values from  50  random orders 781.3100000000004 1058.9899999999996
pvalue of the found order: 0.0196078431372549 [(k+1)/(n+1) corrected; tested order = delivered best order]
similarity values from  50 random orders 0.4 0.9111111111111111
pvalue of the similarity with the input order: 0.0196078431372549 [(k+1)/(n+1) corrected]
```

两处已修正的统计口径：

1. **`(k+1)/(n+1)` 校正**（Phipson & Smyth 2010）。旧实现报 `k/n`，在 N=50 时最小非零 p 被
   压成 `0.0`，制造出“绝对显著”的假象；现在 50 次采样的最小 p 是 `1/51 ≈ 0.0196`。
   打印里显式标注 `[(k+1)/(n+1) corrected]`。
2. **检验对象是最终交付序**（局部搜索之后），而不是搜索前的启发式序——旧实现检验的序
   并不是交付给用户的那棵树。打印里标注 `tested order = delivered best order`。

分布文件是**第 4 个产物**：`<前缀>.mt.random_dist.tsv`（legacy 风格为
`<前缀>_distribution_random`），每行 `value similarity`，共 N 行。
随机序由 `--seed` 驱动，固定种子下逐字节可复现。

> p 值只表示“观测序比 N 个随机合法序更好”，**不**是后验概率，也不构成对排序正确性的证明。

## 8.8 随机化类型 `--random-type`

| 值 | 含义 |
|----|------|
| `0`（默认） | 保持数据原样 |
| `1` | 保留节点、随机化方向 |
| `2` | 完全随机化节点 |

取值域在 argparse 层强制为 `{0,1,2}`：旧实现把 `3`、`-1` 等**静默当作 0**，等于用未随机化的
数据做对照实验。`--random-type 1/2` 会改变数据本身，因此“与默认路径一致”无从谈起，
但同种子仍自洽可复现；仅在明确需要对照时使用。

## 8.9 协作式取消（`stop_check`）

长局部搜索与长 MCMC 链可以被**协作式取消**。这是一条 **API / GUI 层**的能力
（CLI 没有对应开关），形参为 `stop_check: () -> bool`：

| 插入点 | 签名 | 检查位置 |
|--------|------|----------|
| `api.rank(...)` | `stop_check=None` | 透传给 `Ranker.run` |
| `Ranker.run(...)` | `stop_check=None` | 透传给局部搜索与 MCMC 链 |
| `ranking.local_search.optimisation_locale(...)` | `stop_check=None` | **搜索主循环体第一条语句**，即每次迭代边界 |
| `robustness.mcmc.MCMCSampler.sample(...)` | `stop_check=None` | 每次链步之间（burn-in 与记录阶段都检查） |

**它保证什么**：

- 在**迭代边界**停止——不会在一次迭代中间被打断，因此交付的解永远是**自洽**的；
- **返回**截至取消点的历史最优序并照常产出三文件 + HTML 报告，**不抛异常**、
  不毁掉已算出的部分。取消成功时 stdout 末尾追加一条
  `NOTE: 本次运行按请求**取消**：局部搜索在第 N 次迭代边界停止（计划的搜索时长为 T 秒，
  未跑满）……`；
- 是否真的打断可机器判定：`Result.run_metadata["cancelled"]`
  （= 局部搜索 `stats_out["cancelled"]` **或** MCMC 的 `mcmc_stopped_early`）、
  `run_metadata["stop_check_requested"]`（谓词是否被传进来）、
  `run_metadata["mcmc_stopped_early"]`；
- **唯一的例外**：一条被取消且**尚未产出任何样本**的 MCMC 链会 `RuntimeError`
  （“链在产出第一个样本之前就被**取消**（stop_check），本次运行没有可交付的后验近似样本”），
  而不是静默返回“0 个样本的统计”——没有样本就没有可交付的量，程序拒绝编造。

**不传 `stop_check`（默认 `None`）时**，一切行为与不取消时**逐字节一致**，
`cancelled` 恒为 `False`。

```python
import itertools
from maxtic_next import rank

calls = itertools.count()
def stop_after_50():                    # 谓词：第 50 次检查之后返回真
    return next(calls) > 50

r = rank("species.tree", "constraints.tsv",
         local_search=30.0, temperature=60.0, stop_check=stop_after_50)
r.run_metadata["cancelled"]              # True
r.run_metadata["local_search_iterations"]  # 51 —— 计划 30 秒，实际在迭代边界停下
r.values["best"]                         # 截至取消点的最优值，照常交付
```

实测（同一份蓝细菌样例，`--ls 30`、`--t 60`、谓词在第 51 次检查转真）：
`cancelled=True`、迭代数 51、耗时远小于 30 秒、结果照常写出。

**Studio 的「取消 ■」按钮就绑在这条路径上**（Studio 仓库的 `maxtic_studio/cancellation.py`）：按钮 →
`RunEngine.cancel()` → `StopToken.cancel()`；令牌由 `cancellation_hook` 探测形参名后
以 `stop_check=token.should_stop` 注入 `api.rank`，因此取消会在迭代边界**真正打断计算**，
而不是过去那样“等跑完再丢弃结果”。若下层某天不再接受该形参，日志面板会显式打出
`gui/i18n.py` 的 `log_cancel_unsupported` 原句——“当前内核版本不接受取消回调，取消只能
在阶段边界生效（长 local-search 任务会跑完当前阶段后丢弃结果）。”，而不是假装能中断。

> 上一章：[07 · 流程封装与部署](07_workflows_deployment.md) ｜ 下一章：[09 · 常见问题与排错](09_faq_troubleshooting.md)
