"""全局常量与默认参数配置。

这些常量与原版 ``MaxTiC.py`` 保持一致，是算法等价性的基础。
``MAX_NUMBER`` 用作"无穷 / 去信息（uninformative）"哨兵值：一旦某条边被标记为
树种系边（descendant edge）或冲突边，其权重被置为 ``MAX_NUMBER``，从而在 ``value()``
中不会被错误地计入反馈弧集权重（除非排序违反拓扑，那会被正确计入）。
"""

# 等价于原版 MaxTiC.py 的 MAX_NUMBER 哨兵值
MAX_NUMBER = 10000000000

# 偏序产物中"哨兵边（树种系边）"的判定阈值。
# 旧实现写死 ``edge[e] < 100000``，比 MAX_NUMBER 小 5 个数量级，会把权重 >= 1e5 的
# 真实信息性约束静默排除出 partial_order。现统一以 MAX_NUMBER 判定。
PARTIAL_ORDER_SENTINEL_EXCLUSIVE = MAX_NUMBER

# 局部搜索（Metropolis）默认温度；原版 TEMPERATURE 默认值
DEFAULT_TEMPERATURE = 0.001

# 随机化类型：0=保持数据原样，1=保留节点随机化方向，2=完全随机化节点
DEFAULT_RANDOM_TYPE = 0
RANDOM_TYPE_CHOICES = (0, 1, 2)

# 局部搜索时长（秒）；原版 TIME_FOR_SEARCH，默认 0（关闭）
DEFAULT_TIME_FOR_SEARCH = 0.0

# 最小转移距离阈值（按 phylogenetic distance 列过滤）；原版 MIN_TRANSFER_DIST，默认 0
DEFAULT_MIN_TRANSFER_DIST = 0

# 约束权重阈值比例；原版 THRESHOLD_CONSTRAINTS，默认 0.0（不过滤）
# 语义为"删除累计权重占比最低的约束"，故取值域为闭区间 [0, 1]。
DEFAULT_THRESHOLD_CONSTRAINTS = 0.0
THRESHOLD_CONSTRAINTS_MIN = 0.0
THRESHOLD_CONSTRAINTS_MAX = 1.0

# 随机树采样数量；原版 RANDOM，默认 0（关闭）
DEFAULT_RANDOM_TREES = 0

# 近优解收集器容量：局部搜索期间保留多少个**去重**近优排序，用于
# "基于局部搜索访问解的稳健性/敏感性摘要"。原版与本项目的历史默认都是 50；现在
# API（``api.rank(top_k=...)``）与 CLI（``--near-optimal-top-k``）都可覆盖它。
# 单一事实来源在此处：``robustness.sensitivity.DEFAULT_TOP_K`` 即本常量的别名。
DEFAULT_NEAR_OPTIMAL_TOP_K = 50

# 上游适配器端点命中率告警阈值：低于该比例 stderr 警告，
# 命中率为 0 直接报错。CLI ``--min-endpoint-hit-rate`` 与 API
# ``adapter_min_endpoint_hit_rate`` 共用此默认值（放在 config 里以便 CLI/API 无需
# 提前导入适配器层即可引用；``adapters._diagnostics`` 从本模块再导出）。
DEFAULT_MIN_ENDPOINT_HIT_RATE = 0.5

# 默认随机种子，确保无参运行可复现；原版无种子（全局随机），重写后固定为 42
DEFAULT_SEED = 42

# ---- MCMC（Metropolis–Hastings over linear extensions）默认参数----
# 温度 sentinel：0.0 == "auto"，有效温度 = max(total_weight, 1) / MCMC_TEMPERATURE_DIVISOR
MCMC_TEMPERATURE_AUTO = 0.0
MCMC_TEMPERATURE_DIVISOR = 100.0
DEFAULT_MCMC_ITERS = 1000
DEFAULT_MCMC_BURN_IN_FRACTION = 0.5  # burn-in = iters 的一半（分位数取自经验样本前段）
DEFAULT_MCMC_THIN_DIVISOR = 100  # thin = max(1, iters // 100)
MCMC_STATUS_NOTE = (
    "Metropolis–Hastings over linear extensions "
    "(preliminary; convergence diagnostics not validated)"
)

# 输出文件名后缀 —— legacy 风格（严格与原版 MaxTiC.py 一致，便于逐字节复现验证）
SUFFIX_FILTERED = "_MT_output_filtered_list_of_weighted_informative_constraints"
SUFFIX_CONFLICTING = "_MT_output_list_of_constraints_conflicting_with_best_order"
SUFFIX_PARTIAL_ORDER = "_MT_output_partial_order"
SUFFIX_RANDOM_DIST = "_distribution_random"

# 输出文件名后缀 —— short 风格（默认；简洁、带 .tsv 扩展名，便于识别为数据文件）。
# 三文件内容与 legacy **逐字节一致**，仅文件名不同（需要与原版同名产物时用 legacy）。
SUFFIX_FILTERED_SHORT = ".mt.informative.tsv"
SUFFIX_CONFLICTING_SHORT = ".mt.conflicts.tsv"
SUFFIX_PARTIAL_ORDER_SHORT = ".mt.partial_order.tsv"
SUFFIX_RANDOM_DIST_SHORT = ".mt.random_dist.tsv"

# 默认输出命名风格："short"（简洁，默认）或 "legacy"（与原版逐字节等价的长名）
DEFAULT_OUTPUT_STYLE = "short"

# 风格 -> (filtered, conflicting, partial_order, random_dist) 后缀四元组
OUTPUT_SUFFIXES = {
    "legacy": (SUFFIX_FILTERED, SUFFIX_CONFLICTING, SUFFIX_PARTIAL_ORDER, SUFFIX_RANDOM_DIST),
    "short": (
        SUFFIX_FILTERED_SHORT,
        SUFFIX_CONFLICTING_SHORT,
        SUFFIX_PARTIAL_ORDER_SHORT,
        SUFFIX_RANDOM_DIST_SHORT,
    ),
}
