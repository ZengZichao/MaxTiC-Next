"""ALE_undated ``.uml_rec`` 约束适配器。

忠实移植原 ``constraints_from_reconciliations.py``（Python 2）的核心逻辑：

* ``donnor_search`` / ``receptor_search``：在调和树（reconciliation tree）上沿父链/
  子树递归，解析 ``T@``（转移）与 ``D@``（重复）事件，抽取 donor / receptor 物种；
* 按基因家族聚合转移事件与调和事件的计数，做互反减法（reciprocal subtraction），
  再按 ``MINIMUM_SUPPORT_WITHIN_A_FAMILY``（单家族内最小支持度）过滤；
* 按 ``MINIMUM_FAMILY_SIZE``（最小基因家族规模，叶子数）跳过过小的家族；
* 产出统一 ``ConstraintSet``，每条约束 ``metadata`` 含 ``family`` / ``support`` /
  ``distance``（按物种树拓扑距离，等价原版把物种树所有枝长设为 1 后的距离）。

与原版一致的关键约定（语义校验，非正则）：

* 事件标注以 ``.`` 分隔，``T@donor->receptor`` 标记一次转移，``D@species`` 标记重复；
* ``donnor_search`` 沿父链上溯，遇根或 ``@T`` 标记则重置；
* 约束端点（donor / receptor）须为物种树中存在的标签。

可配置参数 ``MINIMUM_SUPPORT_WITHIN_A_FAMILY`` / ``MINIMUM_FAMILY_SIZE`` 作为适配器
属性暴露。为支持数千家族的断点续传，提供按基因家族文件级的 pickle 缓存。

**默认约束来源 = ``trf``（转移事件）—— **
--------------------------------------------------
ALE 官方仓库自己的 MaxTiC 集成文档
（ALE 源码包内 ``maxtic/README.md:131-169``）写得很明确：
"For each transfer detected by ALE, this script reports that the **father of the
donor** branch should be older than the **child of the receptor** branch"，其产物
``constraints_from_transfers`` 就是 MaxTiC 的直接输入（示例命令
``python MaxTiC.py species_tree_dated constraints_from_transfers ls=180``）。
也就是说，ALE↔MaxTiC 的官方口径用的是**转移事件**（``trf``，键
``parent(donor),receptor``），而本实现旧默认值 ``rec``（调和事件）会**静默换掉整个
约束集**（真实数据上 rec 仅 13995 条、总权重 5482.98，trf 与官方产物逐条一致、
总权重 6438.6）。故现默认 ``source="trf"``；两种口径都可选、都保留实现：

* ``"trf"``（默认，与 ALE 官方 MaxTiC 集成一致）：每个 ``T@donor->receptor`` 事件
  产出 ``parent(donor) -> receptor``；
* ``"rec"``：``constraints_from_reconciliations.py`` 里位于 ``T@`` 门控体内的
  "按调和（donnor_search 上溯 + receptor_search 下钻）"分支，产出
  ``donnor_search(node) -> receptor``。

``convert`` 结束时在 stderr 打印一行声明：本次由**哪种口径**产出约束、边数与总权重，
并把口径写进 ``ConstraintSet.diagnostics``（``source`` /
``donor_endpoint_convention = "parent_of_donor"``）。

**与四条"调和报告"路径的层级差**：本适配器取"**供体的父**"
（``parent_of_donor``）为 donor；RANGER-DTLx / ARTra / ecceTERA / AleRax 取
"**供体本身**"（``donor_itself``）。同一转移事件在两条路径下端点相差一个物种层级，
**不可直接互换**；``--from-auto`` / ``registry.merge_constraint_sets`` 混用两种约定时
会给出警告。两者均保持各自工具的原生约定，本实现不做层级换算。
详见 ``constraints/adapters/README_semantics.md``。

**采样分母与丢弃计数**：支持度分母取
``.uml_rec`` 自己打印的 ``<N> reconciled`` 行数（真实采样数）；该缺失时退化为本文件
实际解析到的调和树数，两者都拿不到时权重按**整数计数**输出并在诊断/stderr 声明。
每一次 ``continue``（donor 无法解析、receptor 是现存物种、支持度不足、家族规模不足、
自环）都被计数。文件级缓存键现含物种树**拓扑指纹**与版本号，跨拓扑复用不再可能。

并行化：默认使用 ``ProcessPoolExecutor``（进程级并行，绕过 GIL）；
可选 ``thread`` 回退到 ``ThreadPoolExecutor``。所有并行 worker 均为**模块级函数**，
参数仅含可 pickle 的简单类型（文件路径、配置参数、预计算的 ``parent_map`` 字典等），
结果为纯 dict/str/int 元组。失败时自动回退到顺序执行。
"""

import hashlib
import os
import pickle
from typing import Any, Dict, List, Optional, Tuple

from maxtic_next.constraints.constraint import Constraint, ConstraintSet
from maxtic_next.constraints.adapters._diagnostics import (
    CACHE_KEY_VERSION,
    CONVENTION_PARENT_OF_DONOR,
    WEIGHT_BLOCK_SUPPORT,
    WEIGHT_DECLARED_SUPPORT,
    WEIGHT_INTEGER_COUNT,
    finalize_adapter_diagnostics,
    new_diagnostics,
    reciprocal_subtract,
)
from maxtic_next.constraints.adapters._species import (
    species_labels_digest as _species_digest,
)
from maxtic_next.io.compression import logical_stem, read_text_lines
from maxtic_next.tree.tree import Tree


#: 默认约束来源：``trf``＝转移事件（ALE 官方 MaxTiC 集成口径）
DEFAULT_SOURCE = "trf"


# 与原版一致的默认阈值，同时作为可配置适配器参数暴露
MINIMUM_SUPPORT_WITHIN_A_FAMILY = 0.05
MINIMUM_FAMILY_SIZE = 5


# ----------------------------------------------------------------------
# 模块级函数（可 pickle，供 ProcessPoolExecutor 子进程调用）
#
# ``donnor_search`` / ``receptor_search`` 操作的是**调和树**（reconciliation
# tree），不依赖 ``ALEAdapter`` 实例状态，故可安全提取为模块级函数。
# ``_build_parent_map`` 从物种树预计算 ``label -> parent_label`` 字典，
# 使 worker 无需携带不可 pickle 的 ``Tree`` 对象。
# ----------------------------------------------------------------------


def _donnor_search(t: Tree, n: int):
    """模块级 ``donnor_search``：沿父链上溯，找到首个非 ``T@``/``D@`` 事件。

    等价 ``ALEAdapter.donnor_search``，但不依赖 ``self``，可被子进程直接调用。

    Args:
        t: 调和树。
        n: 起始节点 id。

    Returns:
        donor 物种标签（字符串）；若到达根或遇 ``@T`` 标记仍未找到，返回 -1。
    """
    current_node = n
    # 移植自原版的返回值形状：找到 donor 时是 str，找不到时是 -1（哨兵）。
    result: Any = -1
    while result == -1 and not t.is_root(current_node):
        current_node = t.get_parent(current_node)
        annot = t.get_bootstrap(current_node)
        event = annot[1:].split(".")[0]
        if event[:2] != "T@" and event[:2] != "D@":
            result = event
        if annot.find("@T") >= 0:
            current_node = t.get_root()
            result = -1
    return result


def _receptor_search(t: Tree, n: int) -> List[str]:
    """模块级 ``receptor_search``：递归收集受体物种标签。

    等价 ``ALEAdapter.receptor_search``，但不依赖 ``self``，可被子进程直接调用。

    Args:
        t: 调和树。
        n: 节点 id。

    Returns:
        受体物种标签列表。
    """
    if t.is_leaf(n):
        return []
    annot = t.get_bootstrap(n)
    if annot.find("T@") >= 0:
        return []
    elif annot.split(".")[1][:2] == "D@":
        children = t.get_children(n)
        return _receptor_search(t, children[0]) + _receptor_search(t, children[1])
    else:
        return [annot.split(".")[1]]


def _build_parent_map(species_tree: Tree) -> Dict[str, object]:
    """从物种树预计算 ``label -> parent_label`` 映射（可 pickle，供子进程使用）。

    等价 ``ALEAdapter.parent(label)`` 的查表版：根 / 不存在标签映射到 -1。
    叶子以 ``name`` 为键，内部节点以 ``bootstrap`` 为键，与 ``_find_id`` 的查找逻辑
    一致；缺失内部标签（空字符串）不进入映射，防止空端点通过校验。
    """
    parent_map: Dict[str, object] = {}
    for nid in species_tree.get_nodes():
        if species_tree.is_leaf(nid):
            label = species_tree.get_name(nid)
        else:
            label = species_tree.get_bootstrap(nid)
        if not label:
            continue
        if species_tree.is_root(nid):
            parent_map[label] = -1
        else:
            parent_map[label] = species_tree.get_bootstrap(species_tree.get_parent(nid))
    return parent_map


def _parent_map_identity(parent_map: Dict[str, object]) -> str:
    """由 ``parent_map`` 直接推出的物种树标识。

    当调用方未显式传 ``species_identity`` 时使用：指纹仍取"排序后的
    ``parent->child`` 标签对"摘要，而非标签集合摘要。
    """
    pairs = sorted(f"{p}->{c}" for c, p in parent_map.items() if p != -1)
    joined = "|".join(pairs)
    return (
        f"v{CACHE_KEY_VERSION}|{_species_digest(parent_map.keys())}|"
        f"{hashlib.sha256(joined.encode('utf-8')).hexdigest()[:16]}"
    )


def _process_family_worker(args: Tuple):
    """模块级 worker 函数（可 pickle，供 ``ProcessPoolExecutor`` 调用）。

    解析单个 ``.uml_rec`` 文件，返回
    ``(constraints_rec, constraints_trf, number, family, trees_seen, stats)``。
    所有参数均为简单可 pickle 类型；结果为纯 dict/str/int 元组。

    与 ``ALEAdapter._process_family`` 逐行一致，确保顺序/并行结果相同。

    Args:
        args: ``(rec_path, min_family_size, parent_map, extant_species,
               cache_dir, min_support[, species_identity])``

    Returns:
        ``(constraints_rec, constraints_trf, number, family, trees_seen, stats)``。
        ``number`` 为文件自己声明的 ``<N> reconciled`` 采样数（可能为 0），
        ``trees_seen`` 为实际解析到的调和树数，``stats`` 为诊断字典。
    """
    (rec_path, min_family_size, parent_map, extant_species, cache_dir, min_support) = tuple(
        args[:6]
    )
    species_identity = args[6] if len(args) > 6 else None

    # ---- 文件级缓存（断点续传，与 ALEAdapter._load_cache 一致）----
    cache_payload = None
    cache_file = None
    if cache_dir:
        st = os.stat(rec_path)
        # 缓存载荷依赖物种树的**拓扑**（parent_map → donor 端点），故键必须含
        # 拓扑指纹：仅标签摘要会让"标签集合相同、拓扑不同"的两棵树复用缓存，
        # 静默产出错误的 donor（真值 E、缓存给 F）。版本号
        # CACHE_KEY_VERSION 使旧的、不含拓扑的条目永不被读到。
        identity = species_identity or _parent_map_identity(parent_map)
        raw = (
            f"v{CACHE_KEY_VERSION}|{os.path.abspath(rec_path)}|{st.st_mtime}|"
            f"{st.st_size}|{min_support}|{min_family_size}|{identity}"
        )
        key = hashlib.sha256(raw.encode("utf-8")).hexdigest()
        cache_file = os.path.join(cache_dir, key + ".pkl")
        if os.path.exists(cache_file):
            try:
                with open(cache_file, "rb") as fh:
                    cache_payload = pickle.load(fh)
            except Exception:
                cache_payload = None
        if isinstance(cache_payload, tuple) and len(cache_payload) < 6:
            # 缺字段（4/5 元组）载荷：无采样数/诊断信息。键已升版，正常不会被读到；
            # 若用户手工混入，则据 transfers 反推并明确标注口径不可信。
            cache_payload = tuple(cache_payload) + (cache_payload[2], _legacy_stats(cache_payload))
    if cache_payload is not None:
        return cache_payload

    # ---- 解析（与 ALEAdapter._process_family 逐行一致）----
    # ：读取走 io.compression 的共享文本读取器，.gz / 单成员 .tar.gz 可用
    constraints_rec: Dict[str, int] = {}
    constraints_trf: Dict[str, int] = {}
    stats = new_diagnostics("ale", CONVENTION_PARENT_OF_DONOR)
    stats["transfers_unpaired"] = 0
    file_rec = read_text_lines(rec_path, encoding="utf-8")
    family = logical_stem(rec_path)
    stats["files_seen"] = 1

    number = 0
    trees_seen = 0
    family_skipped = 0
    i = 0
    while i < len(file_rec):
        words = file_rec[i].split()
        # 原版格式：``<N> reconciled`` 后跟 N 棵调和树（中间隔一行被跳过）
        if len(words) > 1 and words[1] == "reconciled":
            number = int(words[0])
            stats["declared_sample_counts"] = [number]
            for j in range(i + 2, i + 2 + number):
                if j >= len(file_rec):
                    break
                tree = Tree()
                tree.read_newick(file_rec[j])
                root = tree.get_root()
                leaves = tree.get_leaves(root)
                trees_seen += 1
                # 家族规模过滤（原版 len(leaves) > MINIMUM_FAMILY_SIZE）
                if len(leaves) > min_family_size:
                    for nid in tree.get_nodes():
                        if tree.is_leaf(nid):
                            annot = tree.get_name(nid)
                        else:
                            annot = tree.get_bootstrap(nid)
                        events = annot.split(".")
                        if events[0] == "":
                            del events[0]
                        # 原版控制流（constraints_from_reconciliations.py L152-211）：
                        # for e in events: if e[:2]=="T@": -> 1/ trf 约束 + 2/ rec 约束。
                        # 即 rec 约束块位于 T@ 门控**体内**，仅对标注了转移事件的节点执行
                        # （每个 T@ 事件执行一次），纯物种形成/重复节点不产生 rec 约束。
                        for e in events:
                            if e[:2] == "T@":
                                stats["transfers_seen"] += 1
                                # 1) 转移事件 -> constraints_trf（来自转移）
                                donnor = e[2:].split("->")[0]
                                # 原版复用同一个变量名：下面"调和事件"分支里它是 list[str]，
                                # 在"转移事件"分支里是 str。这里如实标注为 Any，不改控制流。
                                receptor: Any = e[2:].split("->")[1]
                                p = parent_map.get(donnor, -1)
                                if p != -1 and receptor not in extant_species:
                                    c = str(p) + "," + str(receptor)
                                    if str(p) == str(receptor):
                                        # 自环：仍计数，不静默消失
                                        stats["self_loop_transfers"] += 1
                                    stats["transfers_resolved"] += 1
                                    constraints_trf[c] = constraints_trf.get(c, 0) + 1
                                else:
                                    # 原版宽松校验：donor 无父（根/不存在）或受体是
                                    # 现存物种 ⇒ 丢弃，但**必须可归因**
                                    stats["dropped_by_label_miss"] += 1
                                    if (
                                        donnor not in parent_map
                                        and len(stats["unresolved_endpoint_samples"]) < 8
                                    ):
                                        stats["unresolved_endpoint_samples"].append(donnor)
                                # 2) 调和事件 -> constraints_rec（来自调和）
                                if tree.is_leaf(nid):
                                    receptor = []
                                else:
                                    if events[0][:2] != "T@":
                                        if events[0][:2] == "D@":
                                            species = events[0][2:]
                                        else:
                                            species = events[0]
                                        receptor = [species]
                                    else:
                                        children = tree.get_children(nid)
                                        child1 = children[0]
                                        if tree.is_leaf(child1):
                                            annot_child1 = (
                                                tree.get_name(child1).split(".")[-1].split("_")[0]
                                            )
                                        else:
                                            annot_child1 = tree.get_bootstrap(child1).split(".")[-1]
                                        if annot_child1[:2] == "T@":
                                            annot_child1 = annot_child1[2:].split("->")[0]
                                        if annot_child1[:2] == "D@":
                                            annot_child1 = annot_child1[2:]
                                        child2 = children[1]
                                        if tree.is_leaf(child2):
                                            annot_child2 = (
                                                tree.get_name(child2).split(".")[-1].split("_")[0]
                                            )
                                        else:
                                            annot_child2 = tree.get_bootstrap(child2).split(".")[-1]
                                        if annot_child2[:2] == "T@":
                                            annot_child2 = annot_child2[2:].split("->")[0]
                                        if annot_child2[:2] == "D@":
                                            annot_child2 = annot_child2[2:]
                                        if receptor == annot_child1 and donnor == annot_child2:
                                            receptor = _receptor_search(tree, child1)
                                        elif donnor == annot_child1 and receptor == annot_child2:
                                            receptor = _receptor_search(tree, child2)
                                        else:
                                            receptor = []
                                            stats["transfers_unpaired"] += 1
                                donnor = _donnor_search(tree, nid)
                                if donnor == -1:
                                    stats["dropped_by_label_miss"] += 1
                                for r in receptor:
                                    if donnor != -1 and r not in extant_species:
                                        key = str(donnor) + "," + str(r)
                                        constraints_rec[key] = constraints_rec.get(key, 0) + 1
                else:
                    family_skipped += 1
                    stats["dropped_by_family_size"] += 1
            i = i + number + 1
        i = i + 1

    stats["blocks_seen"] = 1
    stats["trees_seen"] = trees_seen
    stats["families_skipped"] = family_skipped
    # 家族规模探测：.uml_rec 总有叶子信息（get_leaves），故探测成功
    stats["family_size_probe_ok"] = trees_seen > 0 or number > 0

    payload = (constraints_rec, constraints_trf, number, family, trees_seen, stats)

    # ---- 缓存写入（与 ALEAdapter._save_cache 一致）----
    if cache_dir and cache_file is not None:
        os.makedirs(cache_dir, exist_ok=True)
        with open(cache_file, "wb") as fh:
            pickle.dump(payload, fh)

    return payload


def _legacy_stats(payload) -> Dict:
    """缺字段缓存载荷（4/5 元组）的诊断占位：口径按声明采样数推定并标注不可信。"""
    trf, number = payload[1], payload[2]
    stats = new_diagnostics("ale", CONVENTION_PARENT_OF_DONOR)
    stats["blocks_seen"] = 1
    stats["declared_sample_counts"] = [number]
    stats["sample_denominator"] = number or None
    stats["weight_semantics"] = WEIGHT_DECLARED_SUPPORT if number else WEIGHT_INTEGER_COUNT
    stats["transfers_seen"] = sum(trf.values())
    stats["transfers_resolved"] = sum(trf.values())
    stats["files_seen"] = 1
    stats["cache_legacy_payload"] = True
    stats["warnings"].append(
        "读到的是**缺字段 ALE 缓存载荷**（不含采样数/拓扑标识）：其权重口径不可信，"
        "建议换一个 --ale-cache-dir 重跑。"
    )
    return stats


class ALEAdapter:
    """ALE_undated ``.uml_rec`` 适配器：调和结果 -> 统一 ``ConstraintSet``。

    ``source`` 取 ``"trf"``（默认，ALE 官方 MaxTiC 集成的转移事件口径）
    或 ``"rec"``（调和事件口径）。两种口径的 donor 端点都在"**供体的父**"这一层
    （``donor_endpoint_convention = "parent_of_donor"``），与 RANGER-DTLx / ARTra /
    ecceTERA / AleRax 的"供体本身"相差一个物种层级，不可直接互换。
    """

    def __init__(
        self,
        species_tree: Tree,
        min_support: float = MINIMUM_SUPPORT_WITHIN_A_FAMILY,
        min_family_size: int = MINIMUM_FAMILY_SIZE,
        cache_dir: Optional[str] = None,
        source: str = DEFAULT_SOURCE,
        quiet: bool = False,
    ) -> None:
        """初始化适配器。

        Args:
            species_tree: 物种树（内部节点标签在 bootstrap 字段）。
            min_support: 单基因家族内最小支持度阈值（默认 0.05）。
            min_family_size: 参与计算的最小基因家族规模（调和树叶子数，默认 5）。
            cache_dir: 基因家族文件级缓存目录；``None`` 表示不缓存。缓存键含
                物种树**拓扑**指纹与版本号，跨拓扑复用不可能。
            source: 约束来源，``"trf"``（**默认**：转移事件，即 ALE 官方
                ``constraints_from_transfers`` 口径）或 ``"rec"``
                （调和事件，``constraints_from_reconciliations.py`` 的
                ``donnor_search``/``receptor_search`` 分支）。
            quiet: 不在 stderr 打印口径声明与摘要（诊断仍随结果返回）。
        """
        if source not in ("rec", "trf"):
            raise ValueError(f"source 必须为 'rec' 或 'trf'，收到 {source!r}")
        self.species_tree = species_tree
        self.min_support = min_support
        self.min_family_size = min_family_size
        self.cache_dir = cache_dir
        self.source = source
        self.quiet = quiet
        # 现存物种（物种树叶子名），供端点合法性判断
        self.extant_species = set(species_tree.get_leaves_names())
        # 惰性缓存：label -> parent_label 映射（供 ProcessPoolExecutor worker 使用）
        self._parent_map_cache: Optional[Dict[str, object]] = None

    def _get_parent_map(self) -> Dict[str, object]:
        """惰性构建并缓存 ``label -> parent_label`` 映射（供 worker 函数使用）。"""
        if self._parent_map_cache is None:
            self._parent_map_cache = _build_parent_map(self.species_tree)
        return self._parent_map_cache

    # ------------------------------------------------------------------
    # 物种树辅助：标签 -> 节点 id、父标签、拓扑距离
    # （委托共享 ``_species`` 模块：O(1) label_to_id 查表 + 统一距离语义）
    # ------------------------------------------------------------------
    def _find_id(self, label: str) -> Optional[int]:
        """按标签（叶子名或内部 bootstrap）在物种树中查找节点 id。"""
        from maxtic_next.constraints.adapters._species import (
            find_label_id as _find_label_id,
        )

        return _find_label_id(self.species_tree, label)

    def parent(self, x: str):
        """等价原版 ``parent(x)``：返回标签 ``x`` 的父节点 bootstrap 标签，根返回 -1。

        若 ``x`` 在物种树中不存在，返回 -1。
        """
        nid = self._find_id(x)
        if nid is None:
            return -1
        if self.species_tree.is_root(nid):
            return -1
        return self.species_tree.get_bootstrap(self.species_tree.get_parent(nid))

    def distance_from(self, x: str, y: str) -> Optional[float]:
        """等价原版 ``distance_from(x, y)``：物种树上 ``x`` 与 ``y`` 的拓扑距离。

        原版先把物种树所有枝长设为 1.0，故此处直接按边数（unit branch length）计算，
        不修改调用方传入的物种树（避免影响后续 ``order_from_tree`` 的真实距离语义）。

        任一标签不在物种树中时返回 ``None``（"距离未知"），不再返回 0.0
        ——：0.0 与"供体=受体的真距离 0"不可区分，且会被
        ``--min-transfer-distance`` 当成合法距离而把约束静默删除。
        """
        from maxtic_next.constraints.adapters._species import (
            distance_from as _shared_distance_from,
        )

        return _shared_distance_from(self.species_tree, x, y)

    # ------------------------------------------------------------------
    # 调和树事件解析（忠实移植 donnor_search / receptor_search）
    # 注意：以下函数操作的是**调和树** ``t``（reconciliation tree），不是物种树。
    # 委托至模块级函数，使 ProcessPoolExecutor worker 可直接调用同一逻辑。
    # ------------------------------------------------------------------
    def donnor_search(self, t: Tree, n: int):
        """等价原版 ``donnor_search(t, n)``：沿父链上溯，找到首个非 ``T@``/``D@`` 事件。

        委托至模块级 ``_donnor_search``（该函数不依赖 ``self``，可被子进程直接调用）。

        Args:
            t: 调和树。
            n: 起始节点 id。

        Returns:
            donor 物种标签（字符串）；若到达根或遇 ``@T`` 标记仍未找到，返回 -1。
        """
        return _donnor_search(t, n)

    def receptor_search(self, t: Tree, n: int) -> List[str]:
        """等价原版 ``receptor_search(t, n)``：递归收集受体物种标签。

        委托至模块级 ``_receptor_search``（该函数不依赖 ``self``，可被子进程直接调用）。

        Args:
            t: 调和树。
            n: 节点 id。

        Returns:
            受体物种标签列表（叶子返回空列表；``T@`` 节点返回空列表；``D@`` 节点
            递归合并两子树的受体）。
        """
        return _receptor_search(t, n)

    # ------------------------------------------------------------------
    # 单基因家族解析（含文件级缓存）
    # ------------------------------------------------------------------
    def _cache_key(self, rec_path: str) -> str:
        """文件路径 + mtime + size + 参数 + **物种树拓扑标识** → 缓存键。

        ：物种树标识必须含拓扑指纹（而不只是标签集合摘要），否则两棵
        "标签集合相同、父子关系不同"的树会复用同一目录里的缓存，产出错误的 donor；
        ``CACHE_KEY_VERSION`` 前缀使任何无拓扑信息的条目永不被读到。
        """
        from maxtic_next.constraints.adapters._species import (
            species_cache_identity as _species_cache_identity,
        )

        st = os.stat(rec_path)
        raw = (
            f"v{CACHE_KEY_VERSION}|{os.path.abspath(rec_path)}|{st.st_mtime}|"
            f"{st.st_size}|{self.min_support}|{self.min_family_size}|"
            f"{_species_cache_identity(self.species_tree)}"
        )
        return hashlib.sha256(raw.encode("utf-8")).hexdigest()

    def _load_cache(self, rec_path: str) -> Optional[Tuple]:
        if not self.cache_dir:
            return None
        key = self._cache_key(rec_path)
        cache_file = os.path.join(self.cache_dir, key + ".pkl")
        if os.path.exists(cache_file):
            try:
                with open(cache_file, "rb") as fh:
                    payload = pickle.load(fh)
            except Exception:
                return None
            if isinstance(payload, tuple) and len(payload) < 6:
                return tuple(payload) + (payload[2], _legacy_stats(payload))
            return payload
        return None

    def _save_cache(self, rec_path: str, payload: Tuple) -> None:
        if not self.cache_dir:
            return
        os.makedirs(self.cache_dir, exist_ok=True)
        key = self._cache_key(rec_path)
        cache_file = os.path.join(self.cache_dir, key + ".pkl")
        with open(cache_file, "wb") as fh:
            pickle.dump(payload, fh)

    def _process_family(self, rec_path: str) -> Tuple:
        """解析单个 ``.uml_rec`` 文件。

        返回 ``(constraints_rec, constraints_trf, number, family, trees_seen, stats)``
        （``number`` 为文件声明的采样数、``trees_seen`` 为实际解析到的树数）。
        结果经文件级缓存（断点续传）。

        委托至模块级 ``_process_family_worker``，与 ``ProcessPoolExecutor`` 并行路径
        使用完全相同的解析逻辑，确保顺序/并行结果逐字节一致。
        """
        return _process_family_worker(self._worker_args(rec_path))

    def _worker_args(self, rec_path: str) -> Tuple:
        """构造 worker 参数元组（含物种树拓扑标识，供缓存键使用）。"""
        from maxtic_next.constraints.adapters._species import (
            species_cache_identity as _species_cache_identity,
        )

        return (
            rec_path,
            self.min_family_size,
            self._get_parent_map(),
            self.extant_species,
            self.cache_dir,
            self.min_support,
            _species_cache_identity(self.species_tree),
        )

    @staticmethod
    def _reciprocal_subtract(counts: Dict[str, int]) -> None:
        """等价原版互反减法（委托共用实现）。

        与原版唯一的偏离：``"X,X"`` 自环键不再"自己减自己"归零；
        自环由 ``self_loop_transfers`` 计数上报，并交给下游按原版
        "uninformative (to itself)" 分支处理。
        """
        reciprocal_subtract(counts)

    # ------------------------------------------------------------------
    # 单家族结果合并（与顺序执行逐字节一致）
    # ------------------------------------------------------------------
    def _merge_payload(self, cset: ConstraintSet, payload: Tuple) -> None:
        """把一个家族的解析结果合并进 ``cset``（顺序无关：``ConstraintSet.add`` 按边键聚合）。

        支持度分母的确定次序：

        1. 文件自己声明的 ``<N> reconciled`` 采样数；
        2. 退化为本文件实际解析到的调和树数；
        3. 两者皆无 ⇒ 权重按**整数计数**输出（口径 ``integer_count``），
           并在诊断 + stderr 中声明"min_support 此时不是 [0,1] 支持度阈值"。

        Args:
            cset: 待累加的约束集。
            payload: ``(constraints_rec, constraints_trf, number, family
                [, trees_seen[, stats]])``。
        """
        constraints_rec, constraints_trf, number = payload[0], payload[1], payload[2]
        family = payload[3] if len(payload) > 3 else ""
        trees_seen = payload[4] if len(payload) > 4 else number
        if len(payload) > 5 and isinstance(payload[5], dict):
            stats = payload[5]
        else:
            stats = _legacy_stats(payload)
        stats["files_seen"] = 1
        cset.diagnostics.setdefault("files", []).append(stats)
        stats["source"] = self.source

        # 分母与口径（不再"没声明就整族丢弃"）
        denominator = int(number or 0) or int(trees_seen or 0)
        if denominator > 0:
            stats["sample_denominator"] = denominator
            stats["weight_semantics"] = (
                WEIGHT_DECLARED_SUPPORT if int(number or 0) > 0 else WEIGHT_BLOCK_SUPPORT
            )
            if not int(number or 0):
                stats["warnings"].append(
                    f"{family or '某家族'}：``.uml_rec`` 未给出 ``<N> reconciled`` "
                    f"头行，改用实际解析到的 {denominator} 棵调和树作分母。"
                )
        else:
            stats["sample_denominator"] = None
            stats["weight_semantics"] = WEIGHT_INTEGER_COUNT
            stats["warnings"].append(
                f"{family or '某家族'}：既无 ``<N reconciled>`` 声明也无解析到的"
                "调和树 ⇒ 权重按**整数计数**输出（不是 [0,1] 支持度）。"
                "旧实现在此静默丢弃该家族的全部约束。"
            )

        # 互反减法（等价原版，按家族分别处理）
        self._reciprocal_subtract(constraints_trf)
        self._reciprocal_subtract(constraints_rec)
        # 选择约束来源（rec / trf）
        chosen = constraints_rec if self.source == "rec" else constraints_trf
        for key, cnt in chosen.items():
            if cnt <= 0:
                stats["dropped_by_support"] += 1
                continue
            if denominator > 0:
                support = cnt / float(denominator)
                keep = support > self.min_support
            else:
                # 计数口径：阈值不具支持度语义（已在诊断/stderr 声明）
                support = float(cnt)
                keep = cnt > 0
            if not keep:
                stats["dropped_by_support"] += 1
                continue
            donnor, receptor = key.split(",")
            # 自环事件的距离置 None，与其余四个适配器一致：若写 0.0，则在**默认**
            # --min-transfer-distance 0 下 "dist > d" 判假，donor == receptor 的事件
            # 被丢弃，永远进不了原版口径的 "to itself" 统计。
            # None = 无距离信息，与文本输入路径一致。
            distance = None if donnor == receptor else self.distance_from(receptor, donnor)
            cset.add(
                Constraint(
                    donor=donnor,
                    receptor=receptor,
                    weight=float(support),
                    metadata={
                        "family": family,
                        "support": support,
                        "distance": distance,
                        "source": self.source,
                        "weight_semantics": stats.get("weight_semantics"),
                        "sample_denominator": denominator or None,
                    },
                )
            )
            stats["constraints_kept"] += 1

    # ------------------------------------------------------------------
    # 对外主入口
    # ------------------------------------------------------------------
    def convert(
        self,
        rec_files: List[str],
        output_path: Optional[str] = None,
        max_workers: Optional[int] = None,
        parallel: str = "process",
    ) -> ConstraintSet:
        """把若干 ``.uml_rec`` 文件转换为 ``ConstraintSet``。

        各文件解析（``_process_family``）相互独立，默认用进程池并行加速
        （``concurrent.futures.ProcessPoolExecutor``，绕过 GIL 实现真正的多核并行）；
        可选 ``parallel="thread"`` 回退到 ``ThreadPoolExecutor``。结果在**收集后按文件
        顺序**合并，由于 ``ConstraintSet.add`` 按边键聚合，合并顺序不影响最终的
        边键 + 权重集合，故并行结果与顺序执行**逐字节一致**。

        进程池模式下，所有 worker 参数均为可 pickle 的简单类型（文件路径、配置参数、
        预计算的 ``parent_map`` 字典、``extant_species`` 集合），不传递 ``Tree`` 对象
        或绑定方法。若并行执行器抛错，自动回退为顺序执行。

        Args:
            rec_files: ``.uml_rec`` 文件路径列表。
            output_path: 可选，若提供则把约束以逗号格式写入该文件（便于调试/复用）。
            max_workers: 并行进程/线程数；``None`` 时取 ``min(8, len(rec_files), cpu_count)``
                且至少 1。若并行执行器抛错，自动回退为顺序执行。
            parallel: 并行模式，``"process"``（默认，``ProcessPoolExecutor``）或
                ``"thread"``（``ThreadPoolExecutor``）。

        Returns:
            统一 ``ConstraintSet``，每条约束 ``metadata`` 含 ``family`` / ``support`` /
            ``distance`` / ``source``；集合的 ``diagnostics`` 含口径（``source``、
            ``donor_endpoint_convention``、``weight_semantics``、分母）与全部丢弃计数，
            并有一行 stderr 声明"本次由哪种口径产出、边数与总权重"。
        """
        if not rec_files:
            cset = ConstraintSet(diagnostics=new_diagnostics("ale", CONVENTION_PARENT_OF_DONOR))
            if output_path is not None:
                self._write_constraints_file(cset, output_path)
            return cset

        # 预计算可 pickle 的查表数据（进程池 worker 需要）
        parent_map = self._get_parent_map()
        extant_species = self.extant_species
        from maxtic_next.constraints.adapters._species import (
            species_cache_identity as _species_cache_identity,
        )

        identity = _species_cache_identity(self.species_tree)
        worker_args = [
            (
                rf,
                self.min_family_size,
                parent_map,
                extant_species,
                self.cache_dir,
                self.min_support,
                identity,
            )
            for rf in rec_files
        ]

        # 解析阶段：并行优先，失败回退顺序
        payloads: List[Optional[Tuple]] = [None] * len(rec_files)

        executor_cls: Optional[type] = None  # ProcessPoolExecutor / ThreadPoolExecutor
        mp_context: Optional[Any] = None
        if parallel == "process":
            try:
                from concurrent.futures import ProcessPoolExecutor
                import multiprocessing

                executor_cls = ProcessPoolExecutor
                # 在 macOS 上默认 start method 为 spawn，会重新导入 __main__ 模块，
                # 导致在 pytest 等环境下子进程重复执行主模块而死锁。
                # 优先使用 fork（POSIX 可用，不重导入 __main__），Windows 回退默认。
                try:
                    mp_context = multiprocessing.get_context("fork")
                except ValueError:
                    mp_context = None  # Windows 无 fork，使用默认 spawn
            except Exception:
                executor_cls = None
        elif parallel == "thread":
            try:
                from concurrent.futures import ThreadPoolExecutor

                executor_cls = ThreadPoolExecutor
            except Exception:
                executor_cls = None

        if executor_cls is not None:
            try:
                if max_workers is None:
                    max_workers = min(os.cpu_count() or 1, 8, len(rec_files))
                    max_workers = max(1, max_workers)
                kwargs: Dict[str, Any] = {"max_workers": max_workers}
                if mp_context is not None:
                    kwargs["mp_context"] = mp_context
                with executor_cls(**kwargs) as executor:
                    future_to_idx = {
                        executor.submit(_process_family_worker, arg): i
                        for i, arg in enumerate(worker_args)
                    }
                    for future in future_to_idx:
                        idx = future_to_idx[future]
                        try:
                            payloads[idx] = future.result()
                        except Exception:
                            # 单个文件失败：回退为顺序重算该文件
                            payloads[idx] = self._process_family(rec_files[idx])
            except Exception:
                # 进程/线程池不可用：完全顺序解析
                payloads = [self._process_family(rf) for rf in rec_files]
        else:
            # 未知 parallel 值或导入失败：顺序解析
            payloads = [self._process_family(rf) for rf in rec_files]

        # 合并阶段：严格按文件顺序（与顺序执行一致；顺序无关仅指边键聚合不受顺序影响）
        cset = ConstraintSet(diagnostics=new_diagnostics("ale", CONVENTION_PARENT_OF_DONOR))
        for payload in payloads:
            if payload is not None:
                self._merge_payload(cset, payload)

        # ：一行 stderr 声明"本次约束由哪种口径产出"+ 边数 + 总权重
        edges = cset.informative_edges()
        total_weight = float(sum(edges.values()))
        notice = (
            f"[ale] 约束口径：source={self.source}"
            f"（{'trf=转移事件 parent(donor)->receptor，ALE 官方 MaxTiC 集成口径' if self.source == 'trf' else 'rec=调和事件 donnor_search->receptor_search'}）"
            f"，供体端点约定={CONVENTION_PARENT_OF_DONOR}，"
            f"产出边数={len(edges)}，总权重={total_weight}"
        )
        cset.diagnostics["source"] = self.source
        cset.diagnostics["total_weight"] = total_weight
        cset.diagnostics["edges_kept"] = len(edges)
        if output_path is not None:
            self._write_constraints_file(cset, output_path)
        finalize_adapter_diagnostics(
            cset,
            0.0,
            set(self._get_parent_map().keys()) | self.extant_species,
            "ale",
            quiet=self.quiet,
            files=len(rec_files),
            notice=notice,
            raise_on_zero_hit=False,
        )
        return cset

    @staticmethod
    def _write_constraints_file(cset: ConstraintSet, path: str) -> None:
        """把 ``ConstraintSet`` 以逗号格式（family,donor,receptor,weight,distance）写出。"""
        with open(path, "w", encoding="utf-8") as fh:
            fh.write("#family,donor,receptor,weight,distance\n")
            for c in cset.constraints:
                dist = c.metadata.get("distance")
                dist_str = "" if dist is None else str(dist)
                fam = c.metadata.get("family") or ""
                fh.write(f"{fam},{c.donor},{c.receptor},{c.weight},{dist_str}\n")
        return None


def convert_from_ale(
    species_tree: Tree,
    rec_files: List[str],
    min_support: float = MINIMUM_SUPPORT_WITHIN_A_FAMILY,
    min_family_size: int = MINIMUM_FAMILY_SIZE,
    cache_dir: Optional[str] = None,
    source: str = DEFAULT_SOURCE,
    output_path: Optional[str] = None,
    max_workers: Optional[int] = None,
    parallel: str = "process",
    quiet: bool = False,
) -> ConstraintSet:
    """便捷函数：把 ALE ``.uml_rec`` 文件列表转换为 ``ConstraintSet``。

    Args:
        species_tree: 物种树。
        rec_files: ``.uml_rec`` 文件路径列表。
        min_support: 单家族内最小支持度阈值。
        min_family_size: 最小基因家族规模。
        cache_dir: 文件级缓存目录（断点续传；键含物种树拓扑指纹）。
        source: 约束来源 ``"trf"``（**默认**，ALE 官方 MaxTiC 集成口径）/ ``"rec"``（调和事件口径）。
        output_path: 可选约束写出路径。
        max_workers: 并行解析进程/线程数（``None`` 取默认；见 ``ALEAdapter.convert``）。
        parallel: 并行模式，``"process"``（默认）或 ``"thread"``。
        quiet: 不打印 stderr 的口径声明与摘要。

    Returns:
        统一 ``ConstraintSet``。
    """
    adapter = ALEAdapter(
        species_tree,
        min_support=min_support,
        min_family_size=min_family_size,
        cache_dir=cache_dir,
        source=source,
        quiet=quiet,
    )
    return adapter.convert(
        rec_files, output_path=output_path, max_workers=max_workers, parallel=parallel
    )
