"""约束转换适配器子包（多软件兼容框架）。

把 5 个第三方调和/转移推断工具的输出转换为 MaxTiC-Next 统一的
``ConstraintSet``（``Constraint(donor, receptor, weight, metadata)``）。所有适配器
均已按各工具**真实输出格式**（经源码/样本求证）实现，并接入 ``registry`` 的
插件式注册与自动检测：

* ``ale``：ALE_undated ``.uml_rec`` 采样调和（``T@``/``D@`` NHX 事件）；
* ``ranger``：RANGER-DTLx 调和报告（``Mapping-->`` / ``Recipient-->``）；
* ``eccetera``：ecceTERA **recPhyloXML**（``branchingOut`` / ``transferBack``）；
* ``artra``：ARTra ``output.txt``（Replacing / Additive Transfer，与 RANGER 同形）；
* ``alerax``：AleRax ``reconciliations/summaries/`` 下转移频率文件
  （AleRax 早期版本 ``*_transfers.txt``，现行版本 ``*_meanTransfers.txt``；``donor receptor freq``）。

端点校验语义与两个不可通约的口径（重要，混用多工具时须知；详见
``README_semantics.md``）：

* **供体端点层级**：ALE 的 donor 是"**供体的父**"
  （``parent(donnor)`` / ``donnor_search``，诊断值 ``parent_of_donor``）；
  RANGER-DTLx / ARTra / ecceTERA / AleRax 的 donor 是"**供体本身**"
  （``Mapping-->`` / ``branchingOut speciesLocation`` / 表文件第一列，诊断值
  ``donor_itself``）。同一生物学事件在两条路径下相差**一个物种层级**，端点不可
  直接互换；``registry.merge_constraint_sets`` / ``convert_auto(allow_mixed=True)``
  检测到混用会在 stderr 警告并记入 ``diagnostics["convention_conflict"]``。
  各工具均保持**自己的原生约定**，本包不做层级换算。
* **权重口径**：ALE 用其 ``<N> reconciled`` 声明的采样数作分母；
  RANGER-DTLx / ARTra 用报告自己声明的 ``Total number of optimal solutions: N``；
  分母不可确立时权重退化为**整数计数**（``integer_count``）并在诊断/stderr 声明，
  此时 ``min_support`` 不再是 [0,1] 支持度阈值。AleRax 的权重是上游预聚合的
  后验频率（``posterior_frequency``），与前三者量纲不同。
* **端点校验**：RANGER / ecceTERA / ARTra / AleRax 要求 donor 与 receptor **均**可
  解析进用户物种树标签集（数字物种节点 ID 会按 ecceTERA 的自底向上逐层编号规则
  反解），否则丢弃并计数；ALE 忠实移植原版的宽松校验——只要求 donor
  能解析出父节点、且 receptor 不是现存物种（叶子）名。
* **不再静默**：所有丢弃都计入
  ``ConstraintSet.diagnostics``（``blocks_seen`` / ``transfers_seen`` /
  ``constraints_kept`` / ``dropped_by_label_miss`` / ``dropped_by_support`` /
  ``dropped_by_family_size`` / ``self_loop_transfers`` …）并打印一行 stderr 摘要；
  端点命中率为 0 抛 ``LabelMismatchError``，XML 解析失败抛 ``UpstreamParseError``。
* **ALE 默认口径**：``source`` 默认为 ``"trf"``（转移事件），即 ALE
  官方 MaxTiC 集成 ``constraints_from_transfers`` 所规定的口径；``"rec"`` 仍可选。
* **文件级缓存键**：含物种树**拓扑指纹**与版本号，跨拓扑复用不再可能。

``registry`` 提供 ``detect_format`` 自动检测与 ``convert`` / ``convert_auto`` 分发，
支持 ``--from <tool>`` 与 ``--from-auto``。
"""

from maxtic_next.constraints.adapters.ale import (
    ALEAdapter,
    DEFAULT_SOURCE as ALE_DEFAULT_SOURCE,
    convert_from_ale,
)
from maxtic_next.constraints.adapters.ranger_dtl import (
    RangerDTLAdapter,
    convert_from_ranger_dtl,
)
from maxtic_next.constraints.adapters.eccetera import (
    EcceTERAAdapter,
    convert_from_eccetera,
)
from maxtic_next.constraints.adapters.artra import (
    ARTraAdapter,
    convert_from_artra,
)
from maxtic_next.constraints.adapters.alerax import (
    AleRaxAdapter,
    convert_from_alerax,
)
from maxtic_next.constraints.adapters._diagnostics import (
    CONVENTION_DONOR_ITSELF,
    CONVENTION_PARENT_OF_DONOR,
    AdapterDiagnosticError,
    LabelMismatchError,
    UpstreamParseError,
)
from maxtic_next.constraints.adapters.registry import (
    AdapterRegistry,
    AdapterEntry,
    ENDPOINT_CONVENTIONS,
    MixedFormatError,
    REGISTRY,
    detect_format,
    detect_formats,
    endpoint_convention,
    merge_constraint_sets,
    convert,
    convert_auto,
)

__all__ = [
    # 适配器类
    "ALEAdapter",
    "RangerDTLAdapter",
    "EcceTERAAdapter",
    "ARTraAdapter",
    "AleRaxAdapter",
    # 便捷转换函数
    "convert_from_ale",
    "convert_from_ranger_dtl",
    "convert_from_eccetera",
    "convert_from_artra",
    "convert_from_alerax",
    # 口径常量与诊断异常
    "ALE_DEFAULT_SOURCE",
    "CONVENTION_DONOR_ITSELF",
    "CONVENTION_PARENT_OF_DONOR",
    "AdapterDiagnosticError",
    "LabelMismatchError",
    "UpstreamParseError",
    # 注册表与自动检测
    "AdapterRegistry",
    "AdapterEntry",
    "ENDPOINT_CONVENTIONS",
    "MixedFormatError",
    "REGISTRY",
    "detect_format",
    "detect_formats",
    "endpoint_convention",
    "merge_constraint_sets",
    "convert",
    "convert_auto",
]
