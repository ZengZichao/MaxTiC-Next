"""目标类群剪裁。

不可妥协约束：跨类群约束策略必须显式化，**禁止静默引入
虚假偏序**。默认策略（保守）：

* **仅保留两端点均在目标类群内的约束**（``--target-clade <node>`` 指定类群根）；
* 外部节点若为类群根**祖先**（即"更古老"），在显式开启 ``--target-clade-ancestor-map``
  时可映射为物种树根，并继承"外部祖先比目标类群根更古老"的关系：仅当该祖先落在
  donor（较早）侧时保留为 ``根 -> 类群内节点``；若落在 receptor（较晚）侧则丢弃
  （否则会断言"类群内节点早于根"，与根最古老的事实矛盾，属虚假偏序）；
* 旁支 / 更远分支（既不在类群内也非类群祖先）一律**丢弃**。

**标签对齐保证**：本模块仅过滤约束，或对端点字符串做"映射到根标签"的替换，**绝不**
修改树的拓扑或节点标识（bootstrap 标签）。因此剪裁前后内部节点标签完全一致，
不会因重根（ete3/networkx）改变任何标识。

假设（已写入 CLI help 与 README）：
* "目标类群"由单个内部节点标签界定，类群 = 以该节点为根的子树；
* 约束端点必须是物种树中存在的标签（内部节点 bootstrap 或叶子名），否则丢弃；
* 默认不引入任何跨类群偏序；祖先映射为需显式开启的、文档化的保守规则。
"""

from typing import Optional

from maxtic_next.constraints.constraint import Constraint, ConstraintSet
from maxtic_next.tree.tree import Tree


def _label_to_id(tree: Tree, label: str) -> Optional[int]:
    """按标签（内部节点 bootstrap 或叶子名）查找节点 id；不存在返回 None。"""
    try:
        return tree.label_to_node_id(label)
    except KeyError:
        return None


def _in_clade(tree: Tree, target_id: int, label: str) -> bool:
    """判断标签 ``label`` 是否位于目标类群内（target 的后裔，含自身）。"""
    nid = _label_to_id(tree, label)
    if nid is None:
        return False
    return tree.is_ancestor(target_id, nid)


def _is_clade_ancestor(tree: Tree, target_id: int, label: str) -> bool:
    """判断标签 ``label`` 是否为目标类群的祖先（在类群之外但更古老）。"""
    nid = _label_to_id(tree, label)
    if nid is None or nid == target_id:
        return False
    # nid 是 target 的祖先 => nid 在 target 之上（更古老）
    return tree.is_ancestor(nid, target_id)


def prune_constraints(
    tree: Tree, cset: ConstraintSet, target_clade: str, ancestor_map: bool = False
) -> ConstraintSet:
    """按目标类群剪裁约束集，返回新的 ``ConstraintSet``。

    Args:
        tree: 物种树（用于判断谱系关系）。
        cset: 原始约束集。
        target_clade: 目标类群根节点标签（内部节点 bootstrap）。
        ancestor_map: 是否开启"外部祖先映射到根"的保守规则（默认 False）。

    Returns:
        剪裁后的 ``ConstraintSet``。

    Raises:
        KeyError: 当 ``target_clade`` 不是物种树中的有效标签时。
    """
    target_id = _label_to_id(tree, target_clade)
    if target_id is None:
        raise KeyError(f"--target-clade 指定的节点 {target_clade!r} 不在物种树中")
    root_label = tree.get_bootstrap(tree.get_root())

    kept: ConstraintSet = ConstraintSet()
    for c in cset.constraints:
        d_in = _in_clade(tree, target_id, c.donor)
        r_in = _in_clade(tree, target_id, c.receptor)

        if d_in and r_in:
            # 两端点都在类群内：直接保留
            kept.add(c)
            continue

        if not ancestor_map:
            # 默认策略：丢弃任何跨类群约束（保守，禁止虚假偏序）
            continue

        # 祖先映射模式：仅当外部端点为类群祖先时映射到根
        nd, nr = c.donor, c.receptor
        d_is_clade_anc = (not d_in) and _is_clade_ancestor(tree, target_id, c.donor)
        r_is_clade_anc = (not r_in) and _is_clade_ancestor(tree, target_id, c.receptor)

        if d_is_clade_anc:
            nd = root_label
        if r_is_clade_anc:
            nr = root_label

        # 若两端都被映射到根（均为外部祖先）：无信息，丢弃
        if nd == root_label and nr == root_label:
            continue
        # 若 receptor 被映射到根而 donor 在类群内：会断言"类群内节点早于根"，丢弃
        if nr == root_label and nd != root_label:
            continue
        # 若 donor 被映射到根、receptor 在类群内：保留（根更古老），合法
        # 注意：必须附加 r_in 守卫——仅当 receptor 在类群内才保留为
        # ``根 -> 类群内节点``；否则（receptor 是旁支/更远分支，不在类群内、
        # 也非类群祖先）会错误地引入 ``根 -> 旁支`` 虚假偏序，必须丢弃
        # （见本模块 docstring 的「旁支一律丢弃」条目）。
        if nd == root_label and nr != root_label and r_in:
            kept.add(Constraint(donor=nd, receptor=nr, weight=c.weight, metadata=dict(c.metadata)))
            continue
        # 其余（含旁支/更远分支、无法映射的一端）：丢弃
        continue

    return kept
