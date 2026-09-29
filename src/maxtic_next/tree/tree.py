"""自研轻量有根树模块。

本模块是 Python 2 原版 ``script_tree.py`` 的逐函数忠实移植（Python 3 化），
供 ``MaxTiC-Next`` 的排序算法使用。与原版一致：

* 节点以整数 ``id`` 为内部键；
* 内部节点的"节点标识"为 **bootstrap 标签字符串**（等价原版 ``getBootstrap``），
  叶子以 ``getName`` 返回的名为标识；
* 算法层（``ranking/``）使用这些字符串标识构建 ``edge`` / ``graph`` 字典。

移植要点：

* ``dict.keys()`` 视图均用 ``list()`` 包裹后再迭代 / 排序；
* 整数除法改用 ``//``（本模块中距离累加为 ``float``，影响不大，但保持类型清晰）。
"""

from typing import Dict, List

from maxtic_next.tree.node import TreeNode


def _isfloat(value: str) -> bool:
    """判断字符串是否为浮点数（等价原版 ``script_tree.isfloat``）。"""
    try:
        float(value)
        return True
    except ValueError:
        return False


def _strip_quotes(name: str) -> str:
    """去掉类群名外层的成对引号（``'sp-1'`` / ``"sp-1"`` -> ``sp-1``）。"""
    if len(name) >= 2 and name[0] == name[-1] and name[0] in ("'", '"'):
        return name[1:-1]
    return name


class Tree:
    """有根系统发生树的轻量表示。

    内部以 ``{id: TreeNode}`` 字典存储节点，并提供等价原版 ``script_tree`` 的接口。
    """

    def __init__(self) -> None:
        self.nodes: Dict[int, TreeNode] = {}
        self.root: int = -1
        self.label_to_id: Dict[str, int] = {}
        # 原版 ``script_tree`` 残留字段：``@`` 标记的"祖先节点" id。
        # 本项目算法层**不读取**它（显式声明，避免动态属性）。
        self.ancestor: int = -1

    # ------------------------------------------------------------------
    # 解析：等价 script_tree.readTree
    # ------------------------------------------------------------------
    def read_newick(self, treeseq: str) -> "Tree":
        """解析 Newick 字符串，返回 ``self``。等价 ``script_tree.readTree``。

        Args:
            treeseq: Newick 字符串（含末尾 ``;``）。允许**折行书写**：
                若首行不含终止符 ``;``，则把后续各行拼接为一条序列后解析。

        Returns:
            当前 ``Tree`` 实例，便于链式调用。

        Raises:
            ValueError: 输入为空或畸形 Newick（括号不匹配、缺少右括号 / ``]``
                导致下标越界）时，给出清晰的位置说明，而不是裸 ``IndexError``。
        """
        # 原版只取首行；这里在首行缺少终止符 ';' 时拼接后续行（折行 Newick），
        # 否则维持"只读首行"的原版语义。
        if ";" not in treeseq.split("\n", 1)[0]:
            first_line = "".join(seg.strip() for seg in treeseq.split("\n")).strip()
        else:
            first_line = treeseq.split("\n", 1)[0].rstrip("\n")
        if not first_line.strip():
            raise ValueError("空的 Newick 输入：需要形如 '((A:1,B:1)1:1,C:1)2;' 的树字符串。")

        self.nodes = {}
        self.root = -1
        self.label_to_id = {}
        self.ancestor = -1

        id_node = 0
        nb_parenth = 0
        bppnumber = 0
        pile: List[int] = []
        t = 0
        n = len(first_line)

        while t < n:
            c = first_line[t]
            if c == "(":
                id_node += 1
                nb_parenth += 1
                nd = TreeNode(id=id_node)
                nd.name = "N" + str(id_node)
                nd.parent = -1
                nd.children = []
                nd.length = 0
                nd.isdup = ""
                nd.species = ""
                nd.bootstrap = ""
                nd.bppnumber = -1
                self.nodes[id_node] = nd
                if len(pile) > 0:
                    nd.parent = pile[-1]
                pile.append(id_node)
                t += 1
            elif c == ")":
                t += 1
                nb_parenth -= 1
                if not pile:
                    raise ValueError(f"畸形 Newick：右括号多于左括号（位置 {t}），括号不匹配。")
                self.nodes[pile[-1]].bppnumber = bppnumber
                bppnumber += 1
                if t < n and first_line[t] == "@":
                    t += 1
                    self.ancestor = pile[-1]
                while (
                    t < n
                    and first_line[t] != ":"
                    and first_line[t] != ";"
                    and first_line[t] != "["
                    and first_line[t] != ")"
                    and first_line[t] != ","
                ):
                    self.nodes[pile[-1]].bootstrap += first_line[t]
                    t += 1
                if t < n and first_line[t] == ":":
                    debut = t + 1
                    while (
                        t < n
                        and first_line[t] != ","
                        and first_line[t] != ")"
                        and first_line[t] != "["
                        and first_line[t] != ";"
                    ):
                        t += 1
                    longueur = float(first_line[debut:t])
                    self.nodes[pile[-1]].length = longueur
                    while (
                        t < n
                        and first_line[t] != ","
                        and first_line[t] != ")"
                        and first_line[t] != "["
                        and first_line[t] != ";"
                    ):
                        t += 1
                if t < n and first_line[t] == "[":
                    debut = t + 1
                    bracket = first_line[debut:].find("]")
                    if bracket < 0:
                        raise ValueError(f"畸形 Newick：缺少 ']'（位置 {debut}），NHX 注释未闭合。")
                    t = debut + bracket
                    chaine = first_line[debut:t]
                    mots = chaine.split(":")
                    for m in mots:
                        if m == "D=Y" or m == "D=T" or m == "Ev=GDup":
                            self.nodes[pile[-1]].isdup = "D"
                        if m[:2] == "S=":
                            self.nodes[pile[-1]].species = m[2:]
                        if m[:2] == "B=":
                            self.nodes[pile[-1]].bootstrap = m[2:]
                        if m[:3] == "ND=":
                            self.nodes[pile[-1]].nd = m[3:]
                        if _isfloat(m):
                            # 原版写 float：会让同一棵树的内部标签出现 str/float
                            # 混合类型（去重排序、字典键、报错信息全部受影响，
                            # ）。bootstrap 本质是**字符串标识**，此处保留
                            # 其文本形式。
                            self.nodes[pile[-1]].bootstrap = m
                    t += 1
                if t < n and first_line[t] == ":":
                    debut = t + 1
                    while (
                        t < n
                        and first_line[t] != ","
                        and first_line[t] != ")"
                        and first_line[t] != "["
                        and first_line[t] != ";"
                    ):
                        t += 1
                    longueur = float(first_line[debut:t])
                    self.nodes[pile[-1]].length = longueur
                    while (
                        t < n
                        and first_line[t] != ","
                        and first_line[t] != ")"
                        and first_line[t] != "["
                        and first_line[t] != ";"
                    ):
                        t += 1
                del pile[-1]
                if t < n and first_line[t] == ";":
                    t = n
            elif c == ";":
                t = n
            elif c == ",":
                t += 1
            elif c == " ":
                t += 1
            else:  # 叶子名
                id_node += 1
                nd = TreeNode(id=id_node)
                nd.parent = -1
                nd.children = []
                nd.length = 0
                nd.isdup = ""
                nd.species = ""
                nd.bootstrap = ""
                nd.bppnumber = bppnumber
                bppnumber += 1
                if len(pile) > 0:
                    nd.parent = pile[-1]
                self.nodes[id_node] = nd
                pile.append(id_node)
                debut = t
                while (
                    t < n
                    and first_line[t] != ","
                    and first_line[t] != ")"
                    and first_line[t] != ":"
                    and first_line[t] != ";"
                    and first_line[t] != "\n"
                    and first_line[t] != "["
                ):
                    t += 1
                nom = first_line[debut:t].strip()
                self.nodes[pile[-1]].name = nom
                if t < n and first_line[t] == ":":
                    debut = t + 1
                    while (
                        t < n
                        and first_line[t] != ","
                        and first_line[t] != ")"
                        and first_line[t] != "["
                        and first_line[t] != ";"
                    ):
                        t += 1
                    longueur = float(first_line[debut:t])
                    self.nodes[id_node].length = longueur
                if t < n and first_line[t] == "[":
                    debut = t + 1
                    bracket = first_line[debut:].find("]")
                    if bracket < 0:
                        raise ValueError(f"畸形 Newick：缺少 ']'（位置 {debut}），NHX 注释未闭合。")
                    t = debut + bracket
                    chaine = first_line[debut:t]
                    mots = chaine.split(":")
                    for m in mots:
                        if m[:2] == "S=":
                            self.nodes[pile[-1]].species = m[2:]
                        if m[:3] == "ND=":
                            self.nodes[pile[-1]].nd = m[3:]
                    t += 1
                if t < n and first_line[t] == ":":
                    debut = t + 1
                    while (
                        t < n
                        and first_line[t] != ","
                        and first_line[t] != ")"
                        and first_line[t] != "["
                        and first_line[t] != ";"
                    ):
                        t += 1
                    longueur = float(first_line[debut:t])
                    self.nodes[id_node].length = longueur
                del pile[-1]

        # 括号平衡校验：畸形输入（缺右括号等）在此被清晰拦截，避免后续断言崩溃
        if nb_parenth != 0:
            raise ValueError(f"畸形 Newick：括号不匹配（{nb_parenth} 个 '(' 未闭合），位置 {t}。")

        # 标签规范化：内部节点 bootstrap 同样 strip；类群名去掉成对引号。
        # 原版只对叶名 strip、内部标签不 strip，导致 " A" 与 "A" 成为两个不同端点。
        for nd in self.nodes.values():
            if isinstance(nd.bootstrap, str):
                nd.bootstrap = nd.bootstrap.strip()
            if isinstance(nd.name, str):
                nd.name = _strip_quotes(nd.name.strip())

        if not self.nodes:
            raise ValueError(
                f"Newick 未解析出任何节点（输入 {first_line[:60]!r}…）。"
                "请确认文件内容为一条有效的 Newick 序列。"
            )

        # 填充子节点列表（等价原版 readTree 末尾的 enfants 填充）
        for node in sorted(self.nodes.keys()):
            if not self.is_root(node):
                pere = self.get_parent(node)
                self.nodes[pere].children.append(node)

        # 计算根
        self.root = self.get_root()
        # 构建标签 -> id 映射（内部节点用 bootstrap，叶子用 name）
        self.label_to_id = {}
        for nid, nd in self.nodes.items():
            if self.is_leaf(nid):
                self.label_to_id[nd.name] = nid
            else:
                self.label_to_id[nd.bootstrap] = nid
        return self

    # ------------------------------------------------------------------
    # 结构自检辅助
    # ------------------------------------------------------------------
    def non_binary_internal_nodes(self) -> List[tuple]:
        """返回子节点数 != 2 的内部节点：``[(bootstrap 标签, 子节点数), ...]``。

        MaxTiC-Next **不支持**多歧（polytomy）物种树：``opt`` / ``mix`` /
        ``maximum_distance`` / ``random_order`` 只取前两个子节点，多余的子树会被
        静默丢弃。
        """
        out: List[tuple] = []
        for nid in self.get_nodes():
            if not self.is_leaf(nid):
                n = len(self.get_children(nid))
                if n != 2:
                    out.append((str(self.get_bootstrap(nid)), n))
        return out

    def duplicate_leaf_names(self) -> Dict[str, List[str]]:
        """返回重复出现的叶子名 -> 其父节点标签列表。"""
        seen: Dict[str, List[str]] = {}
        for nid in self.get_nodes():
            if self.is_leaf(nid) and not self.is_root(nid):
                name = self.get_name(nid)
                parent = self.get_bootstrap(self.get_parent(nid))
                seen.setdefault(name, []).append(str(parent))
        return {k: v for k, v in seen.items() if len(v) > 1}

    def copy(self) -> "Tree":
        """返回一棵**结构相同、枝长独立**的副本。

        排序流程会按"排序位次差"覆写枝长（``tree_from_order``），故必须在副本上
        进行，否则会污染调用方的树（``run()`` 非幂等）。
        """
        clone = Tree()
        clone.nodes = {}
        for nid, nd in self.nodes.items():
            nd2 = TreeNode(id=nd.id)
            nd2.name = nd.name
            nd2.parent = nd.parent
            nd2.children = list(nd.children)
            nd2.length = nd.length
            nd2.isdup = nd.isdup
            nd2.species = nd.species
            nd2.bootstrap = nd.bootstrap
            nd2.bppnumber = nd.bppnumber
            nd2.nd = nd.nd
            clone.nodes[nid] = nd2
        clone.root = self.root
        clone.label_to_id = dict(self.label_to_id)
        clone.ancestor = self.ancestor
        return clone

    # ------------------------------------------------------------------
    # 基础访问器
    # ------------------------------------------------------------------
    def get_nodes(self) -> List[int]:
        """等价 ``script_tree.getNodes``：返回所有节点 id（升序，等价于原版插入顺序）。"""
        return sorted(self.nodes.keys())

    def get_root(self) -> int:
        """等价 ``script_tree.getRoot``：自最小 id 节点向上回溯到根。"""
        keys = self.get_nodes()
        if not keys:
            raise ValueError("树为空（尚无任何节点），无法确定根节点。")
        start = keys[0]
        while not self.is_root(start):
            start = self.get_parent(start)
        return start

    def is_leaf(self, node: int) -> bool:
        """等价 ``script_tree.isLeaf``。"""
        return len(self.nodes[node].children) == 0

    def is_root(self, node: int) -> bool:
        """等价 ``script_tree.isRoot``。"""
        return self.nodes[node].parent == -1

    def get_parent(self, node: int) -> int:
        """等价 ``script_tree.getParent``。"""
        return self.nodes[node].parent

    def get_children(self, node: int) -> List[int]:
        """等价 ``script_tree.getChildren``。"""
        return self.nodes[node].children

    def get_name(self, node: int) -> str:
        """等价 ``script_tree.getName``。"""
        return self.nodes[node].name

    def get_bootstrap(self, node: int) -> str:
        """等价 ``script_tree.getBootstrap``（内部节点标签）。"""
        return self.nodes[node].bootstrap

    def set_length(self, node: int, length: float) -> None:
        """等价 ``script_tree.setLength``。"""
        self.nodes[node].length = length

    def get_length(self, node: int) -> float:
        """等价 ``script_tree.getLength``。"""
        return self.nodes[node].length

    def is_dup(self, node: int) -> bool:
        """等价 ``script_tree.isDup``。"""
        return self.nodes[node].isdup == "D"

    def internal_node_labels(self) -> List[str]:
        """返回所有内部节点的 bootstrap 标签列表（算法层节点标识集合）。"""
        return [self.get_bootstrap(n) for n in self.get_nodes() if not self.is_leaf(n)]

    def label_to_node_id(self, label: str) -> int:
        """将算法层字符串标识（bootstrap 标签或叶子名）映射回节点 id。"""
        return self.label_to_id[label]

    # ------------------------------------------------------------------
    # 拓扑查询
    # ------------------------------------------------------------------
    def is_ancestor(self, a: int, b: int) -> bool:
        """等价 ``script_tree.isAncestor``：判断 ``a`` 是否为 ``b`` 的祖先（含根特例）。"""
        if self.is_root(a):
            return True
        result = False
        current = b
        while (not result) and (not self.is_root(current)):
            if current == a:
                result = True
            else:
                current = self.get_parent(current)
        return result

    def get_leaves(self, a: int) -> List[int]:
        """等价 ``script_tree.getLeaves``：返回以 ``a`` 为根的子树的所有叶子 id。

        迭代实现（显式栈）：原版递归在极端不平衡（caterpillar）树上会
        ``RecursionError``，输出顺序与原递归实现一致。
        """
        if self.is_leaf(a):
            return [a]
        result: List[int] = []
        stack: List[int] = [a]
        while stack:
            node = stack.pop()
            if self.is_leaf(node):
                result.append(node)
                continue
            for child in reversed(self.get_children(node)):
                stack.append(child)
        return result

    def get_leaves_names(self) -> List[str]:
        """等价 ``script_tree.getLeavesNames``：返回所有叶子名。"""
        result: List[str] = []
        root = self.get_root()
        leaves = self.get_leaves(root)
        for leaf in leaves:
            result.append(self.get_name(leaf))
        return result

    def last_common_ancestor(self, a: int, b: int) -> int:
        """等价 ``script_tree.lastCommonAncestor``。"""
        ancestor = -1
        ancestorsa = [a]
        while not self.is_root(a):
            a = self.get_parent(a)
            ancestorsa.append(a)
        ancestorsb = [b]
        while not self.is_root(b):
            b = self.get_parent(b)
            ancestorsb.append(b)
        while len(ancestorsa) > 0 and len(ancestorsb) > 0 and ancestorsa[-1] == ancestorsb[-1]:
            ancestor = ancestorsa[-1]
            del ancestorsa[-1]
            del ancestorsb[-1]
        return ancestor

    def distance_from(self, a: int, b: int) -> float:
        """等价 ``script_tree.distanceFrom``：``a`` 与 ``b`` 在树上的距离。"""
        ancestor = self.last_common_ancestor(a, b)
        distance = 0.0
        while a != ancestor:
            distance += self.nodes[a].length
            a = self.get_parent(a)
        while b != ancestor:
            distance += self.nodes[b].length
            b = self.get_parent(b)
        return distance

    # ------------------------------------------------------------------
    # 输出：等价 script_tree.writeTree
    # ------------------------------------------------------------------
    def write_newick(self, NHX: bool = False) -> str:
        """等价 ``script_tree.writeTree(tree, root, NHX)``，返回 Newick 字符串。"""
        return self._write(self.get_root(), NHX)

    def _leaf_fragment(self, a: int, NHX: bool) -> str:
        """单个叶子节点的 Newick 片段（与原递归实现的叶子分支逐字符一致）。"""
        nd = self.nodes[a]
        if self.is_root(a):
            chaine = "(" + nd.name
        else:
            chaine = nd.name
        if nd.length != -1:
            chaine = chaine + ":" + str(nd.length)
        if NHX and nd.species != "":
            chaine = chaine + "[&&NHX:S=" + nd.species + "]"
        if self.is_root(a):
            chaine = chaine + ")" + str(nd.bootstrap)
        return chaine

    def _write(self, a: int, NHX: bool) -> str:
        """自 ``a`` 起写出 Newick。

        显式栈迭代实现（极端不平衡树上的递归爆栈），输出与原递归
        实现逐字符一致。
        """
        fragments: Dict[int, str] = {}
        stack: List[tuple] = [(a, False)]
        while stack:
            node, expanded = stack.pop()
            if self.is_leaf(node):
                fragments[node] = self._leaf_fragment(node, NHX)
                continue
            if not expanded:
                stack.append((node, True))
                for child in reversed(self.get_children(node)):
                    stack.append((child, False))
                continue
            nd = self.nodes[node]
            chaine = (
                "("
                + ",".join(fragments[c] for c in self.get_children(node))
                + ")"
                + str(nd.bootstrap)
            )
            if (not self.is_root(node)) and nd.length != -1:
                chaine = chaine + ":" + str(nd.length)
            if NHX and (nd.isdup != "" or nd.species != ""):
                chaine = chaine + "[&&NHX:"
                if nd.species != "":
                    chaine = chaine + "S=" + nd.species
                if nd.isdup == "D" or nd.isdup == "WGD":
                    chaine = chaine + ":D=Y"
                chaine = chaine + "]"
            if self.is_root(node):
                chaine = chaine + ";"
            fragments[node] = chaine
        return fragments[a]
