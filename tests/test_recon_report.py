"""共享调和报告解析核心 ``_recon_report`` 单元测试。

覆盖 RANGER-DTLx / ARTra 同形调和报告的解析：
1. ``parse_transfer_line`` 对 Transfer / Replacing / Additive、含 Edge/Parent 字段、
   非转移事件的抽取；
2. ``parse_reconciliation_report`` 的块计数、边键聚合、家族规模、valid_labels 过滤。
"""

from maxtic_next.constraints.adapters._recon_report import (
    parse_transfer_line,
    parse_reconciliation_report,
)


class TestParseTransferLine:
    """``parse_transfer_line`` 单行抽取测试。"""

    def test_plain_transfer_with_edge_parent(self):
        line = "m3 = LCA[a, b]: Transfer, Mapping --> 61, Edge, Parent = 65, Recipient --> 62"
        assert parse_transfer_line(line) == ("61", "62")

    def test_replacing_transfer(self):
        line = "m3 = LCA[a, b]: Replacing Transfer, Mapping --> 61, Recipient --> 62"
        assert parse_transfer_line(line) == ("61", "62")

    def test_additive_transfer(self):
        line = "m5 = LCA[a, b]: Additive Transfer, Mapping --> 62, Recipient --> 65"
        assert parse_transfer_line(line) == ("62", "65")

    def test_speciation_is_not_transfer(self):
        assert parse_transfer_line("m1 = LCA[a,b]: Speciation, Mapping --> 59") == ("", "")

    def test_duplication_is_not_transfer(self):
        assert parse_transfer_line("m2 = LCA[a,b]: Duplication, Mapping --> 61") == ("", "")

    def test_leaf_is_not_transfer(self):
        assert parse_transfer_line("g1_0: Leaf Node") == ("", "")

    def test_empty_line(self):
        assert parse_transfer_line("") == ("", "")


class TestParseReconciliationReport:
    """``parse_reconciliation_report`` 整块解析测试。"""

    VALID = {"59", "61", "62", "65", "69"}

    def test_single_block_counts_and_aggregates(self):
        lines = [
            "Reconciliation for Gene Tree 1:",
            "g1_0: Leaf Node",
            "g2_1: Leaf Node",
            "m3 = LCA[a,b]: Transfer, Mapping --> 61, Recipient --> 62",
            "m5 = LCA[c,d]: Transfer, Mapping --> 62, Recipient --> 65",
        ]
        transfers, number, fam = parse_reconciliation_report(lines, self.VALID)
        assert transfers == {"61,62": 1, "62,65": 1}
        assert number == 1
        assert fam == 2  # 两个 Leaf Node

    def test_two_blocks_sum_counts(self):
        lines = [
            "Reconciliation for Gene Tree 1:",
            "m3 = LCA[a,b]: Transfer, Mapping --> 61, Recipient --> 62",
            "Reconciliation for Gene Tree 2:",
            "m3 = LCA[a,b]: Transfer, Mapping --> 61, Recipient --> 62",
        ]
        transfers, number, _ = parse_reconciliation_report(lines, self.VALID)
        assert transfers == {"61,62": 2}
        assert number == 2  # 两个块 => support = 2/2 = 1.0

    def test_invalid_labels_filtered(self):
        lines = [
            "Reconciliation for Gene Tree 1:",
            "m3 = LCA[a,b]: Transfer, Mapping --> ZZZ, Recipient --> 62",
            "m4 = LCA[a,b]: Transfer, Mapping --> 61, Recipient --> QQQ",
            "m5 = LCA[a,b]: Transfer, Mapping --> 61, Recipient --> 62",
        ]
        transfers, _, _ = parse_reconciliation_report(lines, self.VALID)
        # 只有两端点都合法的 61->62 被保留
        assert transfers == {"61,62": 1}

    def test_no_block_header_still_counts_one(self):
        lines = ["m3 = LCA[a,b]: Transfer, Mapping --> 61, Recipient --> 62"]
        transfers, number, _ = parse_reconciliation_report(lines, self.VALID)
        assert transfers == {"61,62": 1}
        assert number == 1

    def test_empty_input(self):
        transfers, number, fam = parse_reconciliation_report([], self.VALID)
        assert transfers == {}
        assert number == 0
        assert fam == 0
