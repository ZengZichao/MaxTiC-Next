"""IO 层公开接口。"""

from maxtic_next.io.compression import CompressedArchiveError, open_text, read_text, read_text_lines
from maxtic_next.io.parsing import read_newick_file, read_constraints_file
from maxtic_next.io.output import write_three_files, format_summary

__all__ = [
    "read_newick_file",
    "read_constraints_file",
    "write_three_files",
    "format_summary",
    # ：压缩输入（.gz / 单成员 .tar.gz）透明读取
    "open_text",
    "read_text",
    "read_text_lines",
    "CompressedArchiveError",
]
