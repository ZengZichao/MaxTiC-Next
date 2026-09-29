"""透明读取 gzip / tar.gz 输入。

上游工具的官方示例**全部**以压缩形式分发（MaxTiC 的 ``examples/reconciliations.tgz``
与 ``gene_trees.tgz``、Ranger-DTL / AleRax 输出常被打成 ``.tar.gz``，单文件示例常见
``xxx.gz``）。本模块是**唯一**的"文本输入读取"底层：:func:`open_text` /
:func:`read_text` / :func:`read_text_lines` / :func:`iter_text_chunks`，由

* :mod:`maxtic_next.io.parsing`（约束文件与 Newick 物种树），
* :mod:`maxtic_next.constraints.parsers`（便捷封装 ``parse_constraints_file``），
* 各上游适配器（``constraints/adapters/``：ALE / Ranger-DTL / ecceTERA / ARTra /
  AleRax 的逐文件解析）与 ``registry`` 的格式嗅探，
* :mod:`maxtic_next.dry_run` 的预检

共同调用，因此 ``.gz`` 输入在**任何**入口都不需要先手工解压。

判定规则（两条都认，避免"只看扩展名"或"只看魔数"各自的失效场景）
------------------------------------------------------------------------------
1. **内容魔数**：文件头两字节为 ``1f 8b`` 即按 gzip 解码 —— 于是被改名为
   ``constraints.txt`` 的 gzip 文件仍能读；
2. **扩展名**：``.gz`` / ``.tgz`` / ``.tar.gz`` 作为提示；扩展名声称是 gzip 而内容
   **不是** gzip 时按纯文本读，不再抛 ``BadGzipFile`` —— 误命名的**未压缩**文件
   因此仍然可用。

tar 归档（``.tar.gz`` / ``.tgz`` / ``.tar``，含改名的）同样按内容判定：

* 归档内**恰有一个**文件成员：直接读出该成员内容（透明，等同普通文本文件）；
* 归档内有**多个**成员（例如 ALE 官方 ``reconciliations.tgz`` 的 1000 个
  ``*.uml_rec``）：抛 :class:`CompressedArchiveError`（``ValueError`` 子类），
  消息里给出成员数、前若干成员名，以及**可直接执行**的
  ``tar -xzf 归档 -C 目录`` / ``tar -xzOf 归档 成员 > 文件`` 命令。
  多成员归档**不会**被静默拼接：归档内每个成员是一份独立的上游输出，
  拼接会被当成一个基因家族而按 ``--ale-min-family-size`` / ``--ale-min-support``
  之类的口径整体误算。

本模块只依赖标准库（``gzip`` / ``tarfile``），不引入任何第三方解压依赖。
"""

from __future__ import annotations

import atexit
import contextlib
import gzip
import io
import os
import shutil
import tarfile
import tempfile
from typing import Dict, Iterator, List, Optional, Sequence, Tuple

__all__ = [
    "GZIP_MAGIC",
    "CompressedArchiveError",
    "is_gzip_content",
    "is_gzip_path",
    "logical_path",
    "logical_stem",
    "tarfile_like",
    "is_tar_archive",
    "archive_member_names",
    "extract_archive",
    "expand_archive_inputs",
    "open_text",
    "read_text",
    "read_text_lines",
    "iter_text_chunks",
]

#: gzip 魔数（RFC 1952：两字节 ``1f 8b``）
GZIP_MAGIC = b"\x1f\x8b"

#: tar 归档在偏移 257 处的 ``ustar`` 魔数（POSIX.1-1988 起所有变体都写它）
_TAR_MAGIC = b"ustar"
_TAR_MAGIC_OFFSET = 257

#: 以扩展名识别的压缩/归档后缀
_GZIP_SUFFIXES = (".gz", ".tgz")
_TAR_SUFFIXES = (".tar", ".tar.gz", ".tgz")

#: 多成员归档的报错消息里最多列出多少个成员名
_LISTED_MEMBERS = 12

#: tar 探测时读取的头部字节数（``ustar`` 位于偏移 257，一个 512 字节块足够）
_HEAD_BYTES = 512


class CompressedArchiveError(ValueError):
    """多成员 tar 归档：给出可执行的手工解包指引，而不是裸 traceback。

    继承 :class:`ValueError`，因此 CLI 与 ``api.rank`` 既有的错误通道（退出码 4、
    ``错误：`` 前缀）原样适用，测试也可按 ``ValueError`` 断言。
    """


# ----------------------------------------------------------------------
# 判定
# ----------------------------------------------------------------------
def _path_of(path) -> str:
    """把 ``str`` / :class:`os.PathLike` 统一为字符串路径。"""
    return os.fspath(path) if hasattr(path, "__fspath__") else str(path)


def _head(path: str, n: int = _HEAD_BYTES) -> bytes:
    """读取原始文件头 ``n`` 字节；不可读时返回空串（由上层走"普通打开"）。"""
    try:
        with open(path, "rb") as fh:
            return fh.read(n)
    except OSError:
        return b""


def is_gzip_path(path) -> bool:
    """扩展名是否声称 gzip（``.gz`` / ``.tgz``）。**不**检查内容。"""
    return _path_of(path).lower().endswith(_GZIP_SUFFIXES)


def is_gzip_content(path) -> bool:
    """内容是否为 gzip（魔数 ``1f 8b``），与扩展名无关。"""
    return _head(path, 2) == GZIP_MAGIC


def logical_path(path) -> str:
    """去掉压缩/归档后缀后的"逻辑路径"，供按扩展名嗅探格式的上层复用。

    ``x.uml_rec.gz`` -> ``x.uml_rec``；``x.tar.gz`` / ``x.tgz`` -> ``x.tar``；
    未压缩输入原样返回。
    """
    p = _path_of(path)
    low = p.lower()
    for suffix in (".tar.gz", ".tgz"):
        if low.endswith(suffix):
            return p[: len(p) - len(suffix)] + ".tar"
    if low.endswith(".gz"):
        return p[: len(p) - len(".gz")]
    return p


def logical_stem(path) -> str:
    """文件名的"逻辑词干"：先剥压缩后缀，再剥最后一段扩展名。

    ``/x/gene_1.modif.ale.uml_rec.gz`` -> ``gene_1.modif.ale`` —— 与未压缩输入
    ``gene_1.modif.ale.uml_rec`` 得到**同一个**家族名，故 ``.gz`` 输入不会改变
    适配器写进 ``metadata["family"]`` 的家族标识。
    """
    name = os.path.basename(logical_path(path))
    return os.path.splitext(name)[0]


def tarfile_like(path) -> bool:
    """是否**可能**是 tar 归档（扩展名声明，或流内偏移 257 处有 ``ustar``）。

    只做廉价预判；真正的判定由 :func:`_read_tar_text` 内的 ``tarfile.open`` 完成。
    """
    p = _path_of(path)
    if p.lower().endswith(_TAR_SUFFIXES):
        return True
    if is_gzip_content(p):
        try:
            with gzip.open(p, "rb") as fh:
                head = fh.read(_HEAD_BYTES)
        except OSError:  # 声称 gzip 却读不出：交给普通文本路径去报清晰的错
            return False
    else:
        head = _head(p)
    return head[_TAR_MAGIC_OFFSET : _TAR_MAGIC_OFFSET + len(_TAR_MAGIC)] == _TAR_MAGIC


def archive_member_names(path) -> List[str]:
    """列出归档内的文件成员名（供报错消息与"先看看里面有什么"使用）。"""
    with _tar_context(_path_of(path)) as tf:
        if tf is None:
            return []
        return [m.name for m in tf.getmembers() if m.isfile()]


def is_tar_archive(path) -> bool:
    """该路径（含 ``.tgz`` / ``.tar.gz`` / 改名者）**是否**是 tar 归档。

    与 :func:`tarfile_like` 的区别：这里做真正的 ``tarfile.open`` 判定，
    用于"输入是不是归档"的分支决策（官方 ``.tgz`` 示例）。
    """
    with _tar_context(_path_of(path)) as tf:
        return tf is not None


def extract_archive(path: str, dest_dir: str) -> List[str]:
    """把 tar 归档（含 gzip 包裹）**安全地**解压到 ``dest_dir``，返回解出的文件路径。

    安全边界（等价 ``tarfile`` 的 ``filter="data"``，但兼容 ``requires-python >= 3.9``）：

    * 只解**常规文件**成员；目录 / 符号链接 / 硬链接 / 字符与设备 / FIFO 一律跳过；
    * 绝对路径成员与含 ``..`` 的成员直接跳过（防止解压逃逸出目标目录）；
    * 目录成员（``is_dir()``）会连带创建其内的常规文件成员，故不单独建目录树；
      需要保留目录结构时由成员名内的相对路径决定，同样受上一条约束。

    Args:
        path: 归档路径。
        dest_dir: 目标目录（须已存在）。

    Returns:
        解压得到的**文件**绝对/相对 ``dest_dir`` 的路径列表（按成员顺序）。
    """
    out: List[str] = []
    with _tar_context(_path_of(path)) as tf:
        if tf is None:
            return out
        for member in tf.getmembers():
            if not member.isfile():
                continue  # 目录/链接/设备/FIFO 一律不解
            name = member.name.replace("\\", "/")
            parts = [seg for seg in name.split("/") if seg not in ("", ".")]
            if name.startswith("/") or ".." in parts:
                continue  # 防"归档逃逸"
            try:
                extracted = tf.extractfile(member)
            except (tarfile.TarError, OSError):
                continue
            if extracted is None:
                continue
            target = os.path.join(dest_dir, *parts)
            parent = os.path.dirname(target)
            if parent:
                os.makedirs(parent, exist_ok=True)
            with extracted, open(target, "wb") as fh:
                fh.write(extracted.read())
            out.append(target)
    return out


def _archive_hint(path: str, members: Sequence[str]) -> str:
    """多成员（或空）归档的可执行指引文本（中文，与全项目错误信息同风格）。"""
    shown = ", ".join(members[:_LISTED_MEMBERS])
    more = (
        ""
        if len(members) <= _LISTED_MEMBERS
        else f" …（其余 {len(members) - _LISTED_MEMBERS} 个略）"
    )
    first = members[0] if members else "成员名"
    return (
        f"输入文件是一个含 {len(members)} 个文件成员的 tar 归档（gzip）：{path}\n"
        f"  成员：{shown}{more}\n"
        "  说明：MaxTiC-Next 透明解压**单文件** gzip（.gz / 无扩展名的 gzip 流）与"
        "**恰好一个成员**的 tar 归档；多成员归档不会被自动拼接成一个输入，因为"
        "归档内每个成员是一份独立的上游输出（各自的基因家族/样本），拼接会被当成"
        "一个家族而误算最小家族规模、最小支持度等按家族统计的阈值。\n"
        "  请任选其一，把成员变成真正的输入：\n"
        f"    1) 解包整个归档，再把解出的文件（或目录，AleRax 支持目录）作为输入：\n"
        f"       tar -xzf {path} -C ./unpacked && ls ./unpacked\n"
        f"       maxtic-next species_tree ./unpacked/* [--from auto]\n"
        f"    2) 只取其中一个成员：\n"
        f"       tar -xzOf {path} {first} > {os.path.basename(first)}\n"
        "  （Windows 10+ 自带 tar；--dry-run 会原样复述本提示）\n"
        "  详见用户手册 05_io_formats.md「压缩输入」与 09_faq_troubleshooting.md"
        "「上游示例是 .tgz，怎么直接喂给 MaxTiC-Next？」"
    )


@contextlib.contextmanager
def _tar_context(path: str):
    """打开（可能被 gzip 包裹的）tar 归档；不是 tar 时产出 ``None``。

    非 tar / 不可读一律回 ``None``，由调用方继续按 gzip 或纯文本处理。
    """
    stream = None
    tf = None
    try:
        stream = gzip.open(path, "rb") if is_gzip_content(path) else open(path, "rb")
        tf = tarfile.open(fileobj=stream, mode="r:")
        yield tf
    except (tarfile.TarError, OSError, EOFError):
        yield None
    finally:
        for handle in (tf, stream):
            if handle is None:
                continue
            try:
                handle.close()
            except Exception:  # noqa: BLE001 - 关闭失败不得掩盖真正的读取错误
                pass


def _decode(data: bytes, encoding: str, path: str) -> str:
    """按 ``encoding`` 解码归档成员，失败时给出带文件名与编码的清晰错误。"""
    try:
        return data.decode(encoding)
    except (UnicodeDecodeError, LookupError):
        raise ValueError(
            f"输入文件 {path} 的 tar 成员无法按 {encoding!r} 解码。"
            "  说明：上游工具输出应为 UTF-8 文本；请确认没有把二进制文件当作"
            "文本输入传入。"
        ) from None


def _normalize_newlines(text: str) -> str:
    """统一换行为 ``\\n``，与 ``open(..., newline=None)`` 的通用换行语义一致。"""
    if "\r" in text:
        text = text.replace("\r\n", "\n").replace("\r", "\n")
    return text


def _read_tar_text(path: str, encoding: str) -> Tuple[Optional[str], bool]:
    """读取 tar 归档的文本内容。

    Returns:
        ``(text, handled)``：``handled`` 为 ``False`` 表示这并非 tar 归档（调用方
        继续按 gzip / 纯文本读取）；为 ``True`` 时 ``text`` 是单成员归档的内容
        （多成员或空归档抛 :class:`CompressedArchiveError`）。
    """
    with _tar_context(path) as tf:
        if tf is None:
            return None, False
        members = [m for m in tf.getmembers() if m.isfile()]
        if len(members) != 1:
            raise CompressedArchiveError(_archive_hint(path, [m.name for m in members]))
        with tf.extractfile(members[0]) as fh:
            data = fh.read()
    return _normalize_newlines(_decode(data, encoding, path)), True


def _reader(path, encoding: str, errors: Optional[str]):
    """产出一个定位在**解压后**内容上的文本可读对象（调用方负责关闭）。"""
    p = _path_of(path)
    if tarfile_like(p):
        text, handled = _read_tar_text(p, encoding)
        if handled:
            return io.StringIO(text)
    # 魔数优先：内容以 1f 8b 开头才按 gzip 解码。扩展名声称 .gz 而内容不是 gzip
    # 的（误命名的未压缩文件）按纯文本读，因此仍然可用。
    if is_gzip_content(p):
        return gzip.open(p, "rt", encoding=encoding, errors=errors)
    return open(p, "r", encoding=encoding, errors=errors)


# ----------------------------------------------------------------------
# 归档输入 -> 成员输入（api / CLI / 适配器入口使用）
# ----------------------------------------------------------------------
#: 本进程为"归档输入"创建的临时解包目录（进程退出时回收，避免在 /tmp 里累积）
_TEMP_EXTRACT_DIRS: List[str] = []
_TEMP_CLEANUP_REGISTERED = False


def _cleanup_temp_extract_dirs() -> None:
    """回收 :func:`expand_archive_inputs` 创建的临时解包目录。"""
    while _TEMP_EXTRACT_DIRS:
        shutil.rmtree(_TEMP_EXTRACT_DIRS.pop(), ignore_errors=True)


def _register_temp_cleanup() -> None:
    global _TEMP_CLEANUP_REGISTERED
    if _TEMP_CLEANUP_REGISTERED:
        return
    atexit.register(_cleanup_temp_extract_dirs)
    _TEMP_CLEANUP_REGISTERED = True


def _archive_stem(path: str) -> str:
    """归档名（去掉压缩/归档后缀）作为解包子目录名，避免不同归档互相覆盖。"""
    stem = os.path.basename(logical_path(path))
    for suffix in (".tar", ".TAR"):
        if stem.endswith(suffix):
            stem = stem[: -len(suffix)]
    return stem or "members"


def expand_archive_inputs(
    paths: Sequence[str], tool: Optional[str] = None
) -> Tuple[List[str], Dict[str, List[str]]]:
    """把输入里的 **tar 归档**（``.tgz`` / ``.tar.gz`` / ``.tar``）展开为其成员。

    上游工具的官方示例常以归档分发（MaxTiC 的 ``examples/reconciliations.tgz`` 里是
    1000 个 ``*.uml_rec``；AleRax / Ranger 的结果包同理），
    "用户必须先手工解包、而且没有任何提示"。``api.rank`` / CLI / ``registry.convert``
    先调用本函数把归档**安全解压**到进程临时目录，再把成员当作输入，于是
    ``maxtic-next species.tree examples/reconciliations.tgz --from ale`` 可直接运行。

    ``.gz`` **单文件**不在此处理：那由各读取器（func:`open_text`）透明解压。

    Args:
        paths: 用户给出的输入路径列表。
        tool: 已选定的上游工具名。``"alerax"`` / ``"auto"`` 且归档是一棵 AleRax 输出树
            （成员路径含 ``reconciliations/summaries``）时，返回**解包根目录**，
            由 AleRax 适配器按自己的目录约定去取转移频率文件。

    Returns:
        ``(展开后的输入列表, {归档路径: 该归档展开得到的输入列表})``。第二个返回值供
        调用方打印提示并写入 ``run_metadata``（可追溯性：成员活在临时目录里）。
        非归档路径原样保留；空归档（无可解的常规文件）保留原路径，由读取层报错。
    """
    out: List[str] = []
    archives: Dict[str, List[str]] = {}
    for raw in paths:
        p = _path_of(raw)
        if not is_tar_archive(p):
            out.append(p)
            continue
        _register_temp_cleanup()
        dest_root = tempfile.mkdtemp(prefix="maxtic-input-archive-")
        _TEMP_EXTRACT_DIRS.append(dest_root)
        root = os.path.join(dest_root, _archive_stem(p))
        os.makedirs(root, exist_ok=True)
        members = extract_archive(p, root)
        if not members:
            out.append(p)  # 空归档：交给读取层给出清晰报错，不静默丢输入
            continue
        if tool in ("alerax", "auto") and any(
            "reconciliations/summaries" in m.replace(os.sep, "/") for m in members
        ):
            expanded: List[str] = [root]
        else:
            expanded = members
        archives[p] = expanded
        out.extend(expanded)
    return out, archives


@contextlib.contextmanager
def open_text(path, encoding: str = "utf-8-sig", errors: Optional[str] = None):
    """以文本模式打开 ``path``，透明处理 gzip 与单成员 tar(.gz) 归档。

    与内置 :func:`open` 的差别只有三处：``.gz`` / gzip 流自动解压；单成员 tar 归档
    自动取其成员；多成员归档抛 :class:`CompressedArchiveError`（含
    ``tar -xzOf`` 指引）。换行是通用模式（``\\r\\n`` / ``\\r`` 归一为 ``\\n``），
    与 ``open(...)`` 默认一致。

    Args:
        path: 输入路径（``str`` 或 :class:`os.PathLike`）。
        encoding: 文本编码；默认 ``utf-8-sig``（带 BOM 的文件不污染首列）。
        errors: codec 错误处理（``strict`` / ``replace`` …），默认 ``None``。

    Yields:
        文本可读对象（``read`` / ``readline`` / ``readlines`` 与迭代协议皆可用）。
    """
    handle = _reader(path, encoding, errors)
    try:
        yield handle
    finally:
        try:
            handle.close()
        except Exception:  # noqa: BLE001 - 关闭失败不得掩盖真正的读取错误
            pass


def read_text_lines(path, encoding: str = "utf-8-sig", errors: Optional[str] = None) -> List[str]:
    """读取文本文件为**行列表**（保留行尾 ``\\n``），透明支持 gzip / 单成员归档。"""
    with open_text(path, encoding=encoding, errors=errors) as fh:
        return fh.readlines()


def read_text(path, encoding: str = "utf-8-sig", errors: Optional[str] = None) -> str:
    """读取整个文本文件为一个字符串，透明支持 gzip / 单成员归档。"""
    with open_text(path, encoding=encoding, errors=errors) as fh:
        return fh.read()


def iter_text_chunks(
    path, chunk_size: int, encoding: str = "utf-8", errors: Optional[str] = None
) -> Iterator[str]:
    """按 ``chunk_size`` 个**字符**产出解压后的文本块（供格式嗅探使用）。

    不可读 / 解码失败时不产出任何块（调用方按"未识别"处理）；但多成员归档的
    :class:`CompressedArchiveError` **照常抛出** —— 那是用户必须看到的可执行提示，
    不是"未识别"。
    """
    try:
        with open_text(path, encoding=encoding, errors=errors) as fh:
            while True:
                chunk = fh.read(chunk_size)
                if not chunk:
                    return
                yield chunk
    except (OSError, UnicodeError):
        return
