# 01 · 安装

## 1.1 环境要求

| 项 | 要求 |
|----|------|
| Python | ≥ 3.9（**推荐 3.11**） |
| 核心算法依赖 | 仅标准库（`argparse` / `random` / `math` / `functools` 等），无需联网 |
| HTML 报告（可选） | `jinja2` ≥ 3.0、`plotly` ≥ 5.0 |
| 剪裁 / 图操作（可选） | 无第三方依赖：`prune.py` 仅用标准库，刻意不引入 ete3 / networkx |
| 开发/测试（可选） | `pytest` ≥ 8.0 |

> **重要**：核心排序功能**不依赖任何第三方库**。只有交互式 HTML 报告需要 `jinja2` + `plotly`。

## 1.2 方式一：pip 安装（推荐）

```bash
# 不必事先克隆：直接按仓库地址安装
pip install "git+https://github.com/ZengZichao/MaxTiC-Next.git"

# 需要交互式 HTML 报告（jinja2 + plotly）时，克隆后带 extras 安装
git clone https://github.com/ZengZichao/MaxTiC-Next.git
cd MaxTiC-Next
pip install -e ".[report]"
```

验证：

```bash
maxtic-next --help
python -c "import maxtic_next; print(maxtic_next.__version__)"   # 0.1.0
```

## 1.3 方式二：免安装（PYTHONPATH）

离线或不希望安装时，直接从源码运行：

```bash
git clone https://github.com/ZengZichao/MaxTiC-Next.git
cd MaxTiC-Next
export PYTHONPATH=src
python -m maxtic_next --help
```

## 1.4 方式三：micromamba / conda 环境（本项目推荐）

本项目所有脚本与测试均在 micromamba `python-3.11` 环境下运行、验证。

```bash
# 若尚无该环境
micromamba create -n python-3.11 python=3.11 -y
micromamba run -n python-3.11 pip install -e ".[report]"

# 之后所有命令都用 micromamba run 包裹：
micromamba run -n python-3.11 maxtic-next examples/minitree.tree \
    examples/Cyano_CUTConstraints.tsv --seed 42
```

## 1.5 方式四：容器（Docker / Singularity）

见 [07 · 流程封装与部署](07_workflows_deployment.md)。基础镜像 `python:3.11-slim`：

```bash
docker build -t maxtic-next:0.1.0 .
docker run --rm -v "$PWD/examples":/data maxtic-next:0.1.0 \
    /data/minitree.tree /data/Cyano_CUTConstraints.tsv --seed 42
```

> 镜像名必须全小写（`MaxTiC-Next` 会报 `invalid reference format`）；标签与
> `Dockerfile` / `Singularity.def` / `workflows/nextflow.config` 保持一致。
> 构建上下文由 `.dockerignore` 收敛（`.git`、`tests`、`docs` 等不入镜像）。

## 1.6 运行测试套件（可选）

```bash
python -m pytest tests/ -q
# 期望：全部通过（当前 486 passed, 1 skipped；覆盖适配器、参考实现差分、
#       边界、CLI 开关等；唯一的 skip 需要 python2 与原版 MaxTiC.py 才执行）
```

## 1.7 卸载

```bash
pip uninstall MaxTiC-Next
```

## 1.8 常见安装问题

| 现象 | 原因 / 解决 |
|------|-------------|
| `maxtic-next: command not found` | 未 `pip install`；改用 `PYTHONPATH=src python -m maxtic_next` |
| `ModuleNotFoundError: jinja2` | 未装 report 依赖；`pip install -e ".[report]"` 或加 `--no-html` |
| HTML 报告生成失败但排序正常 | 报告为可选特性，失败不阻断主流程（会打印警告） |
| 中文路径 / 空格路径报错 | 用引号包裹路径；建议在无空格目录下运行 |

> 上一章：[00 · 总览](00_index.md) ｜ 下一章：[02 · 快速开始](02_quickstart.md)
