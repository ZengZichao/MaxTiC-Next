# MaxTiC-Next —— MaxTiC 的 Python 3 重写镜像
#
# 基于 python:3.11-slim，安装本包及运行所需的第三方库（Jinja2 / Plotly）。
# 入口为 `maxtic-next` 命令行；CLI 的退出码即容器退出码（错误时非零，退出码清晰传播）。
#
# 镜像命名：统一为 **小写** `maxtic-next`、标签 `0.1.0`，与
# `Singularity.def` 的产物名、`workflows/nextflow.config` 的 `docker.image` 以及
# `docs/manual/{zh,en}/01_installation.md`、`07_workflows_deployment.md` 中的示例
# 逐字一致。Docker 仓库名必须全小写，大写形式会报 `invalid reference format`。
#
# 构建（在仓库根目录；`docs/`、`tests/`、`.git` 等已由 .dockerignore 排除）：
#   docker build -t maxtic-next:0.1.0 .
#
# 运行（基本用法）：
#   docker run --rm -v "$PWD/out":/data maxtic-next:0.1.0 \
#       /data/minitree.tree /data/Cyano_CUTConstraints.tsv --seed 42
#   （先把 examples/ 里的两个输入拷进 ./out，或直接把宿主 examples 挂到 /data：
#     docker run --rm -v "$PWD/examples":/data maxtic-next:0.1.0 \
#         /data/minitree.tree /data/Cyano_CUTConstraints.tsv --seed 42）
#
# 开启 MCMC 采样（初步实现，收敛诊断未经验证；样本不得当作后验样本使用）：
#   docker run --rm -v "$PWD/examples":/data maxtic-next:0.1.0 \
#       /data/minitree.tree /data/Cyano_CUTConstraints.tsv --seed 42 \
#       --mcmc --mcmc-iters 2000
#
# 产物覆盖：把宿主目录挂载进容器时，若上一次运行已在该目录留下产物，
# CLI 会以退出码 3 拒绝覆盖；此时加 `-f/--force` 或换 `-p` 前缀。
#
# 退出码：maxtic-next 的 sys.exit 值直接成为容器退出码；`set -e` 非必需，因为
# ENTRYPOINT 直接执行 maxtic-next，其返回值即容器返回值。

# 版本与基础镜像参数化：
# tag 必须可校验——`docker inspect` 读得到 LABEL，
# 且 tests/test_architecture_gates.py 会断言它与 pyproject / __init__ / 配方一致。
ARG MAXTIC_VERSION=0.1.0
ARG PYTHON_IMAGE=python:3.11-slim

FROM ${PYTHON_IMAGE}

LABEL org.opencontainers.image.title="MaxTiC-Next" \
      org.opencontainers.image.version="${MAXTIC_VERSION}" \
      org.opencontainers.image.licenses="CECILL-2.0"

# 避免交互式提示与多余缓存
ENV PYTHONUNBUFFERED=1 \
    PYTHONDONTWRITEBYTECODE=1 \
    PIP_NO_CACHE_DIR=1

WORKDIR /app

# 复制项目（含 pyproject.toml / src / examples / LICENSE 等；
# 其余内容由 .dockerignore 排除，见该文件注释）
COPY . /app

# 安装本包（src 布局）及运行所需依赖
# - 基础算法仅依赖标准库；
# - HTML 报告需要 jinja2 + plotly（缺了会降级为基础模板，故一并装上）；
# - pip install . 已按 pyproject 装上核心依赖，此处显式补齐 report 所需的两个。
RUN python -m pip install --upgrade pip \
    && pip install . \
    && pip install jinja2 plotly

# 冒烟测试（默认注释掉以保持构建快速；取消注释即可在构建期验证入口点）
# RUN maxtic-next --version && maxtic-next --help

# 默认命令（可被 `docker run ...` 的命令行参数覆盖）
CMD ["examples/minitree.tree", "examples/Cyano_CUTConstraints.tsv", "--seed", "42"]

# 入口：CLI 退出码即容器退出码
ENTRYPOINT ["maxtic-next"]
