#!/usr/bin/env bash
# MaxTiC-Next 多软件兼容适配框架 —— 各上游工具样例演示脚本
#
# 演示用 5 个上游工具的样例输出（真实格式）驱动 MaxTiC-Next 排序，并展示自动检测。
# 物种树使用 ../minitree.tree（内部节点标签为 bootstrap 数字：42/59/61/62/65/67/69 等）。
#
# 运行（micromamba python-3.11 环境）：
#   micromamba run -n python-3.11 bash examples/adapters/run_adapters.sh
# 或已安装控制台命令后直接：
#   bash examples/adapters/run_adapters.sh
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SP="$HERE/../minitree.tree"

# 用 python -m 方式调用，免依赖控制台脚本；小样例关闭家族规模过滤与 HTML 报告
MT=(python -m maxtic_next --no-html --seed 42 --ale-min-family-size 0)

echo "==================== 1) RANGER-DTLx (--from ranger) ===================="
"${MT[@]}" "$SP" "$HERE/ranger/FAM1.dtl" --from ranger

echo "==================== 2) ecceTERA (--from eccetera, recPhyloXML) ========"
"${MT[@]}" "$SP" "$HERE/eccetera/FAM1.recphyloxml" --from eccetera

echo "==================== 3) ARTra (--from artra, all transfers) ============"
"${MT[@]}" "$SP" "$HERE/artra/FAM1.txt" --from artra

echo "==================== 3b) ARTra (仅 replacing transfer) ================="
"${MT[@]}" "$SP" "$HERE/artra/FAM1.txt" --from artra --artra-transfer-kind replacing

echo "==================== 4) AleRax (--from alerax, 目录自动展开) ==========="
"${MT[@]}" "$SP" "$HERE/alerax/run" --from alerax

echo "==================== 5) 自动检测 (--from-auto) ========================="
"${MT[@]}" "$SP" "$HERE/eccetera/FAM1.recphyloxml" --from-auto

echo "==================== 6) 两阶段：仅生成统一约束文件 ====================="
"${MT[@]}" "$SP" "$HERE/ranger/FAM1.dtl" --from ranger -o "$HERE/ranger_constraints.tsv"
echo "已写出统一约束：$HERE/ranger_constraints.tsv"
head -5 "$HERE/ranger_constraints.tsv"

echo "全部样例运行完成。"
