#!/usr/bin/env nextflow
// MaxTiC-Next 两阶段流程（P2-WF / T16）
//
//   Stage 1（from_ale=true）：把 ALE .uml_rec 文件经适配器转为约束文件；
//   Stage 2：对约束做 MaxTiC 排序，产出三输出文件 + 交互式 HTML 报告。
//
// 参数名与 `maxtic-next` CLI 一致（见 src/maxtic_next/cli.py）。
//
// 运行（在**仓库根目录**，与 docs/manual/{zh,en}/07 章一致）：
//   nextflow -C workflows/nextflow.config run workflows/main.nf \
//       --species_tree examples/minitree.tree \
//       --constraints examples/Cyano_CUTConstraints.tsv
// 或（ALE 模式；注意物种树需与 ALE 物种命名一致）：
//   nextflow -C workflows/nextflow.config run workflows/main.nf \
//       --species_tree examples/adapters/ale/species.tree --from_ale \
//       --ale_rec_files 'examples/adapters/ale/rec.uml_rec'
//
// 本文件相对旧版的修正（审阅 M18）：
//   1. DSL2 不允许 `process` 写在 `if` 块里（报 "process declaration appears
//      inside a control flow block"）—— 两个 process 现均声明在顶层，
//      条件分支只留在 workflow 里决定是否**调用** Stage 1；
//   2. 旧版 `GENERATE_CONSTRAINTS` 声明 `input: path species; each path(rec)`，
//      调用处却传 `(rec_ch, species_ch)`，实参与形参对调，任务一启动即失败。
//      现改为 `tuple path(species), path(rec)`，调用处按同序传参，
//      且每个 .uml_rec 各起一个任务（`each` 对多文件通道并不会扇出）；
//   3. 旧版 RANK 的 output 恒写长文件名，而 params.output_style 默认 legacy、
//      workflows/config.yaml 却写 short —— 一旦按 CLI 默认（short）运行，
//      声明的输出文件不存在，流程报 "output file not found"。
//      现输出名由 params.output_style **推导**，两侧不可能漂移；
//   4. Stage 2 加 `-f/--force`：Nextflow 在 --resume 命中缓存或任务目录复用时，
//      产物可能已存在，而 CLI 自审阅 M17 起默认以退出码 3 拒绝覆盖。
//      任务工作目录归流程管理器所有，此处覆盖是安全的。
//
// 说明：
//   * 物种树与约束文件经 channel 输入自动暂存（staging）到任务目录，
//     相对路径同样可用；
//   * 输出文件名以约束文件的完整文件名（含扩展名）为前缀。
//     注意：工具取**第一个**约束文件作前缀；适配器一次产出多个约束文件时，
//     本流程把每个约束文件作为独立 RANK 输入，故各有一次完整产物。

nextflow.enable.dsl = 2

params {
    species_tree          = "examples/minitree.tree"
    constraints           = "examples/Cyano_CUTConstraints.tsv"
    seed                  = 42
    from_ale              = false
    local_search          = 0
    min_transfer_distance = 0
    threshold_constraints = 0
    random_trees          = 0
    // 与 CLI 默认一致：short（.mt.*.tsv 后缀）。改为 legacy 则得到原版长名。
    output_style          = "short"
    outdir                = "maxtic_next_out"

    // ALE 适配器模式参数（仅 from_ale=true 时使用）
    ale_rec_files         = ['examples/adapters/ale/rec.uml_rec']
    ale_cache_dir         = "ale_cache"
    ale_source            = "trf"     // 'trf'（转移事件，默认）/ 'rec'（调和事件）
    ale_min_support       = 0.05
    ale_min_family_size   = 5
}

// ---------------------------------------------------------------------------
// 输出名按 params.output_style 推导 —— 与 src/maxtic_next/config.py 的
// OUTPUT_SUFFIXES 一一对应，避免声明与实际产物漂移（审阅 M18 第 3 点）。
// ---------------------------------------------------------------------------
def outputSuffixes(String style) {
    if (style == 'legacy') {
        return [
            informative: '_MT_output_filtered_list_of_weighted_informative_constraints',
            conflicting: '_MT_output_list_of_constraints_conflicting_with_best_order',
            partial    : '_MT_output_partial_order',
            random_dist: '_distribution_random',
        ]
    }
    return [
        informative: '.mt.informative.tsv',
        conflicting: '.mt.conflicts.tsv',
        partial    : '.mt.partial_order.tsv',
        random_dist: '.mt.random_dist.tsv',
    ]
}

// ---------------------------------------------------------------------------
// Stage 1：ALE .uml_rec -> 约束文件（仅在 workflow 中按 from_ale 决定是否调用）
// ---------------------------------------------------------------------------
process GENERATE_CONSTRAINTS {
    tag { rec }
    publishDir "${params.outdir}", mode: 'copy', overwrite: true

    input:
        tuple path(species), path(rec)

    output:
        path "${rec.baseName}_constraints.tsv", emit: constraints

    script:
    """
    mkdir -p "${params.ale_cache_dir}"
    maxtic-next "${species}" "${rec}" \
        --from ale \
        --ale-cache-dir "${params.ale_cache_dir}" \
        --ale-source "${params.ale_source}" \
        --ale-min-support ${params.ale_min_support} \
        --ale-min-family-size ${params.ale_min_family_size} \
        -o "${rec.baseName}_constraints.tsv"
    """
}

// ---------------------------------------------------------------------------
// Stage 2：约束 -> MaxTiC 排序（三输出文件 + HTML 报告）
// ---------------------------------------------------------------------------
process RANK {
    tag { constr }
    publishDir "${params.outdir}", mode: 'copy', overwrite: true

    // 工具以约束文件的完整文件名（含扩展名）作为输出前缀（getName 保留扩展名）
    input:
        tuple path(species), path(constr)

    output:
        path("${constr.getName()}${outputSuffixes(params.output_style).partial}")
        path("${constr.getName()}${outputSuffixes(params.output_style).informative}")
        path("${constr.getName()}${outputSuffixes(params.output_style).conflicting}")
        path("${constr.getName()}.html")
        // 第 4 个产物：仅在 --random-trees > 0 时存在（审阅 m13）
        path("${constr.getName()}${outputSuffixes(params.output_style).random_dist}"), optional: true

    script:
    """
    maxtic-next "${species}" "${constr.getName()}" \\
        --seed ${params.seed} \\
        --local-search ${params.local_search} \\
        --min-transfer-distance ${params.min_transfer_distance} \\
        --threshold-constraints ${params.threshold_constraints} \\
        --random-trees ${params.random_trees} \\
        --output-style ${params.output_style} \\
        --force
    """
}

// ---------------------------------------------------------------------------
// 工作流编排
// ---------------------------------------------------------------------------
workflow {
    // 物种树被 Stage 1 与 Stage 2 **共用**。`Channel.fromPath` 返回热通道，
    // 只能订阅一次，复用它会在 Stage 2 报 "channel already subscribed"；
    // 故先显式检查存在性，再用 `Channel.value`（冷通道、可多次订阅并自动广播）。
    def species_file = file(params.species_tree)
    if (!species_file.exists()) {
        error "物种树文件不存在：${params.species_tree}"
    }
    species_ch = Channel.value(species_file)

    if (params.from_ale) {
        Channel.fromPath(params.ale_rec_files, checkIfExists: true)
               .set { rec_ch }
        GENERATE_CONSTRAINTS(species_ch, rec_ch)
            .constraints.collect()
            .set { constr_ch }
    } else {
        Channel.fromPath(params.constraints, checkIfExists: true)
               .set { constr_ch }
    }

    RANK(species_ch, constr_ch)
}
