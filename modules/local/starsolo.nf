// modules/local/starsolo.nf
process STAR_SOLO {

    tag "$meta.id"
    errorStrategy 'ignore'

    publishDir "$params.outdir/starsolo/${meta.id}/logs", \
        mode: 'copy', \
        pattern: "${meta.id}.Log.final.out", \
        enabled: true

    publishDir "$params.outdir/starsolo/${meta.id}/bam", \
        mode: 'copy', \
        pattern: "${meta.id}.Aligned.sortedByCoord.out.bam", \
        enabled: true

    publishDir "$params.outdir/starsolo/${meta.id}", \
        mode: 'copy', \
        pattern: "${meta.id}.Solo.out", \
        enabled: false

    input:
    // for re-sequence samples 260301
    tuple val(meta), path(reads1), path(reads2), path(genomeDir)

    output:
    tuple val(meta), path("${meta.id}.Solo.out"), emit: solo_out_dir
    path "${meta.id}.Log.final.out", emit: log
    path "${meta.id}.Aligned.sortedByCoord.out.bam", emit: bam

    script:
    // 【核心修改】：判断如果是同名多行的文件（List），则用逗号合并；如果是单个文件，直接使用。
    def r1 = reads1 instanceof List ? reads1.join(',') : reads1
    def r2 = reads2 instanceof List ? reads2.join(',') : reads2

    def solo_features = params.run_mode.toLowerCase() == 'velocity' ? 'Gene GeneFull Velocyto' : 'GeneFull'

    """
    STAR \\
        --runThreadN ${task.cpus} \\
        --genomeDir ${genomeDir} \\
        --readFilesIn ${r2} ${r1} \\
        --readFilesCommand zcat \\
        --outFileNamePrefix ${meta.id}. \\
        --soloType ${params.soloType} \\
        --soloCBwhitelist ${params.soloCBwhitelist} \\
        --soloUMIfiltering ${params.soloUMIfiltering} \\
        --soloUMIdedup ${params.soloUMIdedup} \\
        --soloUMIlen ${params.soloUMIlen} \\
        --outSAMtype BAM SortedByCoordinate \\
        --clipAdapterType CellRanger4 \\
        --outFilterScoreMin 30 \\
        --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts \\
        --soloCellFilter EmptyDrops_CR 10000 0.99 10 45000 90000 500 0.01 20000 0.01 10000 \\
        --soloFeatures ${solo_features} \\
        --soloMultiMappers EM \\
        --outSAMattributes NH HI nM AS CR UR CB UB GX GN sS sQ sM \\
        --outFilterMultimapNmax 20 --outBAMsortingBinsN 200 \
        ${task.ext.args ?: ''}
    """
}
