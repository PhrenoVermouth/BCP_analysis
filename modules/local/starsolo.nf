// modules/local/starsolo.nf
process STAR_SOLO {
    tag "$meta.id"
    publishDir "$params.outdir/starsolo/${meta.id}", mode: 'copy'

    input:
    tuple val(meta), path(reads), path(genomeDir)

    output:
    tuple val(meta), path("${meta.id}.Solo.out")
    path "${meta.id}.Log.final.out", emit: log

    script:
    """
    STAR \\
        --runThreadN ${task.cpus} \\
        --genomeDir ${genomeDir} \\
        --readFilesIn ${reads[1]} ${reads[0]} \\
        --readFilesCommand zcat \\
        --outFileNamePrefix ${meta.id}. \\
        --soloType ${params.soloType} \\
        --soloCBwhitelist ${params.soloCBwhitelist} \\
        --soloUMIfiltering ${params.soloUMIfiltering} \\
        --soloUMIdedup ${params.soloUMIdedup} \\
        --soloUMIlen ${params.soloUMIlen} \\
        --outSAMtype BAM SortedByCoordinate \\ //@Prateek: parameters below are fixed based on Notion, might add more flexibility in the future
        --clipAdapterType CellRanger4 \\
        --outFilterScoreMin 30 \\ 
        --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts \\
        --soloCellFilter EmptyDrops_CR 10000 0.99 10 45000 90000 500 0.01 20000 0.01 10000 \\
        --soloFeatures GeneFull \\
        --soloMultiMappers EM \\
        --outSAMattributes NH HI nM AS CR UR CB UB GX GN sS sQ sM \\
        --outFilterMultimapNmax 20 
    """
} 
