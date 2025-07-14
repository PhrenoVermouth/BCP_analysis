// modules/local/starsolo.nf
process STAR_SOLO {
    tag "$meta.id"
    publishDir "$params.outdir/starsolo/${meta.id}", mode: 'copy'

    input:
    tuple val(meta), path(reads)
    path genomeDir

    output:
    tuple val(meta), path("${meta.id}.Solo.out")
    path "Log.final.out", emit: log

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
        --outSAMtype BAM SortedByCoordinate
    """
} 
