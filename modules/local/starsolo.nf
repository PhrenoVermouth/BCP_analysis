// modules/local/starsolo.nf
process STAR_SOLO {

    tag "$meta.id"

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

    // — send the STAR log to …/logs
   // publishDir "$params.outdir/starsolo/${meta.id}/logs", \
    //    mode: 'copy', \
     //   pattern: "${meta.id}.Log.final.out"

    // — send the sorted BAM to …/bam
    // publishDir "$params.outdir/starsolo/${meta.id}/bam", \
      //  mode: 'copy', \
       // pattern: "${meta.id}.Aligned.sortedByCoord.out.bam"

    input:
    tuple val(meta), path(reads), path(genomeDir)

    output:
    tuple val(meta), path("${meta.id}.Solo.out")
    path "${meta.id}.Log.final.out",                 emit: log
    path "${meta.id}.Aligned.sortedByCoord.out.bam", emit: bam

    script:
    """
    STAR \
        --runThreadN ${task.cpus} \
        --genomeDir ${genomeDir} \
        --readFilesIn ${reads[1]} ${reads[0]} \
        --readFilesCommand zcat \
        --outFileNamePrefix ${meta.id}. \
        --soloType ${params.soloType} \
        --soloCBwhitelist ${params.soloCBwhitelist} \
        --soloUMIfiltering ${params.soloUMIfiltering} \
        --soloUMIdedup ${params.soloUMIdedup} \
        --soloUMIlen ${params.soloUMIlen} \
        --outSAMtype BAM SortedByCoordinate \
        --clipAdapterType CellRanger4 \
        --outFilterScoreMin 30 \
        --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts \
        --soloCellFilter EmptyDrops_CR 10000 0.99 10 45000 90000 500 0.01 20000 0.01 10000 \
        --soloFeatures GeneFull \
        --soloMultiMappers EM \
        --outSAMattributes NH HI nM AS CR UR CB UB GX GN sS sQ sM \
        --outFilterMultimapNmax 20
    """
}
