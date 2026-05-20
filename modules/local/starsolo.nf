// modules/local/starsolo.nf
process STAR_SOLO {

    tag "$meta.id"
    errorStrategy 'ignore'
    publishDir "$params.outdir/starsolo/${meta.id}/logs", \
        mode: 'copy', \
        pattern: "${meta.id}.Log.final.out", \
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
    tuple val(meta), path(reads1), path(reads2), path(genomeDir)

    output:
    // tuple val(meta), path("${meta.id}.Solo.out")
    // path "${meta.id}.Log.final.out",                 emit: log
    // path "${meta.id}.Aligned.sortedByCoord.out.bam", emit: bam

    tuple val(meta), path("${meta.id}.Solo.out"), emit: solo_out_dir
    path "${meta.id}.Log.final.out", emit: log

    script:
    def solo_features = params.run_mode.toLowerCase() == 'velocity' ? 'Gene GeneFull Velocyto' : 'GeneFull'
    def active_whitelist = params.run_mode.toLowerCase() == 'multiome' ? params.multiome_rna_whitelist : params.soloCBwhitelist
    """
    STAR \\
        --runThreadN ${task.cpus} \\
        --genomeDir ${genomeDir} \\
        --readFilesIn ${reads2 instanceof List ? reads2.join(',') : reads2} ${reads1 instanceof List ? reads1.join(',') : reads1} \\
        --readFilesCommand zcat \\
        --outFileNamePrefix ${meta.id}. \\
        --soloType ${params.soloType} \\
        --soloCBwhitelist ${active_whitelist} \\
        --soloUMIfiltering ${params.soloUMIfiltering} \\
        --soloUMIdedup ${params.soloUMIdedup} \\
        --soloUMIlen ${params.soloUMIlen} \\
        --outSAMtype None \\
        --clipAdapterType CellRanger4 \\
        --outFilterScoreMin 30 \\
        --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts \\
        --soloCellFilter EmptyDrops_CR 10000 0.99 10 45000 90000 500 0.01 20000 0.01 10000 \\
        --soloFeatures ${solo_features} \\
        --soloMultiMappers Uniform \\
        --outFilterMultimapNmax 20
    """
}
