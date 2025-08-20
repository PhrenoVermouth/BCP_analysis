#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Previously compiled scripts
include { STAR_SOLO } from './modules/local/starsolo'
include { GZIP_SOLO_OUTPUT } from './modules/local/gzip_soloout'
include { SCRUBLET } from './modules/local/scrublet'
include { SOUPX } from './modules/local/soupx'
include { SAM_QC } from './modules/local/sam_qc'
include { MULTIQC } from './modules/local/multiqc'

workflow {
    // 1. Sample input
    Channel
        .fromPath(params.input)
        .splitCsv(header:true)
        .map { row ->
            def meta = [:]
            meta.id = row.sample

            def reads = [ file(row.fastq_1), file(row.fastq_2) ]
            return [ meta, reads, file(row.genomeDir) ]
        }
        .set { ch_input_reads }

    // 2. Run STARsolo
    STAR_SOLO(ch_input_reads)

    // 2.1 gyang0721:gzip STARsolo output
    GZIP_SOLO_OUTPUT(STAR_SOLO.out.solo_out_dir)
    // 2.2 Run Scrublet before SoupX
    SCRUBLET(GZIP_SOLO_OUTPUT.out.gzipped_dir)
    // 2.3 SoupX ambient RNA removal
    SOUPX(GZIP_SOLO_OUTPUT.out.gzipped_dir)

    // 3. Run SAM QC using whitelist from Scrublet
    SAM_QC(SOUPX.out.corrected_h5ad.join(SCRUBLET.out.whitelist))

    // 4.1 Collect all the report files that need to be summarized
    ch_for_multiqc = Channel.empty()
        .mix(STAR_SOLO.out.log)
        .mix(SCRUBLET.out.qc_cells_metrics)
        .mix(SCRUBLET.out.qc_counts_metrics)
        .mix(SCRUBLET.out.qc_genes_metrics)
        .mix(SCRUBLET.out.qc_plots.map { it[1] })
        .mix(SAM_QC.out.qc_plots.map { it[1] })
        .mix(SOUPX.out.ambient_plot)
        .mix(SOUPX.out.contamination_plot)
        .collect()

    // 4.2 Create a channel pointing to the configuration file
    ch_multiqc_config = Channel.fromPath("${baseDir}/multiqc_config.yaml")

    // 4.3 Call MULTIQC
    MULTIQC(
        ch_for_multiqc,
        ch_multiqc_config
    )

}
