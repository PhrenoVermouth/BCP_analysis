#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Previously compiled scripts
include { STAR_SOLO } from './modules/local/starsolo'
include { GZIP_SOLO_OUTPUT } from './modules/local/gzip_soloout' 
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
    // 2.2 gyang0721: SoupX
    SOUPX(GZIP_SOLO_OUTPUT.out.gzipped_dir) 

    // 3. Run SAM QC
    SAM_QC(SOUPX.out.corrected_h5ad)

    // 4.1 Collect all the report files that need to be summarized
    ch_for_multiqc = Channel.empty()
        .mix(STAR_SOLO.out.log)
        .mix(SAM_QC.out.qc_cells_metrics)
        .mix(SAM_QC.out.qc_counts_metrics)
        .mix(SAM_QC.out.qc_genes_metrics)
        .mix(SAM_QC.out.qc_plots)
        .collect()

    // 4.2 Create a channel pointing to the configuration file
    ch_multiqc_config = Channel.fromPath("${baseDir}/multiqc_config.yaml")

    // 4.3 Call MULTIQC
    MULTIQC(
        ch_for_multiqc,
        ch_multiqc_config
    )
    

}

