#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Previously compiled scripts
include { STAR_SOLO } from './modules/local/starsolo'
include { GZIP_SOLO_OUTPUT } from './modules/local/gzip_solout' 
include { SCANPY_QC } from './modules/local/scanpy_qc'

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
    // gyang0716:gzip STARsolo output 
    GZIP_SOLO_OUTPUT(STAR_SOLO.out[0])
    // 3. Run Scanpy QC
    SCANPY_QC(GZIP_SOLO_OUTPUT.out.gzipped_dir)

    // 4. Run MultiQC for visualization
    Channel
        .from(STAR_SOLO.out.log)
        .mix(SCANPY_QC.out.qc_metrics)
        .collect()
        .set{ ch_multiqc_files }
    
    // MultiQC commands
    process MULTIQC {
        publishDir "$params.outdir/multiqc", mode: 'copy'

        input:
        path '*'

        output:
        path "multiqc_report.html"

        script:
        """
        multiqc .
        """
    }

    MULTIQC(ch_multiqc_files)
}

