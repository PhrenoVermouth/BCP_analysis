#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Previously compiled scripts
include { STAR_SOLO } from './modules/local/starsolo'
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

    // 3. Run Scanpy QC
    SCANPY_QC(STAR_SOLO.out[0])

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

