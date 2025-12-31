#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Previously compiled scripts
include { STAR_SOLO } from './modules/local/starsolo'
include { GZIP_SOLO_OUTPUT } from './modules/local/gzip_soloout'
include { SCRUBLET } from './modules/local/scrublet'
include { SOUPX } from './modules/local/soupx'
include { SAM_QC } from './modules/local/sam_qc'
include { ADD_VELOCITY_LAYERS } from './modules/local/add_velocity_layers'
include { MULTIQC } from './modules/local/multiqc'
include { METAQC_MERGE } from './modules/local/metaqc_merge'

def loadMitoMaxOverrides(String mitoMapPath) {
    def overrides = [:]
    if (!mitoMapPath) {
        return overrides
    }

    def mitoFile = new File(mitoMapPath)
    if (!mitoFile.exists()) {
        log.error "mito_max_map file not found: ${mitoFile}"
        System.exit(1)
    }

    mitoFile.text.readLines().eachWithIndex { line, idx ->
        def trimmed = line.trim()
        if (!trimmed || trimmed.startsWith('#')) {
            return
        }

        def parts = trimmed.split('=', 2)
        if (parts.size() != 2) {
            log.warn "Skipping malformed mito_max_map line ${idx + 1}: ${line}"
            return
        }

        def valueStr = parts[1].trim()
        Double mitoValue
        try {
            mitoValue = valueStr.toDouble()
        } catch (Exception e) {
            log.warn "Skipping mito_max_map line ${idx + 1}: invalid mito max '${valueStr}'"
            return
        }

        parts[0].split(',').each { token ->
            def sampleName = token.trim()
            if (sampleName) {
                overrides[sampleName] = mitoValue
            }
        }
    }

    log.info "Loaded mito_max overrides for ${overrides.size()} sample IDs from ${mitoFile}."
    return overrides
}

workflow {
    def runMode = params.run_mode.toLowerCase()
    if (!['genefull', 'velocity'].contains(runMode)) {
        log.error "Invalid params.run_mode: ${runMode}. Choose either 'genefull' or 'velocity'."
        System.exit(1)
    }

    def mitoMaxOverrides = loadMitoMaxOverrides(params.mito_max_map as String)
    // 1. Sample input
    Channel
        .fromPath(params.input)
        .splitCsv(header:true)
        .map { row ->
            def meta = [:]
            meta.id = row.sample
            def sampleMaxMito = mitoMaxOverrides.containsKey(meta.id) ? mitoMaxOverrides[meta.id] : params.max_mito
            meta.max_mito = sampleMaxMito
            if (mitoMaxOverrides.containsKey(meta.id)) {
                log.info "Using mito_max ${sampleMaxMito} for sample ${meta.id} from mito_max_map."
            }

            def reads = [ file(row.fastq_1), file(row.fastq_2) ]
            def defaultCounts = file("${params.outdir}/soupx/${row.sample}/${row.sample}_rm_ambient.h5ad")
            def countsAdata = row.counts_h5ad ? file(row.counts_h5ad) : defaultCounts
            return [ meta, reads, file(row.genomeDir), countsAdata ]
        }
        .set { ch_samples }

    ch_samples
        .map { meta, reads, genomeDir, countsAdata -> [ meta, reads, genomeDir ] }
        .set { ch_input_reads }

    // 2. Run STARsolo
    STAR_SOLO(ch_input_reads)

    // 2.1 gyang0721:gzip STARsolo output
    GZIP_SOLO_OUTPUT(STAR_SOLO.out.solo_out_dir)
        def ch_for_multiqc = STAR_SOLO.out.log

    if (runMode == 'genefull') {
        // 2.2 Run Scrublet before SoupX
        SCRUBLET(GZIP_SOLO_OUTPUT.out.gzipped_dir)
        // 2.3 SoupX ambient RNA removal
        SOUPX(GZIP_SOLO_OUTPUT.out.gzipped_dir.join(SCRUBLET.out.whitelist))

        METAQC_MERGE(SCRUBLET.out.metaqc_partial.join(SOUPX.out.rho))

        // 3. Run SAM QC using whitelist from Scrublet
        def soupx_h5ad_output = params.bypass_soupX ? SOUPX.out.pre_h5ad : SOUPX.out.ambient_h5ad
        SAM_QC(soupx_h5ad_output.join(SCRUBLET.out.whitelist))

        ch_for_multiqc = ch_for_multiqc
            .mix(METAQC_MERGE.out.metaqc_table.map { it[1] })
            .mix(SCRUBLET.out.qc_plots.map { it[1] })
            .mix(SAM_QC.out.qc_plots.map { it[1] })
            .mix(SOUPX.out.ambient_plot)
            .mix(SOUPX.out.contamination_plot)

    } else {
        ch_samples
            .map { meta, reads, genomeDir, countsAdata ->
                if (!countsAdata.exists()) {
                    throw new IllegalArgumentException(
                        "counts_h5ad is required for velocity mode for sample ${meta.id}. " +
                        "Provide counts_h5ad in samples.csv or ensure GeneFull output exists at " +
                        "${params.outdir}/soupx/${meta.id}/${meta.id}_rm_ambient.h5ad"
                    )
                }
                [ meta, countsAdata ]
            }
            .set { ch_counts_h5ad }

        ADD_VELOCITY_LAYERS(GZIP_SOLO_OUTPUT.out.gzipped_dir.join(ch_counts_h5ad))
    }
    // 4.1 Collect all the report files that need to be summarized
    ch_for_multiqc = ch_for_multiqc.collect()

    // 4.2 Create a channel pointing to the configuration file
    ch_multiqc_config = Channel.fromPath("${baseDir}/multiqc_config.yaml")

    // 4.3 Call MULTIQC
    MULTIQC(
        ch_for_multiqc,
        ch_multiqc_config
    )

}
