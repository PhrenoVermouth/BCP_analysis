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
include { QC_PANEL } from './modules/local/qc_panel'
include { DEBUG_ANALYSIS } from './modules/local/debug_analysis'
import java.io.FilenameFilter

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

def validateRunMode() {
    def runMode = params.run_mode.toLowerCase()
    if (!['genefull', 'velocity'].contains(runMode)) {
        log.error "Invalid params.run_mode: ${runMode}. Choose either 'genefull' or 'velocity'."
        System.exit(1)
    }
    return runMode
}

def prepareSamples() {
    def runMode = validateRunMode()
    def mitoMaxOverrides = loadMitoMaxOverrides(params.mito_max_map as String)

    def ch_samples = Channel
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
            [ meta, reads, file(row.genomeDir), countsAdata ]
        }

    [ runMode: runMode, ch_samples: ch_samples ]
}

def createInputReads(ch_samples) {
    ch_samples.map { meta, reads, genomeDir, countsAdata -> [ meta, reads, genomeDir ] }
}

def createCountsH5ad(ch_samples) {
    ch_samples.map { meta, reads, genomeDir, countsAdata ->
        if (!countsAdata.exists()) {
            throw new IllegalArgumentException(
                "counts_h5ad is required for velocity mode for sample ${meta.id}. " +
                "Provide counts_h5ad in samples.csv or ensure GeneFull output exists at " +
                "${params.outdir}/soupx/${meta.id}/${meta.id}_rm_ambient.h5ad"
            )
        }
        [ meta, countsAdata ]
    }
}

def starLogsFromPublish(ch_samples) {
    ch_samples.map { meta, reads, genomeDir, countsAdata ->
        def logPath = file("${params.outdir}/starsolo/${meta.id}/logs/${meta.id}.Log.final.out")
        if (!logPath.exists()) {
            throw new IllegalArgumentException("STARsolo log not found for sample ${meta.id} at ${logPath}")
        }
        logPath
    }
}

def gzippedDirsFromPublish(ch_samples) {
    ch_samples.map { meta, reads, genomeDir, countsAdata ->
        def gzDir = file("${params.outdir}/gzipped_matrix/${meta.id}/${meta.id}.Solo.out")
        def fallback = file("${params.outdir}/gzipped_matrix/${meta.id}")
        def chosen = gzDir.exists() ? gzDir : fallback
        if (!chosen.exists()) {
            throw new IllegalArgumentException("Gzipped STARsolo output not found for sample ${meta.id} under ${params.outdir}/gzipped_matrix")
        }
        [ meta, chosen ]
    }
}

def scrubletOutputsFromPublish(ch_samples) {
    def ch_whitelist = ch_samples.map { meta, reads, genomeDir, countsAdata ->
        def whitelist = file("${params.outdir}/scrublet/${meta.id}/${meta.id}_whitelist.txt")
        if (!whitelist.exists()) {
            throw new IllegalArgumentException("Scrublet whitelist not found for sample ${meta.id} at ${whitelist}")
        }
        [ meta, whitelist ]
    }

    def ch_metaqc_partial = ch_samples.map { meta, reads, genomeDir, countsAdata ->
        def pattern = "${params.outdir}/scrublet/${meta.id}/${meta.id}_total_metaqc_partial.tsv"
        def metaqcPartial = file(pattern)
        if (!metaqcPartial.exists()) {
            throw new IllegalArgumentException("Scrublet metaqc_partial not found for sample ${meta.id} at ${pattern}")
        }
        [ meta, metaqcPartial ]
    }

    def ch_qc_plots = ch_samples
        .map { meta, reads, genomeDir, countsAdata ->
            def plotDir = new File("${params.outdir}/scrublet/${meta.id}")
            def qcPlots = plotDir.listFiles({ d, name -> name.toLowerCase().endsWith('.png') } as FilenameFilter) ?: []
            if (!qcPlots) {
                throw new IllegalArgumentException("Scrublet QC plots not found for sample ${meta.id} in ${plotDir}")
            }
            qcPlots.collect { plot -> [ meta, plot ] }
        }
        .flatten()

    [ whitelist: ch_whitelist, metaqc_partial: ch_metaqc_partial, qc_plots: ch_qc_plots ]
}

def soupxOutputsFromPublish(ch_samples) {
    def buildPath = { meta, suffix ->
        def path = file("${params.outdir}/soupx/${meta.id}/${suffix}")
        if (!path.exists()) {
            throw new IllegalArgumentException("SoupX output ${suffix} not found for sample ${meta.id} at ${path}")
        }
        path
    }

    def ch_pre_h5ad = ch_samples.map { meta, reads, genomeDir, countsAdata -> [ meta, buildPath(meta, "${meta.id}_pre_soupx.h5ad") ] }
    def ch_ambient_h5ad = ch_samples.map { meta, reads, genomeDir, countsAdata -> [ meta, buildPath(meta, "${meta.id}_rm_ambient.h5ad") ] }
    def ch_rho = ch_samples.map { meta, reads, genomeDir, countsAdata -> [ meta, buildPath(meta, "${meta.id}_soupx_rho.tsv") ] }
    def ch_ambient_plot = ch_samples.map { meta, reads, genomeDir, countsAdata -> buildPath(meta, "0.${meta.id}_ambient_RNA_removed_mqc.png") }
    def ch_contamination_plot = ch_samples.map { meta, reads, genomeDir, countsAdata -> buildPath(meta, "0.${meta.id}_soupx_contamination_estimation_mqc.png") }
    def ch_combined_plot = ch_samples.map { meta, reads, genomeDir, countsAdata -> buildPath(meta, "0.${meta.id}_soupx_combined_mqc.png") }

    [
        pre_h5ad: ch_pre_h5ad,
        ambient_h5ad: ch_ambient_h5ad,
        rho: ch_rho,
        ambient_plot: ch_ambient_plot,
        contamination_plot: ch_contamination_plot,
        combined_plot: ch_combined_plot
    ]
}

def runGenefullFromGzip(gzippedDir, ch_for_multiqc) {
    def scrubletResults = SCRUBLET(gzippedDir)
    def soupxResults = SOUPX(gzippedDir.join(scrubletResults.out.whitelist))

    def metaqcMergeResults = METAQC_MERGE(scrubletResults.out.metaqc_partial.join(soupxResults.out.rho))

    def soupx_h5ad_output = params.bypass_soupX ? soupxResults.out.pre_h5ad : soupxResults.out.ambient_h5ad
    def samQcResults = SAM_QC(soupx_h5ad_output.join(scrubletResults.out.whitelist))

    def ch_scrublet_plots = scrubletResults.out.qc_plots.groupTuple().map { meta, plots ->
        def finder = { List plotList, String token ->
            def match = plotList.find { it.getName().toString().contains(token) }
            if (!match) {
                throw new IllegalStateException("Missing ${token} plot for sample ${meta.id}")
            }
            match
        }
        def knee = finder(plots, '_knee_plot_mqc')
        def histogram = finder(plots, '_doublet_score_histogram_QC1_mqc')
        def violin = finder(plots, '_violin_comparison_QC1_mqc')
        def mito = finder(plots, '_violin_mito_filtering_QC2_mqc')
        [ meta, knee, histogram, violin, mito ]
    }

    def ch_samqc_plots = samQcResults.out.qc_plots.groupTuple().map { meta, plots ->
        def finder = { List plotList, String token ->
            def match = plotList.find { it.getName().toString().contains(token) }
            if (!match) {
                throw new IllegalStateException("Missing ${token} plot for sample ${meta.id}")
            }
            match
        }
        def umap = finder(plots, '_umap_leiden_QC2_mqc')
        def dotplot = finder(plots, '_marker_genes_dotplot_mqc')
        [ meta, umap, dotplot ]
    }

    QC_PANEL(
        ch_scrublet_plots
            .join(soupxResults.out.combined_plot)
            .join(ch_samqc_plots)
            .map { meta, knee, histogram, violin, mito, soupx_combined, umap, dotplot ->
                [ meta, knee, histogram, violin, mito, soupx_combined, umap, dotplot ]
            }
    )

    ch_for_multiqc
        .mix(metaqcMergeResults.out.metaqc_table.map { it[1] })
        .mix(QC_PANEL.out.panel.map { it[1] })
}

def finalizeMultiqc(ch_for_multiqc) {
    def collected = ch_for_multiqc.collect()
    def ch_multiqc_config = Channel.fromPath("${baseDir}/multiqc_config.yaml")

    MULTIQC(
        collected,
        ch_multiqc_config
    )

    [
        report: MULTIQC.out.report,
        data: MULTIQC.out.data
    ]
}

def runDebugAnalysis(ch_multiqc_report) {
    if (params.debug_off) {
        log.info "Skipping debug and AI analysis because --debug_off is set."
        return
    }

    def ch_trace = Channel.value(file(params.trace_file))
    def ch_workdir = Channel.value(workflow.workDir.toString())
    def ch_base_url = Channel.value(params.debug_base_url)
    def ch_api_key = Channel.value(params.debug_api_key ?: '')

    DEBUG_ANALYSIS(
        ch_multiqc_report,
        ch_trace,
        ch_workdir,
        ch_base_url,
        ch_api_key
    )
}

workflow {
    run_from_start()
}

workflow run_from_start {
    def context = prepareSamples()
    def ch_samples = context.ch_samples
    def runMode = context.runMode

    def ch_input_reads = createInputReads(ch_samples)

    STAR_SOLO(ch_input_reads)
    GZIP_SOLO_OUTPUT(STAR_SOLO.out.solo_out_dir)
    def ch_for_multiqc = STAR_SOLO.out.log

    if (runMode == 'genefull') {
        ch_for_multiqc = runGenefullFromGzip(GZIP_SOLO_OUTPUT.out.gzipped_dir, ch_for_multiqc)
    } else {
        def ch_counts_h5ad = createCountsH5ad(ch_samples)
        ADD_VELOCITY_LAYERS(GZIP_SOLO_OUTPUT.out.gzipped_dir.join(ch_counts_h5ad))
    }

    def multiqcOutputs = finalizeMultiqc(ch_for_multiqc)
    runDebugAnalysis(multiqcOutputs.report)
}

workflow post_gzip_entry {
    def context = prepareSamples()
    def ch_samples = context.ch_samples
    def runMode = context.runMode

    def ch_for_multiqc = starLogsFromPublish(ch_samples)
    def ch_gzipped = gzippedDirsFromPublish(ch_samples)

    if (runMode == 'genefull') {
        ch_for_multiqc = runGenefullFromGzip(ch_gzipped, ch_for_multiqc)
    } else {
        def ch_counts_h5ad = createCountsH5ad(ch_samples)
        ADD_VELOCITY_LAYERS(ch_gzipped.join(ch_counts_h5ad))
    }

    def multiqcOutputs = finalizeMultiqc(ch_for_multiqc)
    runDebugAnalysis(multiqcOutputs.report)
}

workflow post_soupx_entry {
    def context = prepareSamples()
    def ch_samples = context.ch_samples
    def runMode = context.runMode

    if (runMode != 'genefull') {
        log.error "post_soupx_entry is only available in genefull run_mode."
        System.exit(1)
    }

    def ch_for_multiqc = starLogsFromPublish(ch_samples)
    def scrubletChannels = scrubletOutputsFromPublish(ch_samples)
    def soupxChannels = soupxOutputsFromPublish(ch_samples)

    METAQC_MERGE(scrubletChannels.metaqc_partial.join(soupxChannels.rho))

    def soupx_h5ad_output = params.bypass_soupX ? soupxChannels.pre_h5ad : soupxChannels.ambient_h5ad
    SAM_QC(soupx_h5ad_output.join(scrubletChannels.whitelist))

    def ch_scrublet_plots = scrubletChannels.qc_plots.groupTuple().map { meta, plots ->
        def finder = { List plotList, String token ->
            def match = plotList.find { it.getName().toString().contains(token) }
            if (!match) {
                throw new IllegalStateException("Missing ${token} plot for sample ${meta.id}")
            }
            match
        }
        def knee = finder(plots, '_knee_plot_mqc')
        def histogram = finder(plots, '_doublet_score_histogram_QC1_mqc')
        def violin = finder(plots, '_violin_comparison_QC1_mqc')
        def mito = finder(plots, '_violin_mito_filtering_QC2_mqc')
        [ meta, knee, histogram, violin, mito ]
    }

    def ch_samqc_plots = SAM_QC.out.qc_plots.groupTuple().map { meta, plots ->
        def finder = { List plotList, String token ->
            def match = plotList.find { it.getName().toString().contains(token) }
            if (!match) {
                throw new IllegalStateException("Missing ${token} plot for sample ${meta.id}")
            }
            match
        }
        def umap = finder(plots, '_umap_leiden_QC2_mqc')
        def dotplot = finder(plots, '_marker_genes_dotplot_mqc')
        [ meta, umap, dotplot ]
    }

    QC_PANEL(
        ch_scrublet_plots
            .join(soupxChannels.combined_plot)
            .join(ch_samqc_plots)
            .map { meta, knee, histogram, violin, mito, soupx_combined, umap, dotplot ->
                [ meta, knee, histogram, violin, mito, soupx_combined, umap, dotplot ]
            }
    )

    ch_for_multiqc
        .mix(METAQC_MERGE.out.metaqc_table.map { it[1] })
        .mix(QC_PANEL.out.panel.map { it[1] })

    def multiqcOutputs = finalizeMultiqc(ch_for_multiqc)
    runDebugAnalysis(multiqcOutputs.report)
}
