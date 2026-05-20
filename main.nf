#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// ── RNA / shared modules ───────────────────────────────────────────────────
include { STAR_SOLO           } from './modules/local/starsolo'
include { GZIP_SOLO_OUTPUT    } from './modules/local/gzip_soloout'
include { SCRUBLET            } from './modules/local/scrublet'
include { SOUPX               } from './modules/local/soupx'
include { SAM_QC              } from './modules/local/sam_qc'
include { ADD_VELOCITY_LAYERS } from './modules/local/add_velocity_layers'
include { MULTIQC             } from './modules/local/multiqc'
include { METAQC_MERGE        } from './modules/local/metaqc_merge'
include { QC_PANEL            } from './modules/local/qc_panel'
include { DEBUG_ANALYSIS      } from './modules/local/debug_analysis'
// ── Multiome-specific modules ─────────────────────────────────────────────
include { CHROMAP             } from './modules/local/chromap'
include { ARCHR_INTEGRATION   } from './modules/local/archr_integration'

import java.io.FilenameFilter

// ──────────────────────────────────────────────────────────────────────────
// Helper: load per-sample numeric overrides from a "key=value" file.
// Supports comma-separated keys (e.g.  pcf_1,pcf_2=0.15) and
// substring-match fallback for convenience.
// Shared by --mito_max_map and --scrublet_threshold_map.
// ──────────────────────────────────────────────────────────────────────────
def loadKeyValueOverrides(String mapPath) {
    def overrides = [:]
    if (!mapPath) return overrides

    def mapFile = new File(mapPath)
    if (!mapFile.exists()) {
        log.error "Override file not found: ${mapFile}"
        System.exit(1)
    }

    mapFile.text.readLines().eachWithIndex { line, idx ->
        def trimmed = line.trim()
        if (!trimmed || trimmed.startsWith('#')) return

        def parts = trimmed.split('=', 2)
        if (parts.size() != 2) {
            log.warn "Skipping malformed line ${idx + 1}: ${line}"
            return
        }
        def valueStr = parts[1].trim()
        Double val
        try { val = valueStr.toDouble() }
        catch (Exception e) {
            log.warn "Skipping line ${idx + 1}: invalid numeric value '${valueStr}'"
            return
        }
        parts[0].split(',').each { token ->
            def k = token.trim()
            if (k) overrides[k] = val
        }
    }
    log.info "Loaded ${overrides.size()} overrides from ${mapFile}."
    return overrides
}

// Resolve a per-sample value with exact → substring-longest → null fallback.
def resolveOverride(Map overrides, String sampleId) {
    if (overrides.containsKey(sampleId)) return [ sampleId, overrides[sampleId] ]
    def hits = overrides.findAll { k, v -> sampleId.contains(k) }
    if (!hits) return [ null, null ]
    def bestKey = hits.keySet().toList().sort { -it.size() }[0]
    return [ bestKey, hits[bestKey] ]
}

// ──────────────────────────────────────────────────────────────────────────
// Validate run_mode at startup
// ──────────────────────────────────────────────────────────────────────────
def validateRunMode() {
    def runMode = params.run_mode.toLowerCase()
    if (!['genefull', 'velocity', 'multiome'].contains(runMode)) {
        log.error "Invalid params.run_mode: '${runMode}'. Choose 'genefull', 'velocity', or 'multiome'."
        System.exit(1)
    }
    return runMode
}

// ──────────────────────────────────────────────────────────────────────────
// Build the RNA sample channel from samples.csv.
// Tuple emitted: [meta, fq1s, fq2s, genomeDir, countsAdata]
//   fq1s / fq2s are lists (supports per-lane multi-FASTQ).
// ──────────────────────────────────────────────────────────────────────────
def prepareSamples() {
    def runMode                    = validateRunMode()
    def mitoMaxOverrides           = loadKeyValueOverrides(params.mito_max_map as String)
    def scrubletThresholdOverrides = loadKeyValueOverrides(params.scrublet_threshold_map as String)

    def ch_samples = Channel
        .fromPath(params.input)
        .splitCsv(header: true)
        .map { row ->
            // Required column validation
            if (!row.sample?.toString()?.trim() ||
                !row.fastq_1?.toString()?.trim() ||
                !row.fastq_2?.toString()?.trim() ||
                !row.genomeDir?.toString()?.trim()) {
                throw new IllegalArgumentException(
                    "Invalid samples.csv row (required columns: sample, fastq_1, fastq_2, genomeDir): ${row}")
            }

            def meta   = [:]
            meta.id    = row.sample.toString().trim()

            // Per-sample mito_max (exact → substring → global default)
            def (mitoKey, mitoVal) = resolveOverride(mitoMaxOverrides, meta.id)
            meta.max_mito = (mitoVal != null) ? mitoVal : params.max_mito
            if (mitoKey != null)
                log.info "mito_max ${meta.max_mito} for ${meta.id} (key '${mitoKey}')."

            // Per-sample scrublet threshold (map → global param → null = auto)
            def (scrubKey, scrubVal) = resolveOverride(scrubletThresholdOverrides, meta.id)
            Double scrubletThreshold = scrubVal
            if (scrubletThreshold == null && params.scrublet_manual_threshold != null)
                scrubletThreshold = params.scrublet_manual_threshold as Double
            if (scrubletThreshold != null) {
                meta.scrublet_threshold = scrubletThreshold
                log.info "scrublet_threshold ${scrubletThreshold} for ${meta.id} (key '${scrubKey}')."
            }

            def defaultCounts = file("${params.outdir}/sam_qc/${meta.id}/${meta.id}_filtered_QC2.h5ad")
            def countsAdata   = row.counts_h5ad ? file(row.counts_h5ad) : defaultCounts

            [ meta,
              file(row.fastq_1.toString().trim()),
              file(row.fastq_2.toString().trim()),
              file(row.genomeDir.toString().trim()),
              countsAdata ]
        }
        .groupTuple(by: 0)
        // After groupTuple: [meta, [fq1…], [fq2…], [genomeDir…], [countsAdata…]]
        .map { meta, fq1s, fq2s, genomeDirs, countsAdatas ->
            [ meta, fq1s, fq2s, genomeDirs[0], countsAdatas[0] ]
        }

    [ runMode: runMode, ch_samples: ch_samples ]
}

// ──────────────────────────────────────────────────────────────────────────
// Build ATAC sample channel from samples.csv (multiome only).
// Requires columns: atac_R1, atac_R2, atac_R3
//   R1 = genomic read, R2 = 16-bp barcode, R3 = genomic paired-end
// Genome-level resources (ref_fasta, chromap_index) come from global params.
// ──────────────────────────────────────────────────────────────────────────
def prepareATACSamples() {
    Channel
        .fromPath(params.input)
        .splitCsv(header: true)
        .map { row ->
            if (!row.atac_fastq_1?.toString()?.trim() ||
                !row.atac_fastq_2?.toString()?.trim() ||
                !row.atac_fastq_3?.toString()?.trim()) {
                throw new IllegalArgumentException(
                    "Multiome mode requires atac_fastq_1, atac_fastq_2, atac_fastq_3 columns in samples.csv. Missing for: ${row.sample}")
            }
            if (!row.atac_ref?.toString()?.trim() ||
                !row.atac_index?.toString()?.trim()) {
                throw new IllegalArgumentException(
                    "Multiome mode requires atac_ref and atac_index columns in samples.csv. Missing for: ${row.sample}")
            }
            def meta = [id: row.sample.toString().trim()]
            [ meta,
              file(row.atac_fastq_1.toString().trim()),
              file(row.atac_fastq_2.toString().trim()),
              file(row.atac_fastq_3.toString().trim()),
              file(row.atac_ref.toString().trim()),
              file(row.atac_index.toString().trim()) ]
        }
        .groupTuple(by: 0)
        // After groupTuple: [meta, [R1s...], [R2s...], [R3s...], [refs...], [idxs...]]
        // Collect per-read lists for multi-lane support (chromap accepts comma-separated inputs)
        .map { meta, r1s, r2s, r3s, refs, idxs -> [ meta, r1s, r2s, r3s, refs[0], idxs[0] ] }
}

// ──────────────────────────────────────────────────────────────────────────
// Channel extractors: convert ch_samples 5-tuple to task-specific shapes
// ──────────────────────────────────────────────────────────────────────────
def createInputReads(ch_samples) {
    // [meta, fq1s, fq2s, genomeDir]
    ch_samples.map { meta, fq1s, fq2s, genomeDir, countsAdata ->
        [ meta, fq1s, fq2s, genomeDir ]
    }
}

def createCountsH5ad(ch_samples) {
    ch_samples.map { meta, fq1s, fq2s, genomeDir, countsAdata ->
        if (!countsAdata.exists()) {
            throw new IllegalArgumentException(
                "counts_h5ad required for velocity mode for sample ${meta.id}. " +
                "Provide it in samples.csv or ensure SAM_QC output exists at " +
                "${params.outdir}/sam_qc/${meta.id}/${meta.id}_filtered_QC2.h5ad")
        }
        [ meta, countsAdata ]
    }
}

// ──────────────────────────────────────────────────────────────────────────
// Re-entry helpers: rebuild channels from already-published outputs
// ──────────────────────────────────────────────────────────────────────────
def starLogsFromPublish(ch_samples) {
    ch_samples.flatMap { meta, fq1s, fq2s, genomeDir, countsAdata ->
        def p = file("${params.outdir}/starsolo/${meta.id}/logs/${meta.id}.Log.final.out")
        if (!p.exists()) {
            log.warn "STARsolo log not found for ${meta.id} at ${p}. Skipping..."
            return []
        }
        [ p ]
    }
}

def gzippedDirsFromPublish(ch_samples) {
    ch_samples.flatMap { meta, fq1s, fq2s, genomeDir, countsAdata ->
        def gzDir   = file("${params.outdir}/gzipped_matrix/${meta.id}/${meta.id}.Solo.out")
        def fallback = file("${params.outdir}/gzipped_matrix/${meta.id}")
        def chosen  = gzDir.exists() ? gzDir : fallback
        if (!chosen.exists()) {
            log.warn "Gzipped STARsolo output not found for ${meta.id}. Skipping..."
            return []
        }
        [ [ meta, chosen ] ]
    }
}

def scrubletOutputsFromPublish(ch_samples) {
    def ch_whitelist = ch_samples.flatMap { meta, fq1s, fq2s, genomeDir, countsAdata ->
        def p = file("${params.outdir}/scrublet/${meta.id}/${meta.id}_whitelist.txt")
        if (!p.exists()) { log.warn "Scrublet whitelist not found for ${meta.id}. Skipping..."; return [] }
        [ [ meta, p ] ]
    }
    def ch_metaqc_partial = ch_samples.flatMap { meta, fq1s, fq2s, genomeDir, countsAdata ->
        def p = file("${params.outdir}/scrublet/${meta.id}/${meta.id}_total_metaqc_partial.tsv")
        if (!p.exists()) { log.warn "Scrublet metaqc_partial not found for ${meta.id}. Skipping..."; return [] }
        [ [ meta, p ] ]
    }
    def ch_qc_plots = ch_samples.flatMap { meta, fq1s, fq2s, genomeDir, countsAdata ->
        def plotDir = new File("${params.outdir}/scrublet/${meta.id}")
        def qcPlots = plotDir.listFiles({ d, name -> name.toLowerCase().endsWith('.png') } as FilenameFilter) ?: []
        if (!qcPlots) { log.warn "Scrublet QC plots not found for ${meta.id}. Skipping..."; return [] }
        qcPlots.collect { plot -> [ meta, file(plot.absolutePath) ] }
    }
    [ whitelist: ch_whitelist, metaqc_partial: ch_metaqc_partial, qc_plots: ch_qc_plots ]
}

def soupxOutputsFromPublish(ch_samples) {
    def tryPath = { meta, suffix ->
        def p = file("${params.outdir}/soupx/${meta.id}/${suffix}")
        if (!p.exists()) { log.warn "SoupX output ${suffix} not found for ${meta.id}. Skipping..."; return null }
        p
    }
    def ch_pre_h5ad         = ch_samples.flatMap { meta, fq1s, fq2s, gd, ca ->
        def p = tryPath(meta, "${meta.id}_pre_soupx.h5ad");           p ? [ [meta, p] ] : [] }
    def ch_ambient_h5ad     = ch_samples.flatMap { meta, fq1s, fq2s, gd, ca ->
        def p = tryPath(meta, "${meta.id}_rm_ambient.h5ad");          p ? [ [meta, p] ] : [] }
    def ch_rho              = ch_samples.flatMap { meta, fq1s, fq2s, gd, ca ->
        def p = tryPath(meta, "${meta.id}_soupx_rho.tsv");            p ? [ [meta, p] ] : [] }
    def ch_ambient_plot     = ch_samples.flatMap { meta, fq1s, fq2s, gd, ca ->
        def p = tryPath(meta, "0.${meta.id}_ambient_RNA_removed_mqc.png"); p ? [ p ] : [] }
    def ch_contamination_plot = ch_samples.flatMap { meta, fq1s, fq2s, gd, ca ->
        def p = tryPath(meta, "0.${meta.id}_soupx_contamination_estimation_mqc.png"); p ? [ p ] : [] }
    def ch_combined_plot    = ch_samples.flatMap { meta, fq1s, fq2s, gd, ca ->
        def p = tryPath(meta, "0.${meta.id}_soupx_combined_mqc.png"); p ? [ [meta, p] ] : [] }
    [
        pre_h5ad: ch_pre_h5ad, ambient_h5ad: ch_ambient_h5ad,
        rho: ch_rho, ambient_plot: ch_ambient_plot,
        contamination_plot: ch_contamination_plot, combined_plot: ch_combined_plot
    ]
}

// ──────────────────────────────────────────────────────────────────────────
// Shared utility: find a named plot in a (possibly nested) list
// ──────────────────────────────────────────────────────────────────────────
def makeFinder(meta) {
    return { List plotList, String token ->
        def match = plotList.flatten().find { plot ->
            def name = plot.respondsTo('getFileName') ? plot.getFileName().toString() : plot.toString()
            name.contains(token)
        }
        if (!match) throw new IllegalStateException("Missing plot token '${token}' for sample ${meta.id}")
        match
    }
}

// ──────────────────────────────────────────────────────────────────────────
// Sub-workflow: RNA QC from gzipped STARsolo output.
// Runs SCRUBLET → SOUPX → METAQC_MERGE → SAM_QC → QC_PANEL.
// Also used as the RNA-QC leg of the multiome workflow.
// ──────────────────────────────────────────────────────────────────────────
workflow runGenefullFromGzip {
    take:
        ch_gzipped      // [meta, solo_out_dir]
        ch_for_multiqc  // upstream files (e.g. STAR logs) to include in MultiQC

    main:
        SCRUBLET(ch_gzipped)
        SOUPX(ch_gzipped.join(SCRUBLET.out.whitelist))
        METAQC_MERGE(SCRUBLET.out.metaqc_partial.join(SOUPX.out.rho))

        def soupx_h5ad_output = params.bypass_soupX ? SOUPX.out.pre_h5ad : SOUPX.out.ambient_h5ad
        SAM_QC(soupx_h5ad_output.join(SCRUBLET.out.whitelist))

        ch_scrublet_plots = SCRUBLET.out.qc_plots.groupTuple().map { meta, plots ->
            def f = makeFinder(meta)
            [ meta,
              f(plots, '_knee_plot_mqc'),
              f(plots, '_doublet_score_histogram_QC1_mqc'),
              f(plots, '_violin_comparison_QC1_mqc'),
              f(plots, '_violin_mito_filtering_QC2_mqc') ]
        }

        ch_samqc_plots = SAM_QC.out.qc_plots.groupTuple().map { meta, plots ->
            def f = makeFinder(meta)
            [ meta,
              f(plots, '_umap_leiden_QC2_mqc'),
              f(plots, '_marker_genes_dotplot_mqc') ]
        }

        QC_PANEL(
            ch_scrublet_plots
                .join(SOUPX.out.combined_plot)
                .join(ch_samqc_plots)
                .map { meta, knee, histogram, violin, mito, soupx_combined, umap, dotplot ->
                    [ meta, knee, histogram, violin, mito, soupx_combined, umap, dotplot ]
                }
        )

        ch_multiqc_out = ch_for_multiqc
            .mix(METAQC_MERGE.out.metaqc_table.map { it[1] })
            .mix(QC_PANEL.out.panel.map { it[1] })

    emit:
        ch_multiqc = ch_multiqc_out
}

// ──────────────────────────────────────────────────────────────────────────
// Sub-workflow: ATAC alignment via Chromap.
// Encapsulates genome-resource channel construction so that this step can be
// launched at the top-level entry points, in parallel with STAR_SOLO.
// ──────────────────────────────────────────────────────────────────────────
workflow runChromapWorkflow {
    take:
        ch_atac   // [meta, [R1s], [R2s], [R3s], ref_fasta, chromap_index]

    main:
        ch_atac_whitelist = Channel.value(
            params.multiome_atac_whitelist
                ? file(params.multiome_atac_whitelist)
                : file('NO_FILE')
        )

        CHROMAP(
            ch_atac
                .combine(ch_atac_whitelist)
                .map { meta, r1s, r2s, r3s, ref, idx, wl -> [ meta, r1s, r2s, r3s, ref, idx, wl ] }
        )

    emit:
        fragments     = CHROMAP.out.fragments
        fragments_tbi = CHROMAP.out.fragments_tbi
        mqc_json      = CHROMAP.out.mqc_json
        chromap_log   = CHROMAP.out.log
}

// ──────────────────────────────────────────────────────────────────────────
// Sub-workflow: ArchR integration (RNA GeneFull dir + pre-built fragments).
//
// Receives already-computed fragment channels so it can be called only after
// both GZIP_SOLO_OUTPUT and runChromapWorkflow finish.
// The R script handles SoupX (ATAC-cluster-guided) and multi-modal doublet
// removal internally — no separate Scrublet/SoupX/SAM_QC is needed.
// ──────────────────────────────────────────────────────────────────────────
workflow runMultiome {
    take:
        ch_gzipped        // [meta, solo_out_dir]  from GZIP_SOLO_OUTPUT
        ch_fragments      // [meta, fragments.tsv.gz]  from runChromapWorkflow
        ch_fragments_tbi  // [meta, fragments.tsv.gz.tbi]  from runChromapWorkflow
        ch_star_logs      // STAR log files for MultiQC

    main:
        ch_barcode_translation = Channel.value(
            params.multiome_barcode_translation
                ? file(params.multiome_barcode_translation)
                : file('NO_FILE')
        )

        // Use meta.id as the explicit join key so minor differences in meta
        // map content between channels don't cause silent join failures.
        ch_rna_keyed  = ch_gzipped.map { meta, solo_dir ->
            [ meta.id, meta, file("${solo_dir}/GeneFull") ]
        }
        ch_frag_keyed = ch_fragments.map     { meta, frag -> [ meta.id, frag ] }
        ch_tbi_keyed  = ch_fragments_tbi.map { meta, tbi  -> [ meta.id, tbi  ] }

        ARCHR_INTEGRATION(
            ch_rna_keyed
                .join(ch_frag_keyed)
                .join(ch_tbi_keyed)
                .combine(ch_barcode_translation)
                .map { id, meta, rna_dir, frag, frag_tbi, wl ->
                    [ meta, rna_dir, frag, frag_tbi, wl ]
                }
        )

        // Collect STAR logs + ARCHR MultiQC outputs (plots, QC table)
        // CHROMAP mqc_json is mixed in at the entry-point level
        ch_multiqc_out = ch_star_logs
            .mix(ARCHR_INTEGRATION.out.mqc_plots.flatten())
            .mix(ARCHR_INTEGRATION.out.mqc_table)

    emit:
        ch_multiqc = ch_multiqc_out
}

// ──────────────────────────────────────────────────────────────────────────
// Sub-workflow: run MultiQC and produce final HTML report
// ──────────────────────────────────────────────────────────────────────────
workflow finalizeMultiqc {
    take:
        ch_for_multiqc

    main:
        collected         = ch_for_multiqc.collect()
        ch_multiqc_config = Channel.fromPath("${baseDir}/multiqc_config.yaml")
        MULTIQC(collected, ch_multiqc_config)

    emit:
        report = MULTIQC.out.report
        data   = MULTIQC.out.data
}

// ──────────────────────────────────────────────────────────────────────────
// Sub-workflow: optional AI-assisted debug analysis after MultiQC
// ──────────────────────────────────────────────────────────────────────────
workflow runDebugAnalysis {
    take:
        ch_multiqc_report

    main:
        if (!params.debug_off) {
            ch_trace    = Channel.value(file(params.trace_file))
            ch_workdir  = Channel.value(workflow.workDir.toString())
            ch_base_url = Channel.value(params.debug_base_url)
            ch_api_key  = Channel.value(params.debug_api_key ?: '')

            DEBUG_ANALYSIS(
                ch_multiqc_report,
                ch_trace,
                ch_workdir,
                ch_base_url,
                ch_api_key
            )
        }
}

// ══════════════════════════════════════════════════════════════════════════
// Entry-point: full run from raw FASTQs
// ══════════════════════════════════════════════════════════════════════════
workflow {
    run_from_start()
}

workflow run_from_start {
    main:
        context    = prepareSamples()
        ch_samples = context.ch_samples
        runMode    = context.runMode

        ch_input_reads = createInputReads(ch_samples)

        // ── RNA alignment (always) ─────────────────────────────────────
        STAR_SOLO(ch_input_reads)
        GZIP_SOLO_OUTPUT(STAR_SOLO.out.solo_out_dir)

        if (runMode == 'genefull') {
            runGenefullFromGzip(GZIP_SOLO_OUTPUT.out.gzipped_dir, STAR_SOLO.out.log)
            finalizeMultiqc(runGenefullFromGzip.out.ch_multiqc)

        } else if (runMode == 'multiome') {
            // ── ATAC alignment: runs in parallel with STAR_SOLO ────────
            // prepareATACSamples() reads the CSV immediately at launch;
            // runChromapWorkflow has no dependency on STAR_SOLO outputs,
            // so CHROMAP and STAR_SOLO execute concurrently.
            ch_atac = prepareATACSamples()
            runChromapWorkflow(ch_atac)

            // ── ArchR integration: waits for both alignments ───────────
            runMultiome(
                GZIP_SOLO_OUTPUT.out.gzipped_dir,
                runChromapWorkflow.out.fragments,
                runChromapWorkflow.out.fragments_tbi,
                STAR_SOLO.out.log
            )
            finalizeMultiqc(runMultiome.out.ch_multiqc.mix(runChromapWorkflow.out.mqc_json))

        } else {
            // velocity mode
            ch_counts_h5ad = createCountsH5ad(ch_samples)
            ADD_VELOCITY_LAYERS(GZIP_SOLO_OUTPUT.out.gzipped_dir.join(ch_counts_h5ad))
            ch_velocity_signal = ADD_VELOCITY_LAYERS.out.velocity_h5ad.map { it[1] }
            finalizeMultiqc(STAR_SOLO.out.log.mix(ch_velocity_signal))
        }

        runDebugAnalysis(finalizeMultiqc.out.report)
}

// ══════════════════════════════════════════════════════════════════════════
// Re-entry: resume from after STARsolo + gzip (skip STAR re-run)
// ══════════════════════════════════════════════════════════════════════════
workflow post_gzip_entry {
    main:
        context    = prepareSamples()
        ch_samples = context.ch_samples
        runMode    = context.runMode

        ch_star_logs = starLogsFromPublish(ch_samples)
        ch_gzipped   = gzippedDirsFromPublish(ch_samples)

        if (runMode == 'genefull') {
            runGenefullFromGzip(ch_gzipped, ch_star_logs)
            finalizeMultiqc(runGenefullFromGzip.out.ch_multiqc)

        } else if (runMode == 'multiome') {
            // CHROMAP starts immediately (no dependency on ch_gzipped)
            ch_atac = prepareATACSamples()
            runChromapWorkflow(ch_atac)

            runMultiome(
                ch_gzipped,
                runChromapWorkflow.out.fragments,
                runChromapWorkflow.out.fragments_tbi,
                ch_star_logs
            )
            finalizeMultiqc(runMultiome.out.ch_multiqc.mix(runChromapWorkflow.out.mqc_json))

        } else {
            ch_counts_h5ad = createCountsH5ad(ch_samples)
            ADD_VELOCITY_LAYERS(ch_gzipped.join(ch_counts_h5ad))
            ch_velocity_signal = ADD_VELOCITY_LAYERS.out.velocity_h5ad.map { it[1] }
            finalizeMultiqc(ch_star_logs.mix(ch_velocity_signal))
        }

        runDebugAnalysis(finalizeMultiqc.out.report)
}

// ══════════════════════════════════════════════════════════════════════════
// Re-entry: resume from after SoupX (RNA QC + MultiQC only)
// For multiome: assumes CHROMAP + ArchR are already done externally.
// ══════════════════════════════════════════════════════════════════════════
workflow post_soupx_entry {
    main:
        context    = prepareSamples()
        ch_samples = context.ch_samples
        runMode    = context.runMode

        if (runMode != 'genefull' && runMode != 'multiome') {
            log.error "post_soupx_entry is only available in genefull or multiome run_mode."
            System.exit(1)
        }

        ch_star_logs     = starLogsFromPublish(ch_samples)
        scrubletChannels = scrubletOutputsFromPublish(ch_samples)
        soupxChannels    = soupxOutputsFromPublish(ch_samples)

        METAQC_MERGE(scrubletChannels.metaqc_partial.join(soupxChannels.rho))

        soupx_h5ad_output = params.bypass_soupX ? soupxChannels.pre_h5ad : soupxChannels.ambient_h5ad
        SAM_QC(soupx_h5ad_output.join(scrubletChannels.whitelist))

        ch_scrublet_plots = scrubletChannels.qc_plots.groupTuple().map { meta, plots ->
            def f = makeFinder(meta)
            [ meta,
              f(plots, '_knee_plot_mqc'),
              f(plots, '_doublet_score_histogram_QC1_mqc'),
              f(plots, '_violin_comparison_QC1_mqc'),
              f(plots, '_violin_mito_filtering_QC2_mqc') ]
        }

        ch_samqc_plots = SAM_QC.out.qc_plots.groupTuple().map { meta, plots ->
            def f = makeFinder(meta)
            [ meta,
              f(plots, '_umap_leiden_QC2_mqc'),
              f(plots, '_marker_genes_dotplot_mqc') ]
        }

        QC_PANEL(
            ch_scrublet_plots
                .join(soupxChannels.combined_plot)
                .join(ch_samqc_plots)
                .map { meta, knee, histogram, violin, mito, soupx_combined, umap, dotplot ->
                    [ meta, knee, histogram, violin, mito, soupx_combined, umap, dotplot ]
                }
        )

        ch_final_multiqc = ch_star_logs
            .mix(METAQC_MERGE.out.metaqc_table.map { it[1] })
            .mix(QC_PANEL.out.panel.map { it[1] })

        finalizeMultiqc(ch_final_multiqc)
        runDebugAnalysis(finalizeMultiqc.out.report)
}
