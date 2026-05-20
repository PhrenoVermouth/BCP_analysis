#!/usr/bin/env Rscript
suppressPackageStartupMessages({
    library(ArchR)
    library(SoupX)
    library(SingleCellExperiment)
    library(zellkonverter)
    library(Seurat)
})

# Setup basic argument parsing (to avoid argparse dependency)
args_raw <- commandArgs(trailingOnly = TRUE)
args <- list()
for (i in seq(1, length(args_raw), by = 2)) {
    key <- gsub("--", "", args_raw[i])
    args[[key]] <- args_raw[i+1]
}

# Cast numerical args
args$tss_enrichment <- as.numeric(args$tss_enrichment)
args$min_frags <- as.numeric(args$min_frags)
args$rna_min_umi <- as.numeric(args$rna_min_umi)
args$threads <- as.numeric(args$threads)
if (is.null(args$genome)) args$genome <- "hg38"
if (is.null(args$barcode_translation)) args$barcode_translation <- "NULL"
if (is.null(args$sam_python)) args$sam_python <- "python"
if (is.null(args$project_dir)) args$project_dir <- "."

# Setup ArchR globals
addArchRThreads(threads = args$threads)
addArchRGenome(args$genome)
set.seed(1)

# QC tracking table: record cell counts and per-modality stats at each step
qc_track <- data.frame(
    Step=character(),
    Cells=integer(),
    Median_ATAC_Frags=numeric(),
    Median_ATAC_TSS=numeric(),
    Median_RNA_UMI=numeric(),
    Median_RNA_Genes=numeric(),
    stringsAsFactors=FALSE
)
track_step <- function(step_name, proj_obj, rna_available=FALSE) {
    n <- length(proj_obj$cellNames)
    med_frags <- median(proj_obj$nFrags, na.rm=TRUE)
    med_tss   <- round(median(proj_obj$TSSEnrichment, na.rm=TRUE), 2)
    if (rna_available && "Gex_nUMI" %in% colnames(getCellColData(proj_obj))) {
        med_umi   <- median(proj_obj$Gex_nUMI, na.rm=TRUE)
        med_genes <- median(proj_obj$Gex_nGenes, na.rm=TRUE)
    } else {
        med_umi   <- NA
        med_genes <- NA
    }
    qc_track[nrow(qc_track)+1, ] <<- list(step_name, as.integer(n),
                                           med_frags, med_tss, med_umi, med_genes)
}

out_archr_dir <- file.path(args$outdir, paste0(args$sample_id, "_ArchR_Project"))

# 1. Generate Arrow File from ATAC Fragments
cat("\n[1] Creating Arrow Files from ATAC fragments...\n")
arrow_files <- createArrowFiles(
    inputFiles = setNames(args$atac_fragments, args$sample_id),
    sampleNames = args$sample_id,
    filterTSS = args$tss_enrichment,
    filterFrags = args$min_frags,
    addTileMat = TRUE,
    addGeneScoreMat = TRUE,
    force = TRUE,
    logFile = createLogFile("createArrowFiles")
)

# 2. Build ArchR Project (ATAC Only)
cat("\n[2] Building ATAC ArchR Project and Initial Clustering...\n")
proj <- ArchRProject(
    ArrowFiles = arrow_files,
    outputDirectory = out_archr_dir,
    copyArrows = TRUE
)
track_step("After ATAC QC (TSS/Frags filter)", proj, rna_available=FALSE)

# ATAC Initial LSI & Clustering (crucial for SoupX)
proj <- addIterativeLSI(
    ArchRProj = proj,
    useMatrix = "TileMatrix", 
    name = "IterativeLSI_ATAC", 
    iterations = 2, 
    clusterParams = list(resolution = c(2), sampleCells = 10000, n.start = 10), 
    varFeatures = 25000, 
    dimsToUse = 1:30,
    force = TRUE
)
proj <- addClusters(
    input = proj,
    reducedDims = "IterativeLSI_ATAC",
    method = "Seurat",
    name = "ATAC_SeuratClusters",
    resolution = 1.0, 
    force = TRUE
)

# 3. Load RAW RNA and Translate Barcodes
cat("\n[3] Loading RNA raw matrix and preparing SoupX...\n")
rna_mtx_files <- list.files(args$rna_matrix_dir, pattern="matrix.mtx(\\.gz)?$", recursive=TRUE, full.names=TRUE)
# For SoupX, we STRICTLY want the 'raw' matrix, not the filtered one
raw_mtx_file <- rna_mtx_files[grepl("raw", rna_mtx_files)][1]
if(is.na(raw_mtx_file)) stop(paste("Could not find raw matrix.mtx(.gz) inside", args$rna_matrix_dir))

rna_mtx_file <- raw_mtx_file
rna_dir_exact <- dirname(rna_mtx_file)
rna_bc_file <- list.files(rna_dir_exact, pattern="barcodes.tsv(\\.gz)?$", full.names=TRUE)[1]
rna_feat_file <- list.files(rna_dir_exact, pattern="features.tsv(\\.gz)?$", full.names=TRUE)[1]

raw_rna_sparse <- Matrix::readMM(rna_mtx_file)
rna_bcs <- read.table(rna_bc_file, header=FALSE, stringsAsFactors=FALSE)$V1
# STARsolo features.tsv has 3 cols: Ensembl_ID, GeneSymbol, FeatureType
# ArchR's geneAnnotation uses gene symbols (V2), NOT Ensembl IDs (V1).
# Using V1 causes zero intersection downstream — always use V2 here.
rna_feats_tbl <- read.table(rna_feat_file, header=FALSE, stringsAsFactors=FALSE)
rna_feats <- if (ncol(rna_feats_tbl) >= 2) rna_feats_tbl$V2 else rna_feats_tbl$V1

# 10x references include ~20 duplicate gene symbols (e.g. TBCE) with distinct
# Ensembl IDs. SoupX errors on duplicate rownames, so sum-collapse them here.
if (any(duplicated(rna_feats))) {
    n_dup <- sum(duplicated(rna_feats))
    cat(sprintf("Collapsing %d duplicate gene symbols by summing counts...\n", n_dup))
    rowname_fac <- factor(rna_feats, levels = unique(rna_feats))
    agg <- Matrix::sparseMatrix(
        i = as.integer(rowname_fac),
        j = seq_along(rna_feats),
        x = 1,
        dims = c(length(levels(rowname_fac)), length(rna_feats))
    )
    raw_rna_sparse <- agg %*% raw_rna_sparse
    rna_feats <- levels(rowname_fac)
}

rownames(raw_rna_sparse) <- rna_feats
colnames(raw_rna_sparse) <- rna_bcs

# Barcode translation logic (if whitelist provided)
if (args$barcode_translation != "NULL" && file.exists(args$barcode_translation)) {
    cat("Translating RNA barcodes to ATAC barcodes based on 10x whitelist...\n")
    # For 10x Multiome, usually the whitelist is 10x ATAC barcode in col 1, RNA in col 2, or vice versa.
    # We load it assuming standard 10x format (could be comma separated)
    wl <- read.csv(args$barcode_translation, header=FALSE, stringsAsFactors=FALSE)
    # Assumed structure: ATAC, RNA, MULTIOME or just find matching columns.
    # We will build a named vector: rna_to_atac[ RNA_bc ] = ATAC_bc
    # Since whitelists vary, let's assume standard V2 is RNA, V1 is ATAC
    # Please adapt if custom whitelist mapping is used!
    # For many arcs, col 2 is RNA, col 1 is ATAC.
    if(ncol(wl) >= 2) {
        rna_to_atac <- setNames(wl$V1, wl$V2)
        # Apply translation. Keep only those that map
        mapped_bcs <- rna_to_atac[colnames(raw_rna_sparse)]
        valid_idx <- !is.na(mapped_bcs)
        raw_rna_sparse <- raw_rna_sparse[, valid_idx]
        colnames(raw_rna_sparse) <- mapped_bcs[valid_idx]
    }
}

# Align ArchR barcode format (SampleName#Barcode)
colnames(raw_rna_sparse) <- paste0(args$sample_id, "#", colnames(raw_rna_sparse))

# Auto-detect ATAC barcode orientation: Chromap may write fragments with RC'd
# barcodes depending on sequencer/i5 chemistry. Try forward first, RC if it
# yields more overlap with ArchR cell names.
rc_bc <- function(names_vec) {
    vapply(names_vec, function(x) {
        p <- strsplit(x, "#", fixed=TRUE)[[1]]
        bc_rc <- chartr("ACGT", "TGCA", paste(rev(strsplit(p[2], "")[[1]]), collapse=""))
        paste0(p[1], "#", bc_rc)
    }, character(1), USE.NAMES=FALSE)
}
n_fwd <- length(intersect(colnames(raw_rna_sparse), proj$cellNames))
rc_names <- rc_bc(colnames(raw_rna_sparse))
n_rc  <- length(intersect(rc_names, proj$cellNames))
cat(sprintf("Barcode orientation check: fwd overlap=%d, RC overlap=%d\n", n_fwd, n_rc))
if (n_rc > n_fwd) {
    cat("Using RC'd ATAC barcodes (Chromap fragments appear reverse-complemented).\n")
    colnames(raw_rna_sparse) <- rc_names
}

# Subset RNA to cells that survived ATAC QC
valid_cells <- intersect(colnames(raw_rna_sparse), proj$cellNames)
if(length(valid_cells) == 0){
    cat("\n--- DEBUG INFO ---\n")
    cat("Top 5 RNA matrix barcodes after translation:\n")
    print(head(colnames(raw_rna_sparse)))
    cat("Top 5 ATAC Arrow barcodes:\n")
    print(head(proj$cellNames))
    cat("------------------\n")
    stop("ERROR: Zero overlapping cells between RNA matrix and ATAC Arrow files! Check your Barcode formatting/translation.")
}
raw_rna_aligned <- raw_rna_sparse[, valid_cells, drop=FALSE]
# Temporarily subset proj to intersected cells for QC stats
proj_intersected <- proj[proj$cellNames %in% valid_cells, ]
track_step("After RNA-ATAC barcode intersection", proj_intersected, rna_available=FALSE)
rm(proj_intersected)

# 4. Cross-Modality SoupX Correction
cat("\n[4] Running ATAC-guided SoupX Correction...\n")
# Setup SoupX using full raw matrix for ambient profile, and the aligned valid cells for counts
soup <- SoupChannel(tod = raw_rna_sparse, toc = raw_rna_aligned)

# Enforce ATAC clusters onto SoupX
atac_clusters <- proj$ATAC_SeuratClusters[match(colnames(raw_rna_aligned), proj$cellNames)]
names(atac_clusters) <- colnames(raw_rna_aligned)

soup <- setClusters(soup, atac_clusters)
soup <- tryCatch({
    autoEstCont(soup)
}, error = function(e) {
    cat("\n[WARNING] SoupX autoEstCont failed: ", e$message, "\n")
    cat("This is common on small/test datasets. Falling back to a conservative 5% contamination rate.\n")
    setContaminationFraction(soup, 0.05)
})
adjSoup <- adjustCounts(soup, roundToInt=TRUE)

# Build SummarizedExperiment for ArchR
# We match rowRanges with ArchR's geneAnnotation
geneAnno <- getGeneAnnotation(proj)
# Map gene names
matched_genes <- intersect(rownames(adjSoup), geneAnno$genes$symbol)
adjSoup_matched <- adjSoup[matched_genes, ]
se_genes <- geneAnno$genes[match(rownames(adjSoup_matched), geneAnno$genes$symbol)]
sceSoup <- SummarizedExperiment(list(counts=adjSoup_matched), rowRanges=se_genes)

# 5. Integrate and Filter RNA
cat("\n[5] Integrating Cleaned RNA and Filtering UMI...\n")
proj <- addGeneExpressionMatrix(input=proj, seRNA=sceSoup, force=TRUE)
track_step("After SoupX + RNA integration", proj, rna_available=TRUE)
# RNA Min UMI Filter
proj <- proj[!is.na(proj$Gex_nUMI) & proj$Gex_nUMI >= args$rna_min_umi, ]
track_step(paste0("After RNA UMI filter (>=", args$rna_min_umi, ")"), proj, rna_available=TRUE)

# 6. Doublet Removal
cat("\n[6] Removing Trans-Modality Doublets...\n")
proj <- addDoubletScores(proj)
proj <- filterDoublets(proj)

track_step("After doublet removal", proj, rna_available=TRUE)

# 7. Joint LSI & Final Clustering
cat("\n[7] Executing Joint LSI and Final Clustering...\n")
# RNA LSI
proj <- addIterativeLSI(
    ArchRProj = proj,
    useMatrix = "GeneExpressionMatrix", 
    depthCol = "Gex_nUMI",
    name = "RNALSI",
    iterations = 2,
    clusterParams = list(resolution = c(2), sampleCells = 10000, n.start = 10),
    varFeatures = 2500,
    binarize = FALSE,
    firstSelection = "variable",
    force = TRUE
)

# Combined Dims
proj <- addCombinedDims(proj, reducedDims = c("IterativeLSI_ATAC", "RNALSI"), name = "JointLSI")
# UMAP and Clustering
proj <- addUMAP(proj, reducedDims = "JointLSI", name = "UMAP_Joint", minDist = 0.4, force=TRUE)
proj <- addClusters(
    input = proj,
    reducedDims = "JointLSI",
    method = "Seurat",
    name = "Joint_SeuratClusters",
    resolution = 0.8,
    force = TRUE
)

# 8. MACS2 Peak Calling
cat("\n[8] Calling Peaks via MACS2...\n")
has_peaks <- FALSE
tryCatch({
    proj <- addGroupCoverages(ArchRProj=proj, groupBy="Joint_SeuratClusters")
    proj <- addReproduciblePeakSet(
        ArchRProj = proj, 
        groupBy = "Joint_SeuratClusters", 
        pathToMacs2 = Sys.which("macs2")
    )
    proj <- addPeakMatrix(proj)
    proj <- addBgdPeaks(proj, method="ArchR")
    has_peaks <- TRUE
}, error = function(e) {
    cat("\n[WARNING] MACS2 peak calling failed:", e$message, "\n")
    cat("This is likely a MACS2 installation issue in the container.\n")
    cat("Peak calling skipped. You can rerun this step separately after fixing MACS2.\n")
})

track_step("Final (after clustering)", proj, rna_available=TRUE)

# 8b. Export RNA counts h5ad + clusters → call SAM marker script → read back → draw heatmap
cat("\n[8b] SAM-based marker calling and side-by-side heatmaps...\n")
tryCatch({
    # ── Export RNA counts h5ad (X = raw counts, no logCPM) ──
    cat("  Exporting RNA counts for SAM...\n")
    rna_se <- getMatrixFromProject(proj, useMatrix = "GeneExpressionMatrix")
    rna_sce <- SingleCellExperiment(list(counts = assay(rna_se)))
    rownames(rna_sce) <- rowData(rna_se)$name
    colData(rna_sce) <- getCellColData(proj, select = "Joint_SeuratClusters")
    sam_input_h5ad <- file.path(args$outdir, paste0(args$sample_id, "_rna_counts_for_sam.h5ad"))
    writeH5AD(rna_sce, sam_input_h5ad)

    # ── Export cluster assignment CSV ──
    cluster_csv <- file.path(args$outdir, paste0(args$sample_id, "_clusters.csv"))
    cl_data <- getCellColData(proj, select = "Joint_SeuratClusters")
    write.csv(data.frame(cluster = cl_data$Joint_SeuratClusters, row.names = rownames(cl_data)),
              cluster_csv, quote = FALSE)

    # ── Call SAM marker script (BCP2 conda python) ──
    sam_script <- file.path(args$project_dir, "bin", "run_sam_markers.py")
    sam_cmd <- sprintf(
        '%s %s --input %s --clusters %s --sample_id %s --outdir %s --top_n 5',
        args$sam_python, sam_script, sam_input_h5ad, cluster_csv, args$sample_id, args$outdir
    )
    cat("  Running SAM:", sam_cmd, "\n")
    ret <- system(sam_cmd)
    if (ret != 0) stop("SAM marker calling exited with code ", ret)

    # ── Read marker results ──
    marker_csv <- file.path(args$outdir, paste0(args$sample_id, "_sam_markers.csv"))
    marker_df <- read.csv(marker_csv, stringsAsFactors = FALSE)
    if (nrow(marker_df) == 0) {
        cat("[WARNING] SAM returned no markers. Skipping heatmap.\n")
    } else {
        # Get unique top genes preserving per-cluster order
        top_genes <- unique(marker_df$gene)
        cat("  SAM markers:", length(top_genes), "genes across",
            length(unique(marker_df$cluster)), "clusters\n")

        # ── Build side-by-side heatmap using ArchR matrices (still in memory) ──
        rna_mat <- assay(rna_se)
        rownames(rna_mat) <- rowData(rna_se)$name
        clusters <- getCellColData(proj, select = "Joint_SeuratClusters")[, 1]
        # Natural sort: C1, C2, ..., C10, C11 (not lexicographic C1, C10, C11, C2)
        raw_ids <- unique(clusters)
        nums <- as.numeric(gsub("^C", "", raw_ids))
        if (all(!is.na(nums))) {
            cluster_ids <- raw_ids[order(nums)]
        } else {
            cluster_ids <- sort(raw_ids)
        }

        gs_se <- getMatrixFromProject(proj, useMatrix = "GeneScoreMatrix")
        gs_mat <- assay(gs_se)
        rownames(gs_mat) <- rowData(gs_se)$name

        # Shared gene set (present in both matrices)
        shared_genes <- intersect(intersect(top_genes, rownames(rna_mat)), rownames(gs_mat))
        if (length(shared_genes) < 2) {
            cat("[WARNING] <2 shared genes between markers and GeneScore. Skipping heatmap.\n")
        } else {
            # Per-cluster mean z-scores
            rna_means <- sapply(cluster_ids, function(cl) {
                Matrix::rowMeans(rna_mat[shared_genes, clusters == cl, drop = FALSE])
            })
            rna_z <- t(scale(t(rna_means)))

            gs_means <- sapply(cluster_ids, function(cl) {
                Matrix::rowMeans(gs_mat[shared_genes, clusters == cl, drop = FALSE])
            })
            gs_z <- t(scale(t(gs_means)))

            # Row order: staircase by marker cluster assignment
            gene_order_df <- marker_df[marker_df$gene %in% shared_genes, ]
            gene_order_df <- gene_order_df[!duplicated(gene_order_df$gene), ]
            gene_order_df$cl_num <- as.numeric(gsub("^C", "", gene_order_df$cluster))
            gene_order_df <- gene_order_df[order(gene_order_df$cl_num, -gene_order_df$score), ]
            ordered_genes <- c(gene_order_df$gene, setdiff(shared_genes, gene_order_df$gene))

            rna_z <- rna_z[ordered_genes, ]
            gs_z  <- gs_z[ordered_genes, ]

            # Clamp z-scores
            rna_z[rna_z >  2] <-  2; rna_z[rna_z < -2] <- -2
            gs_z[gs_z  >  2] <-  2; gs_z[gs_z  < -2] <- -2

            col_fun <- circlize::colorRamp2(c(-2, 0, 2), c("#2166AC", "white", "#B2182B"))
            ht_rna <- ComplexHeatmap::Heatmap(
                rna_z,
                name = "RNA z-score", col = col_fun,
                cluster_rows = FALSE, cluster_columns = FALSE,
                column_title = "RNA Expression",
                row_names_gp = grid::gpar(fontsize = 7),
                column_names_gp = grid::gpar(fontsize = 9),
                show_row_names = TRUE, show_column_names = TRUE
            )
            ht_gs <- ComplexHeatmap::Heatmap(
                gs_z,
                name = "GeneScore z-score", col = col_fun,
                cluster_rows = FALSE, cluster_columns = FALSE,
                column_title = "Gene Activity Score",
                row_names_gp = grid::gpar(fontsize = 7),
                column_names_gp = grid::gpar(fontsize = 9),
                show_row_names = TRUE, show_column_names = TRUE
            )

            ht_list <- ht_rna + ht_gs
            png(file.path(args$outdir, paste0("3.", args$sample_id, "_Marker_Heatmaps_mqc.png")),
                width = 14, height = max(6, length(shared_genes) * 0.18 + 2),
                units = "in", res = 300)
            ComplexHeatmap::draw(ht_list,
                column_title = paste0(args$sample_id, " \u2014 SAM Markers: RNA Expression vs Gene Activity Score"),
                column_title_gp = grid::gpar(fontsize = 14, fontface = "bold"))
            dev.off()
            cat("  Heatmap exported (", length(shared_genes), "genes).\n")
        }
    }
}, error = function(e) {
    cat("[WARNING] SAM marker calling/heatmap failed:", e$message, "\n")
    cat("Continuing without marker heatmaps.\n")
})

# Write QC tracking table as MultiQC-compatible TSV
qc_tsv_path <- file.path(args$outdir, paste0(args$sample_id, "_multiome_qc_mqc.tsv"))
qc_con <- file(qc_tsv_path, "w")
writeLines("# plot_type: 'table'", qc_con)
writeLines("# section_name: 'Multiome QC Tracking'", qc_con)
writeLines("# description: 'Cell counts at each filtering step in multiome integration'", qc_con)
writeLines("# pconfig:", qc_con)
writeLines("#     sortRows: false", qc_con)
close(qc_con)
write.table(qc_track, qc_tsv_path, sep="\t", quote=FALSE, row.names=FALSE, append=TRUE)
cat("QC tracking table written.\n")

# 9. Outputs and Export
cat("\n[9] Exporting diagnostic plots and H5AD matrices...\n")
dir.create(file.path(args$outdir, "plots"), showWarnings=FALSE, recursive=TRUE)

# A. TSS vs nFrags
png(file.path(args$outdir, paste0("0.", args$sample_id, "_TSS_vs_Frags_mqc.png")), width=6, height=6, units="in", res=300)
df <- getCellColData(proj, select = c("log10(nFrags)", "TSSEnrichment"))
p1 <- ggPoint(x=df[,1], y=df[,2], colorDensity=TRUE, continuousSet="sambaNight",
              xlabel="Log10 Unique Fragments", ylabel="TSS Enrichment") +
      geom_hline(yintercept=args$tss_enrichment, lty="dashed") +
      geom_vline(xintercept=log10(args$min_frags), lty="dashed")
print(p1)
dev.off()

# B. Fragment Size
png(file.path(args$outdir, paste0("1.", args$sample_id, "_Fragment_Sizes_mqc.png")), width=6, height=6, units="in", res=300)
p2 <- plotFragmentSizes(ArchRProj = proj)
print(p2)
dev.off()

# C. Joint UMAP
png(file.path(args$outdir, paste0("2.", args$sample_id, "_Joint_UMAP_mqc.png")), width=8, height=8, units="in", res=300)
p3 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Joint_SeuratClusters", embedding = "UMAP_Joint")
print(p3)
dev.off()

# D. Export Gene Scores & Peaks to H5AD (with RDS fallback if anndata is missing)
export_matrix <- function(proj_obj, matrix_name, out_base) {
    mat <- getMatrixFromProject(proj_obj, useMatrix=matrix_name)
    sce <- SingleCellExperiment(list(counts=assay(mat)))
    feat_names <- rowData(mat)$name
    if (is.null(feat_names) || length(feat_names) == 0) {
        rr <- rowRanges(mat)
        feat_names <- paste0(as.character(seqnames(rr)), ":", start(rr), "-", end(rr))
    }
    rownames(sce) <- feat_names
    colData(sce) <- getCellColData(proj_obj, select='Joint_SeuratClusters')
    tryCatch({
        # zellkonverter::writeH5AD() accepts a SingleCellExperiment directly.
        # Do NOT pass an AnnData object returned by SCE2AnnData() — that would fail.
        writeH5AD(sce, paste0(out_base, ".h5ad"))
        cat("Exported", matrix_name, "as H5AD\n")
    }, error = function(e) {
        cat("[WARNING] H5AD export failed:", e$message, "\n")
        cat("Falling back to RDS format...\n")
        saveRDS(sce, paste0(out_base, ".rds"))
    })
}

export_matrix(proj, "GeneScoreMatrix", file.path(args$outdir, paste0(args$sample_id, "_genescores")))
if (has_peaks) {
    export_matrix(proj, "PeakMatrix", file.path(args$outdir, paste0(args$sample_id, "_peaks")))
    peaks <- getPeakSet(proj)
    write.csv(as.data.frame(unname(peaks)), file.path(args$outdir, paste0(args$sample_id, "_peaks.csv")), quote=FALSE, row.names=FALSE)
} else {
    cat("Skipping PeakMatrix/Peak export (MACS2 failed).\n")
    write.csv(data.frame(note="Peak calling skipped due to MACS2 error"), file.path(args$outdir, paste0(args$sample_id, "_peaks.csv")), quote=FALSE, row.names=FALSE)
}

cat("\n[10] Saving final ArchR Project...\n")
saveArchRProject(ArchRProj = proj, load = FALSE)
cat("Integration complete.\n")
