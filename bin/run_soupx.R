#!/usr/bin/env Rscript

if (Sys.getenv("BASILISK_USE_SYSTEM_BINARIES") == "") {
  Sys.setenv(BASILISK_USE_SYSTEM_BINARIES = "TRUE")
}

python_bin <- Sys.which("python")
if (Sys.getenv("RETICULATE_PYTHON") == "" && python_bin != "") {
  Sys.setenv(RETICULATE_PYTHON = python_bin)
}

# 1. library loading
suppressPackageStartupMessages(library(SoupX))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(DropletUtils))
suppressPackageStartupMessages(library(zellkonverter))
suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(ggplot2))

# 2. parse paras
parser <- ArgumentParser(description="Run SoupX to correct ambient RNA")
parser$add_argument("--raw_dir", required=TRUE, help="Path to the raw matrix directory")
parser$add_argument("--filtered_dir", required=TRUE, help="Path to the filtered matrix directory")
parser$add_argument("--sample_id", required=TRUE, help="Sample identifier")
parser$add_argument("--whitelist", required=TRUE, help="Path to whitelist barcodes")
args <- parser$parse_args()

# 3. load input
toc <- Read10X(data.dir = args$raw_dir)
tod <- Read10X(data.dir = args$filtered_dir)

whitelist_barcode <- readLines(args$whitelist, warn = FALSE)
whitelist_barcode <- whitelist_barcode[whitelist_barcode %in% colnames(tod)]
if (length(whitelist_barcode) == 0) {
  stop("No whitelist barcodes found in filtered matrix.")
}
tod <- tod[, whitelist_barcode, drop = FALSE]
sc <- SoupChannel(toc, tod)

# 4. clustering
srat_tmp <- CreateSeuratObject(counts = tod)
srat_tmp <- NormalizeData(srat_tmp, verbose = FALSE)
srat_tmp <- FindVariableFeatures(srat_tmp, verbose = FALSE)
srat_tmp <- ScaleData(srat_tmp, verbose = FALSE)
srat_tmp <- RunPCA(srat_tmp, verbose = FALSE)
srat_tmp <- FindNeighbors(srat_tmp, dims = 1:20, verbose = FALSE)
srat_tmp <- FindClusters(srat_tmp, verbose = FALSE)

seurat_clusters <- srat_tmp@meta.data$seurat_clusters
names(seurat_clusters) <- rownames(srat_tmp@meta.data)
sc <- setClusters(sc, clusters = seurat_clusters)

# 5. auto-correction
# Capture contamination estimation plot
contamination_plot_file <- paste0('0.',args$sample_id, "_soupx_contamination_estimation_mqc.png")
png(contamination_plot_file, width = 600, height = 400)
sc <- autoEstCont(sc, forceAccept = TRUE)
dev.off()
#Hard code 260101
if (args$sample_id == 'q_em') {
    sc <- setContaminationFraction(sc, 0.2)
    warning("Hard-coded for sample ", args$sample_id, "Fraction = 0.2")
} 


cat("Contamination estimation plot saved to:", contamination_plot_file, "\n")

# write uncorrected h5ad for comparison
pre_soupx_file <- paste0(args$sample_id, "_pre_soupx.h5ad")
raw_counts_mat <- as(tod, "dgCMatrix")
sce_pre <- SingleCellExperiment(assays = list(counts = raw_counts_mat))
writeH5AD(sce_pre, file = pre_soupx_file)
cat("Uncorrected matrix saved to:", pre_soupx_file, "\n")


adj_counts <- as(adjustCounts(sc), "dgCMatrix")

raw_counts <- Matrix::colSums(raw_counts_mat)
corr_counts <- Matrix::colSums(adj_counts)
removed_fraction <- (raw_counts - corr_counts) / raw_counts
plot_df <- data.frame(fraction_removed = removed_fraction)

# Create ambient RNA removal plot as ggplot object
p_ambient <- ggplot(plot_df, aes(x = fraction_removed)) +
  geom_histogram(binwidth = 0.01, fill = "steelblue", color = "black") +
  labs(title = "Fraction of Counts Removed by SoupX",
       x = "Fraction Removed",
       y = "Number of Cells")

# Save individual ambient RNA removal plot
plot_file <- paste0("0.",args$sample_id, "_ambient_RNA_removed_mqc.png")
ggsave(plot_file, plot = p_ambient, width = 6, height = 4)
cat("Ambient RNA removal plot saved to:", plot_file, "\n")

# Create combined side-by-side plot using gridExtra
suppressPackageStartupMessages(library(gridExtra))
suppressPackageStartupMessages(library(png))
suppressPackageStartupMessages(library(grid))

# Read the contamination estimation plot as an image
img_contam <- readPNG(contamination_plot_file)
g_contam <- rasterGrob(img_contam, interpolate = TRUE)

# Combine both plots side by side
combined_plot_file <- paste0("0.",args$sample_id, "_soupx_combined_mqc.png")
png(combined_plot_file, width = 1200, height = 400)
grid.arrange(g_contam, ggplotGrob(p_ambient), ncol = 2)
dev.off()
cat("Combined SoupX plot saved to:", combined_plot_file, "\n")


# 6. output .h5ad and rho
output_file <- paste0(args$sample_id, "_rm_ambient.h5ad")
sce_to_write <- SingleCellExperiment(assays = list(counts = adj_counts))
writeH5AD(sce_to_write, file = output_file)


rho_file <- paste0(args$sample_id, "_soupx_rho.tsv")
rho_df <- data.frame(Sample = args$sample_id, Rho = sc$metaData$rho)
write.table(rho_df, file = rho_file, sep = '\t', row.names = FALSE, quote = FALSE)

cat("SoupX correction complete. Corrected matrix saved to:", output_file, "\n")
cat("Ambient-removed matrix saved to:", output_file, "\n")
cat("SoupX rho saved to:", rho_file, "\n")
