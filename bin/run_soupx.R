#!/usr/bin/env Rscript

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
args <- parser$parse_args()

# 3. load input
toc <- Read10X(data.dir = args$raw_dir)
tod <- Read10X(data.dir = args$filtered_dir)
sc <- SoupChannel(toc, tod)

# 4. clustering
srat_tmp <- CreateSeuratObject(counts = sc$toc)
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
contamination_plot_file <- paste0('0.',args$sample_id, "_soupx_contamination_estimation_mqc.png")
png(contamination_plot_file, width = 600, height = 400)
sc <- autoEstCont(sc)
dev.off()
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
plot_file <- paste0("0.",args$sample_id, "_ambient_RNA_removed_mqc.png")
p <- ggplot(plot_df, aes(x = fraction_removed)) +
  geom_histogram(binwidth = 0.01, fill = "steelblue", color = "black") +
  labs(title = "Fraction of Counts Removed by SoupX",
       x = "Fraction Removed",
       y = "Number of Cells")
ggsave(plot_file, plot = p, width = 6, height = 4)
cat("Ambient RNA removal plot saved to:", plot_file, "\n")


# 6. output .h5ad and rho
output_file <- paste0(args$sample_id, "_corrected.h5ad")
sce_to_write <- SingleCellExperiment(assays = list(counts = adj_counts))
writeH5AD(sce_to_write, file = output_file)


rho_file <- paste0(args$sample_id, "_soupx_rho.tsv")
rho_df <- data.frame(Sample=args$sample_id, Rho=sc$rho)
write.table(rho_df, file = rho_file, sep='\t', row.names = FALSE, quote = FALSE)

cat("SoupX correction complete. Corrected matrix saved to:", output_file, "\n")
cat("Ambient-removed matrix saved to:", ambient_output_file, "\n")
cat("SoupX rho saved to:", rho_file, "\n")
