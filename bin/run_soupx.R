#!/usr/bin/env Rscript

# 1. library loading
suppressPackageStartupMessages(library(SoupX))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(DropletUtils))
suppressPackageStartupMessages(library(zellkonverter))
suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(argparse))

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
sc <- autoEstCont(sc)
adj_counts <- adjustCounts(sc)

# 6. output .h5ad
output_file <- paste0(args$sample_id, "_corrected.h5ad")
sce_to_write <- SingleCellExperiment(assays = list(counts = adj_counts))
writeH5AD(sce_to_write, file = output_file)

cat("SoupX correction complete. Corrected matrix saved to:", output_file, "\n")
