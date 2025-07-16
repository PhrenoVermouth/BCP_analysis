

#required libraries
library(Seurat)

# ImportH5ad
# import data to Seurat format
# Input: folder. Folder containing the processed files from h5ad. Should include a .mtx sparce matrix count table, cellID, geneID, metadata. 
# org: organism code.
# Output: a seurat object.
# note: in the output file, the .mtx is loaded to RNA assay, count slot. Optimally it is umi count or raw counts.

ImportH5ad <- function (folder, org) {
  # load cell barcode and gene ID
  cellID <- read.csv(file = paste0(folder, "cellID_", org, ".csv"), header=T, row.names = 1)
  geneID <- read.csv(file = paste0(folder, "geneID_", org, ".csv"), header=T, row.names = 1)
  #load metadata
  metadata <- read.csv(file = paste0(folder, "metadata_",org, ".csv"), header=T, row.names = 1)
  #load count table
  counts <- readMM(file = paste0(folder, "scdata_", org, ".mtx"))
  rownames(counts) <- cellID[,1]
  colnames(counts) <- geneID[,1]
  counts<-t(counts)
  scdata <- CreateSeuratObject(counts = counts, assay = "RNA", meta.data = metadata, project = paste(org))
  scdata <- DietSeurat(scdata)
}


