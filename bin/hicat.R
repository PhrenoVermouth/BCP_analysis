#require libraries
library(parallel) # parallelization
library(scrattch.hicat) #hicat
library(dendextend) #dendrogram
library(dplyr) # %in%
library(matrixStats) #matrix operation
library(Matrix) #matrix operation
library(Seurat) # Seurat single cell tool kit.
library(ggplot2) # ggplot
library(RColorBrewer) #for making color palette
library(stringr)#change gene name


# Function predict_remaining_cell_cluster
# Input file, file name of the .rda file produced during bootstraping. e.g., "subsample_PCA/result.2.rda"
# Output: a vector of clusters. factor name is the cell barcode. example as bellow:
# AAACCTGAGACCCACC AAACCTGAGCTAAGAT AAACCTGAGGGTCGAT AAACCTGGTTCTGGTA AAACGGGAGCCTATGT AAACGGGAGCTTCGCG 
#             "1"              "5"              "7"              "7"              "5"              "7" 
predict_remaining_cell_cluster <- function(file, all.cells, norm.dat){
  print(file)
  tmp=load(file)
  if(is.null(result)){
    return(NULL)
  }
  cl= result$cl
  test.cells = setdiff(all.cells, names(cl))
  markers=unique(result$markers)
  map.df = map_by_cor(norm.dat[markers,names(cl)],cl, norm.dat[markers,test.cells],method="mean")$pred.df
  test.cl = setNames(map.df$pred.cl, row.names(map.df))
  all.cl = c(setNames(as.character(cl),names(cl)), setNames(as.character(test.cl), names(test.cl)))
  return(all.cl[all.cells])
}



# Function, cocluster_matrix
# Input: a cluster object, vector of cluster name, factor name is cell barcode
# Output: a cocluster matrix

cocluster_matrix<- function(cl) {
  n.cells <- length(cl)
  co.matrix= matrix(0, nrow=n.cells, ncol=n.cells)
  row.names(co.matrix)=colnames(co.matrix)=names(cl)
  clusters <- unique(cl)
  for (i in clusters) {
    cells_i <- names(cl[cl==i])
    co.matrix[cells_i, cells_i]<-1
    
  }
  return(co.matrix)
}


# Function collect_subsample_cl_matrix_yw
# rewrite collect_subsample_cl_matrix. it runs into an error
# norm.dat is the log2p1 result of raw count, or SCT
# result.files is a vector of all bootstraping results as .rda, e.g., c(round1_result.rda, round2_result.rda)
# all.cells is the name of all cell barcodes in the norm.dat.
# max.cl.size is the maximal number of cells in a cluster. for ploting purpose, for example, a big cluster will obscure the visualization. has to be between 1-1000
# mc.cores is the cores used for parallelization.
# parallelization change to forking. WARNING, this makes the code only run on linux system.
# niter is the number of iteration.
# Output, a list containing cl.list and cl.mat
# cl.list is the cluster labels.
# cl.mat is the cocluster matrix.


infer_concensus_cluster <- function(norm.dat,result.files,all.cells,max.cl.size=NULL,mc.cores=1) {
  #pre-test
  #Do these files exist?
  if (!(TRUE %in% file.exists(result.files))) {
    print("Warning: bootstrap result files don't exist.")
  }
  #Are there more than 1 files in result.file?
  if (!length(result.files>1)) {
    print("Warning: bootstrap result files less than 2.")
  }
  #is norm.dat a sparce matrix?
  if (!is(SCT.dat, 'sparseMatrix')) {
    print("Warning: data is not sparce.")
  }
  #Do all.cells match norm.dat?
  if (!(TRUE %in% c(all.cells==colnames(SCT.dat)))) {
    print("Warning: cell names are inconsistant between inputs.")
  }
  #is mc.cores beyond the machine limit?
  if (detectCores()<mc.cores) {
    print("Warning: less cores than expected.")
  }
  #is max.cl.size within 1-1000 range?
  if (!is.null(max.cl.size)) {
    if (max.cl.size<0 | max.cl.size >1000) {
      print("Warning: maximal cells in cluster out of range")
    }
  }
  #step 1.
  #predict remaining 20% cells clusters on all iterations. 
  #output a data frame of the cells's cluster at each iteration
  if (mc.cores==1){ #single core behavior
    cl.list = lapply(result.files, predict_remaining_cell_cluster, all.cells = all.cells, norm.dat = norm.dat)
  }
  # multi-core behavior, forking.
  cl.list = mclapply(result.files, predict_remaining_cell_cluster, all.cells = all.cells, norm.dat = norm.dat, mc.cores = mc.cores)
  #step 2.
  # generate per iteration coclstering matrix.
  #if (mc.cores==1) #single core behavior
  cl.mat = lapply(cl.list,  cocluster_matrix)
  #multi core behavior
  #cl.mat = mclapply(cl.list,  cocluster_matrix, mc.cores=mc.cores)
  #step 3.
  #calculate mean frequency
  cl.mat.mean <- Reduce("+", cl.mat) / length(cl.mat)  
  #ste 4.
  #infer concensus cluster
  #cluster cells based on coclustering matrix
  nclust<-sapply(cl.list, unique)
  avg_nclust = mean(apply(nclust,2,length))
  k <- avg_nclust
  set.seed(1234)  # for reproducibility
  cl.cut <- kmeans(cl.mat.mean, centers = k)
  #step 5.
  #return cluster label and coclustering matrix
  return(list(cl=cl.cut$cluster, cl.mat = cl.mat.mean))
  
}

# Function bootstrap_iter_clust
# Input n_iter, number of iteration,
# Input norm.dat, normalized data, log2cpm
# Input de.param, de threshold
# Input directory, folder to save results
# Input mc.cores, how many cores to use
# Output save each round of iteration results as .rda file in folder ./directory/

bootstrap_iter_clust <- function(n_iter, norm.dat, de.param, directory, mc.cores){
  set.seed(12345)
  time1 = Sys.time()
  saveRDS(time1, file = "startime.rds")
  result <- mclapply(1:n_iter, function(i) {
    prefix <- paste0("iter",i)    
    tmp.cells <- sample(all.cells, round(0.8 * length(all.cells)))
    result <- iter_clust(norm.dat, 
                         select.cells = tmp.cells, 
                         prefix = prefix, 
                         dim.method = "pca", 
                         de.param = de.param)
    save(result, file = file.path(directory, paste0("result.", i, ".rda")))
  },
  mc.cores = mc.cores
  )
  time2 = Sys.time()
  saveRDS(time2, file = paste("endtime.rds",sep = ""))
}



# function iter_clust_hicat
# a wrapper around hicat iter_clust function to run with Seurat object and SCT normalization method.
# scdata is a Seurat object with 10X, umi count matrix normalized with SCT.
# de.param is a vector following hicat format
# reduction is between "pca", "umap", "tsne"
iter_clust_hicat <- function(scdata, de.param, reduction) {
  # Convert SCT count to matrix
  SCT.dat <- scdata@assays$SCT$scale.data
  SCT.dat <- Matrix(SCT.dat, sparse = TRUE)
  
  # Setting clustering parameters
  de.param <- de.param
  
  # Dimension Filtering
  # Remove pc that is correlated with depth
  gene.counts <- colSums(SCT.dat > 0)
  rm.eigen <- matrix(gene.counts, ncol = 1)
  row.names(rm.eigen) <- names(gene.counts)
  colnames(rm.eigen) <- "SCT"
  
  # Iterative clustering
  # call iter.result from hicat
  time<-Sys.time()
  print(time)
  iter.result <- iter_clust(SCT.dat, 
                            dim.method = "pca", 
                            de.param = de.param, 
                            rm.eigen = rm.eigen)
  time<-Sys.time()
  print(time)
  # merge clusters with less DE genes
  rd.dat <- t(SCT.dat[iter.result$markers,])
  merge.result <- merge_cl(SCT.dat, 
                           cl = iter.result$cl, 
                           rd.dat = rd.dat,
                           de.param = de.param)
  # save number of clusters
  nclust.hicat.high <- length(unique(iter.result$cl))
  nclust.hicat.merged <- length(unique(merge.result$cl))
  
  # Add cluster label to metadata
  hicat_merged <- merge.result$cl
  scdata.1<-scdata
  scdata.1@meta.data<-cbind(scdata.1@meta.data, hicat_merged)
  Idents(scdata.1) <- "hicat_merged"
  n_col <- length(unique(scdata.1@meta.data$"hicat_merged"))
  color2<-colorRampPalette(rev(brewer.pal(n = 8, name ="Set1")))(n_col)
  plot <- DimPlot(scdata.1, cols = sample(color2), reduction=reduction, group.by = "hicat_merged")+coord_fixed()+NoLegend()
  
  #return output
  cluster.result<-list(scdata.1, iter.result, plot)
  return(cluster.result)
}







