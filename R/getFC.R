#' Get fold-changes.
#'
#' This function is to use differential gene expression analysis to compute the fold-change in gene/feature expression
#' by comparing the cluster profile against all the other identified clusters together.
#'
#' @importFrom Seurat FindMarkers
#' @importFrom stats cutree
#'
#' @param obj A Seurat object with cluster.
#' @param min.pct A decimal value between 0 and 1. Only test genes that are detected in a minimum fraction of min.pct cells in either of the two populations. Meant to speed up the function by not testing genes that are very infrequently expressed.
#' @param logfc.threshold A decimal value between 0 and 1. Limit testing to genes which show, on average, at least X-fold difference (log-scale) between the two groups of cells. Increasing logfc.threshold speeds up the function, but can miss weaker signals.
#' @param test.use Denotes which test to use. Available options are 'wilcox', 'bimod', 'roc', 'negbinom', 'poisson', 'LR', 'MAST' and 'DESeq2'. The defalut is "MAST" for scRNAseq data, we suggest to use 'wilcox' for other data type.
#'
#' @return For each cluster, there will be a ranked gene/feature list.
#'
#' @export
#'
#' @examples
#' pbmc_example <- scqc(pbmc_small, min.cells = 1, min.features = 10, nfeatures = 100, npcs = 10)
#' pbmc_example <- doClustering(pbmc_example, dims = 1:10, k.param = 5, resolution = 0.75)
#' cluster_list <- getFC(pbmc_example, min.pct = 0.25, test.use = "MAST")
#'
getFC <- function(obj, min.pct = 0.25, test.use = "MAST", logfc.threshold = 0.1) {
  # check cluster
  if (is.null(obj@meta.data$seurat_clusters) == TRUE){
    stop(cat("Please go through 'check_cluster' fucntion to do cluster first.\n"))
  }

  # number of cluster
  ncluster <- length(unique(obj@meta.data$seurat_clusters))

  list.names <- paste0("Cluster", 0:(ncluster - 1))
  cluster_list <- vector("list", length(list.names))
  names(cluster_list) <- list.names

  # Get gene list with rank for each cluster
  for (i in 1:ncluster) {

    # find marker gene
    cat(paste0("Find marker genes for cluster"), i - 1, "\n")
    options(warn=-1)
    cluster.markers <- Seurat::FindMarkers(obj, ident.1 = (i - 1), min.pct = min.pct, logfc.threshold = logfc.threshold, test.use = test.use, min.cells.group = 1)

    # Order gene by 'avg_logFC'
    cluster.markers <- cluster.markers[order(cluster.markers[, "avg_logFC"], -cluster.markers[, "p_val_adj"], decreasing = TRUE), ]
    cluster.markers[which(cluster.markers[, "avg_logFC"] == Inf),  "avg_logFC"] = max(cluster.markers[-which(cluster.markers[, "avg_logFC"] == Inf), "avg_logFC"])
    cluster.markers[which(cluster.markers[, "avg_logFC"] == -Inf), "avg_logFC"] = min(cluster.markers[-which(cluster.markers[, "avg_logFC"] == -Inf), "avg_logFC"])
    cluster_list[[i]] <- cluster.markers[, "avg_logFC"]

    names(cluster_list[[i]]) <- rownames(cluster.markers)
    names(cluster_list[[i]]) <- toupper(names(cluster_list[[i]]))
  }
  # return cluster list
  return(cluster_list)
}
