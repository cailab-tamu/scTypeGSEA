#' Get full changes.
#'
#' This function is to do DE for each cluster regarding all other clusters.
#' And then for each cluster, it will return a ranked gene list.
#'
#' @importFrom Seurat FindMarkers
#' @param obj A Seurat object with cluster.
#' @param min.pct only test genes that are detected in a minimum fraction of min.pct cells in either of the two populations. Meant to speed up the function by not testing genes that are very infrequently expressed.
#' @param logfc.threshold Limit testing to genes which show, on average, at least X-fold difference (log-scale) between the two groups of cells. Increasing logfc.threshold speeds up the function, but can miss weaker signals.
#' @param test.use Denotes which test to use. Available options are 'wilcox', 'bimod', 'roc', 'negbinom', 'poisson', 'LR', 'MAST' and 'DESeq2'.
#'
#' @return For each cluster, there will be a ranked gene list.
#' @export
#' @examples
#' pbmc_example <- doClustering(pbmc_test, nfeatures = 100, npcs = 10,
#'                               dims = 1:10, k.param = 5, resolution = 0.75)
#' cluster_list <- getFC(pbmc_example, min.pct = 0.25, test.use = "MAST")
#'
getFC <- function(obj, min.pct = 0.1, test.use = "wilcox", logfc.threshold = 0.25) {

  # number of cluster
  ncluster <- length(unique(obj@meta.data$seurat_clusters))
  if (ncluster == 0) stop(paste0("You should do cluster first, please go through 'check_cluster' fucntion."))

  list.names <- paste0("Cluster", 0:(ncluster - 1))
  cluster_list <- vector("list", length(list.names))
  names(cluster_list) <- list.names

  # Get gene list with rank for each cluster
  for (i in 1:ncluster) {

    # find marker gene
    cat(paste0("Find marker genes for cluster"), i - 1, "\n")
    options(warn=-1)
    cluster.markers <- Seurat::FindMarkers(obj, ident.1 = (i - 1), min.pct = min.pct, logfc.threshold = logfc.threshold, test.use = test.use)

    # Order gene by 'avg_logFC'
    cluster.markers <- cluster.markers[order(cluster.markers[, "avg_logFC"], -cluster.markers[, "p_val_adj"], decreasing = TRUE), ]
    cluster_list[[i]] <- cluster.markers[, "avg_logFC"]
    names(cluster_list[[i]]) <- rownames(cluster.markers)
    names(cluster_list[[i]]) <- toupper(names(cluster_list[[i]]))
  }

  return(cluster_list)
}
