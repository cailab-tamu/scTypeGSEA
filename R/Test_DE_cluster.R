#' Title Gene rank of identity cluster
#'
#' This function is to do DE for each cluster regarding all other clusters.
#' And then for each cluster, it will return a ranked gene list.
#'
#' @importFrom Seurat FindMarkers
#' @param obj A Seurat object with cluster.
#' @param min.pct only test genes that are detected in a minimum fraction of min.pct cells in either of the two populations. Meant to speed up the function by not testing genes that are very infrequently expressed.
#' @param test.use Denotes which test to use. Available options are "wilcox", "bimod", "roc", "negbinom", "poisson", "LR", "MAST" and "DESeq2".
#'
#' @return For each cluster, there will be a ranked gene list.
#' @export
#'
Test_DE_cluster <- function(obj, min.pct = 0.25, test.use = "wilcox"){

  # number of cluster
  ncluster <- length(unique(obj@meta.data$seurat_clusters))

  list.names <- paste0("Cluster", 0:(ncluster - 1))
  cluster_list <- vector("list", length(list.names))
  names(cluster_list) <- list.names

  # Get gene list with rank for each cluster
  for (i in 1:ncluster){

    # find marker gene
    cat(paste0("Find marker genes for cluster"), i - 1)
    cluster.markers <- Seurat::FindMarkers(obj, ident.1 = (i-1), min.pct = min.pct, test.use = test.use)

    # Order gene by "avg_logFC"
    cluster.markers <- cluster.markers[order(cluster.markers[, "avg_logFC"], -cluster.markers[, "p_val_adj"], decreasing = TRUE), ]
    cluster_list[[i]] <- cluster.markers[, "avg_logFC"]
    names(cluster_list[[i]]) <- rownames(cluster.markers)
  }

  return(cluster_list)
}
