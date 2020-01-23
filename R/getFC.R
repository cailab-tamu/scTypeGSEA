#' Get fold-changes.
#'
#' This function is to use differential gene expression analysis to compute the fold-change in gene/feature expression
#' by comparing the cluster profile against all the other identified clusters together.
#'
#' @importFrom Seurat FindMarkers
#' @importFrom stats cutree
#'
#' @param obj A Seurat object with cluster or a list including data and results from hierarchical clustering.
#' @param min.pct A decimal value between 0 and 1. Only test genes that are detected in a minimum fraction of min.pct cells in either of the two populations. Meant to speed up the function by not testing genes that are very infrequently expressed.
#' @param logfc.threshold A decimal value between 0 and 1. Limit testing to genes which show, on average, at least X-fold difference (log-scale) between the two groups of cells. Increasing logfc.threshold speeds up the function, but can miss weaker signals.
#' @param test.use Denotes which test to use. Available options are 'wilcox', 'bimod', 'roc', 'negbinom', 'poisson', 'LR', 'MAST' and 'DESeq2'.
#' @param ncluster An integer, which is the number of cluster when your input including results from hierarchical clustering.
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
getFC <- function(obj, datatype = "sc", min.pct = 0.25, test.use = "MAST", logfc.threshold = 0.1,
                  ncluster = 3) {

  if (datatype == "sc"){
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
      cluster.markers <- Seurat::FindMarkers(obj, ident.1 = (i - 1), min.pct = min.pct, logfc.threshold = logfc.threshold, test.use = test.use)

      # Order gene by 'avg_logFC'
      cluster.markers <- cluster.markers[order(cluster.markers[, "avg_logFC"], -cluster.markers[, "p_val_adj"], decreasing = TRUE), ]
      cluster_list[[i]] <- cluster.markers[, "avg_logFC"]
      names(cluster_list[[i]]) <- rownames(cluster.markers)
      names(cluster_list[[i]]) <- toupper(names(cluster_list[[i]]))
    }
    # return cluster list
    return(cluster_list)
  } else{
    cluster_membership <- cutree(obj$cluster_cell, k = ncluster)

    list.names <- paste0("Cluster", 0:(ncluster - 1))
    cluster_list <- vector("list", length(list.names))
    names(cluster_list) <- list.names

    nfeatures <- nrow(obj[[1]])

    # Get feature list with rank for each cluster
    for (i in 1:ncluster) {
      # calculate fold-change
      index_cell <- which(cluster_membership == i)
      cell1 <- obj$data[, index_cell, drop = F]
      cell2 <- obj$data[, -index_cell, drop = F]


      wilcox_res <- matrix(0, nrow = nfeatures, ncol = 2)
      colnames(wilcox_res) <- c("statistic", "p-value")
      rownames(wilcox_res) <- rownames(obj[[1]])

      for (j in 1:nfeatures){
        test_res <- wilcox.test(cell1[j, ], cell2[j, ])
        wilcox_res[j, 1] <- test_res$statistic
        wilcox_res[j, 2] <- test_res$p.value
      }

      wilcox_res <- wilcox_res[order(wilcox_res[, 1], wilcox_res[, 2], decreasing = FALSE), ]

      # return feature list
      cluster_list[[i]] <- wilcox_res[, 1]
      names(cluster_list[[i]]) <- rownames(wilcox_res)
    }
    # return cluster list
    return(cluster_list)
  }
}
