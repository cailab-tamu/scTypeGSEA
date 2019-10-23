
Test_DE_cluster <- function(obj){
  ncluster <- length(unique(obj@meta.data$seurat_clusters))
  list.names <- paste0("Cluster", 0:(ncluster - 1))
  cluster_list <- vector("list", length(list.names))
  names(cluster_list) <- list.names

  for (i in 1:ncluster){
    cluster.markers <- FindMarkers(obj, ident.1 = (i-1), min.pct = 0.25)
    index_marker <- which(cluster.markers[, "p_val_adj"] < 0.05)
    cluster.markers <- cluster.markers[index_marker, ]
    cluster.markers <- cluster.markers[order(cluster.markers[, "avg_logFC"], decreasing = TRUE), ]
    cluster_list[[i]] <- cluster.markers[, "avg_logFC"]
    names(cluster_list[[i]]) <- rownames(cluster.markers)
  }

  return(cluster_list)
}
