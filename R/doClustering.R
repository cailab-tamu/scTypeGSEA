#' Cluster Determination
#'
#' It will do cluster for any type of data. If data type is single cell, the input must be Seurat object and it will use "Findcluster" function in Seurat pacakge.
#' For any other data types, it will do hierarchical clustering.
#'
#' @importFrom Seurat Project
#' @importFrom Seurat FindNeighbors
#' @importFrom Seurat FindClusters
#' @importFrom stats dist
#' @importFrom stats hclust
#'
#' @param obj A Seurat object or any matrix where each column is a cell.
#' @param datatype Data type to do cluster, which can be "sc" for single cell data and "others" for all other data type.
#' @param dims An integer value. Define dimensions of reduction to use as input. (Do cluster for single cell data.)
#' @param k.param An integer value. Defines k for the k-nearest neighbor algorithm. (Do cluster for single cell data.)
#' @param resolution Value of the resolution parameter, use a value above (below) 1.0 if you want to obtain a larger (smaller) number of communities. (Do cluster for single cell data.)
#' @param method The agglomeration method to be used for hierarchical clustering, defalut is "complete".
#'
#' @return For single cell type of data, it will return a Seurat object with cluster. For other type of data, it will return a list including data and result from hierarchical clustering.
#'
#' @export
#'
#' @examples
#' pbmc_example <- scqc(pbmc_small, min.cells = 1, min.features = 10, nfeatures = 100, npcs = 10)
#' is.null(pbmc_example@meta.data$seurat_clusters)
#' pbmc_example <- doClustering(pbmc_example, dims = 1:10, k.param = 5, resolution = 0.75)
#' head(pbmc_example@meta.data$seurat_clusters)
doClustering <- function(obj, datatype = "sc", cluster_list = NULL, dims = 1:50, k.param = 30, resolution = 0.5,
                         method = "complete") {
  if (datatype == "sc"){

    # check Seurat object
    info <- try(Seurat::Project(obj), silent = TRUE)
    if (grepl("Error", info) == TRUE) {
      stop(cat("The input should be Seurat Object for single cell data type.\n"))
    }

    # check data process
    if (is.null(obj@reductions$pca) == TRUE){
      stop(cat("The Seurat object hasn't been processed, please use 'scqc' function first.\n"))
    }

    # Add cluster if it has its own cluster
    if (is.null(cluster_list) == FALSE){
      if (length(cluster_list) != ncol(obj)){
        stop(cat("The length of 'cluster_list' doesn't match the number of cells, please check it.\n"))
      }
      obj@meta.data$seurat_clusters <- cluster_list
      # return Seurat object
      return(obj)
    } else{
      # Cluster the cells
      obj <- Seurat::FindNeighbors(obj, dims = dims, k.param = 20)
      obj <- Seurat::FindClusters(obj, resolution = resolution)
      # return Seurat object
      return(obj)
    }
  } else{
    obj <- as.matrix(obj)
    dd <- dist(t(obj))
    cluster_cell <- hclust(dd, method = method)
    return(list(data = obj, cluster_cell = cluster_cell))
  }
}
