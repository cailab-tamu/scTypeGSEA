#' Cluster Determination
#'
#' It will do clustering. If the data type is single cell data, the input must be Seurat object and it will use the “Findcluster” function in the Seurat package
#' For any other data type, it will do hierarchical clustering.
#'
#' @importFrom Seurat Project
#' @importFrom Seurat FindNeighbors
#' @importFrom Seurat FindClusters
#' @importFrom Seurat Idents
#' @importFrom stats dist
#' @importFrom stats hclust
#'
#' @param obj A Seurat object.
#' @param datatype Data type for your data, default is 'datatype = "RNA"', which is used for scRNAseq data.
#' @param cluster_cell The cluster result for cells if it is already known.
#' @param dims An integer value. Define dimensions of reduction to use as input. (Do cluster for single cell data.)
#' @param k.param An integer value. Defines k for the k-nearest neighbor algorithm. (Do cluster for single cell data.)
#' @param resolution Value of the resolution parameter, use a value above (below) 1.0 if you want to obtain a larger (smaller) number of communities. (Do cluster for single cell data.)
#' @param hclustmethod The agglomeration method to be used for hierarchical clustering, defalut is "complete".
#' @param ncluster An integer, which is the number of cluster when your input including results from hierarchical clustering.
#'
#' @return It will return a Seurat object with cluster.
#'
#' @export
#'
#' @examples
#' # It may take several seconds to run the example.
#' pbmc_example <- scqc(small_RNA, min.cells = 1, min.features = 10, nfeatures = 100, npcs = 10)
#' pbmc_example <- doClustering(pbmc_example, dims = 1:10, k.param = 5, resolution = 0.75)
#' head(pbmc_example@meta.data$seurat_clusters)
doClustering <- function(obj, datatype = "RNA", cluster_cell = NULL, dims = 1:50, k.param = 30, resolution = 0.5,
                         hclustmethod = "complete", ncluster = 3) {
    # check Seurat object
    info <- try(Seurat::Project(obj), silent = TRUE)
    if (grepl("Error", info) == TRUE) {
      stop(cat("The input should be Seurat Object.\n"))
    }

    # Add cluster if it has its own cluster
    if (is.null(cluster_cell) == FALSE){
      if (length(cluster_cell) != ncol(obj)){
        stop(cat("The length of 'cluster_cell' doesn't match the number of cells, please check it.\n"))
      }
      obj@meta.data$seurat_clusters <- cluster_cell
      # return Seurat object
      return(obj)
    } else{
      # for scRNAseq data
      if (datatype == "RNA"){
        # check data process
        if (is.null(obj@reductions$pca) == TRUE){
          stop(cat("The Seurat object hasn't been processed, please use 'scqc' function first.\n"))
        }
        # Cluster the cells
        obj <- Seurat::FindNeighbors(obj, dims = dims, k.param = 20)
        obj <- Seurat::FindClusters(obj, resolution = resolution)
        # return Seurat object
        return(obj)
      } else{
        # do Hierarchical Clustering
        dta <- as.matrix(obj@assays[[1]]@counts)
        dd <- dist(t(dta))
        cluster_cell <- hclust(dd, method = hclustmethod)
        cluster_membership <- cutree(cluster_cell, k = ncluster)
        Seurat::Idents(obj) <- cluster_membership - 1
        obj@meta.data$seurat_clusters <- cluster_membership - 1
        return(obj)
      }
    }
}
