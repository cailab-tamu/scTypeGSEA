#' Checks if the input is already clustered, if not, the function should do it.
#' @param obj A Seurat object
#' @import Seurat
#' @return If the input is already clustered, return the original object. If not, the function will return a Seurat object will cluster.
#' @export
#'
#' @examples
#' library(Seurat)
#' inputMatrix <- matrix(data = rnbinom(n = 1e6, size = 10, prob = .9), nrow = 5000)
#' rownames(inputMatrix) <- paste0('Gene_', seq_len(nrow(inputMatrix)))
#' colnames(inputMatrix) <- paste0('Cell_', seq_len(ncol(inputMatrix)))
#' dta <- CreateSeuratObject(counts = inputMatrix, min.cells = 3, min.features = 20)
#' dta <- check_cluster(dta)

check_cluster <- function(obj){
  check_clu <- is.null(obj@meta.data$seurat_clusters)
  if (check_clu == TRUE){
    # Normalizing the data
    obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000)

    # Identification of highly variable features (feature selection)
    obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000)

    # Scaling the data
    all.genes <- rownames(obj)
    obj <- ScaleData(obj, features = all.genes)

    # Perform linear dimensional reduction
    obj <- RunPCA(obj, features = VariableFeatures(object = obj), verbose = FALSE)

    #Cluster the cells
    obj <- FindNeighbors(obj, dims = 1:50)
    obj <- FindClusters(obj, resolution = 0.5)

    return(obj)
  } else{
    return(obj)
  }
}
