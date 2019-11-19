#' Cluster Determination
#'
#' If the input is already clustered, just return it.
#' If not, the function will check whether the data has already been pre-processed, for example, PCA. If not, the function will do data pre-process first.
#' After the basic data pre-process (Normalization, Scale data, Find HVG and PCA), the function will do cluster.
#'
#' @importFrom Seurat NormalizeData
#' @importFrom Seurat FindVariableFeatures
#' @importFrom Seurat ScaleData
#' @importFrom Seurat RunPCA
#' @importFrom Seurat VariableFeatures
#' @importFrom Seurat FindNeighbors
#' @importFrom Seurat FindClusters
#'
#' @param obj An Seurat object.
#' @param normalization.method Method for normalization. Include 'LogNormalize', 'CLR' and 'RC'.
#' @param scale.factor Sets the scale factor for cell-level normalization.
#' @param selection.method How to choose top variable features. Include 'vst', 'mean.var.plot' and 'dispersion'.
#' @param nfeatures An integer value. Define the number of features to select as top variable features.
#' @param npcs An integer value. Define total Number of PCs to compute and store (50 by default).
#' @param dims An integer value. Define dimensions of reduction to use as input. (For cluster.)
#' @param k.param An integer value. Defines k for the k-nearest neighbor algorithm. (For cluster.)
#' @param resolution Value of the resolution parameter, use a value above (below) 1.0 if you want to obtain a larger (smaller) number of communities. (For cluster.)
#' @param doit A boolean value (TRUE/FALSE). If true, the function will do the cluster with new parameters even the object contains original cluster. (Default is FALSE.)
#' @param doprocess A boolean value (TRUE/FALSE). If true, the function will do data process with new parameters. (Default is FALSE.)
#'
#' @return If the input is already clustered, just return it. If not, the function do cluster and return obj with cluster.
#' @export
#' @examples
#' pbmc_example <- scqc(pbmc_small, min.cells = 1, min.features = 10)
#' is.null(pbmc_example@meta.data$seurat_clusters)
#' pbmc_example <- doClustering(pbmc_example, nfeatures = 100, npcs = 10,
#'                               dims = 1:10, k.param = 5, resolution = 0.5)
#' head(pbmc_example@meta.data$seurat_clusters)
#' pbmc_example <- doClustering(pbmc_example, nfeatures = 100, npcs = 10,
#'                               dims = 1:10, k.param = 5, resolution = 0.75, doit = TRUE)
#' head(pbmc_example@meta.data$seurat_clusters)
doClustering <- function(obj, normalization.method = "LogNormalize", scale.factor = 10000, selection.method = "vst", nfeatures = 2000, npcs = 50, dims = 1:50, k.param = 30, resolution = 0.5, doit = "FALSE", doprocess = "FALSE") {

  # check cluster
  check_clu <- is.null(obj@meta.data$seurat_clusters)

  if (check_clu == FALSE & doit == FALSE){
    return(obj)
  } else{
    # check PCA dimension reduction
    check_pca <- is.null(obj@reductions$pca)
    if (check_pca == TRUE || doprocess == TRUE) {

      if (check_pca == TRUE) cat("The data has not been pre-processed, let's do it!\n")

      # Normalize object
      obj <- Seurat::NormalizeData(obj, normalization.method = normalization.method, scale.factor = scale.factor)

      # Identification of highly variable features (feature selection)
      obj <- Seurat::FindVariableFeatures(obj, selection.method = selection.method, nfeatures = nfeatures)

      # Scaling the data
      all.genes <- rownames(obj)
      obj <- Seurat::ScaleData(obj, features = all.genes)

      # Perform linear dimensional reduction
      obj <- Seurat::RunPCA(obj, features = Seurat::VariableFeatures(object = obj), verbose = FALSE, npcs = npcs)
    }

    # Cluster the cells
    obj <- Seurat::FindNeighbors(obj, dims = dims, k.param = 20)
    obj <- Seurat::FindClusters(obj, resolution = resolution)

    return(obj)
  }
}
