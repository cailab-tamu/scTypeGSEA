#' Title Cluster Determination
#'
#' If the input is already clustered, just return it.
#' If not, the function will check whether the data has already been pre-processed, for example, PCA. If not, the function will do data pre-process first.
#' After the basic data pre-process (Normalization, Scale data, Find HVG and PCA), the function will do cluster.
#'
#' @importFrom Seurat PercentageFeatureSet
#' @importFrom Seurat NormalizeData
#' @importFrom Seurat FindVariableFeatures
#' @importFrom Seurat ScaleData
#' @importFrom Seurat RunPCA
#' @importFrom Seurat VariableFeatures
#' @importFrom Seurat FindNeighbors
#' @importFrom Seurat FindClusters
#'
#' @param obj An Seurat object.
#' @param percent.mt The highest percentage of reads that map to the mitochondrial genome.
#' @param normalization.method Method for normalization. Include 'LogNormalize', 'CLR' and 'RC'.
#' @param scale.factor Sets the scale factor for cell-level normalization.
#' @param selection.method How to choose top variable features. Include 'vst', 'mean.var.plot' and 'dispersion'.
#' @param nfeatures Number of features to select as top variable features.
#' @param npcs Total Number of PCs to compute and store (50 by default).
#' @param dims Dimensions of reduction to use as input. (For cluster.)
#' @param k.param Defines k for the k-nearest neighbor algorithm. (For cluster.)
#' @param resolution Value of the resolution parameter, use a value above (below) 1.0 if you want to obtain a larger (smaller) number of communities. (For cluster.)
#'
#' @return If the input is already clustered, just return it. If not, the function do cluster and return obj with cluster.
#' @export
#'
check_cluster <- function(obj, percent.mt = 5, normalization.method = "LogNormalize", scale.factor = 10000, selection.method = "vst", nfeatures = 2000, npcs = 50, dims = 1:50, k.param = 20, resolution = 0.5) {

  # check cluster
  check_clu <- is.null(obj@meta.data$seurat_clusters)

  if (check_clu == TRUE) {

    cat("The data has not been pre-processed, let's do it!\n")
    # check PCA dimension reduction
    check_pca <- is.null(obj@reductions$pca)
    if (check_pca == TRUE) {
      # Do quality control
      obj[["percent.mt"]] <- Seurat::PercentageFeatureSet(obj, pattern = "^MT-")
      obj <- subset(obj, subset = percent.mt < percent.mt)

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
  } else {
    return(obj)
  }
}
