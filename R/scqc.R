#' Single cell quality control and basic data pre-process.
#'
#' For single cell type data, it will create a Seurat object and perform single-cell data quality control,
#' including checking for minimum cell library size, mitochondrial ratio, outlier cells, and the fraction of cells where a gene is expressed.
#' And then the function will do basic data pre-process, which includes “Normalization”, “Scale data”, “Find HVG” and “PCA”.
#' For other type of data, it will create a Seurat object to save the matrix.
#'
#' @importFrom Seurat CreateSeuratObject
#' @importFrom Seurat PercentageFeatureSet
#' @importFrom Seurat NormalizeData
#' @importFrom Seurat FindVariableFeatures
#' @importFrom Seurat ScaleData
#' @importFrom Seurat RunPCA
#' @importFrom Seurat VariableFeatures
#' @importFrom stats sd
#'
#' @param obj A seurat object or a raw count gene expression.
#' @param datatype Data type for your data, default is 'datatype = "RNA"', which is used for scRNAseq data.
#' @param percent.mt A decimal value between 0 and 1. Define the highest percentage of reads that map to the mitochondrial genome.
#' @param min.cells An integer value. Include features detected in at least this many cells.
#' @param min.features An integer value. Include cells where at least this many features are detected.
#' @param oversd Remove cells whose library size is greater than mean + oversd * sd. Default is null, which doesn't remove cells.
#' @param normalization.method Method for normalization. Include 'LogNormalize', 'CLR' and 'RC'.
#' @param scale.factor Sets the scale factor for cell-level normalization.
#' @param selection.method How to choose top variable features. Include 'vst', 'mean.var.plot' and 'dispersion'.
#' @param nfeatures An integer value. Define the number of features to select as top variable features.
#' @param npcs An integer value. Define total Number of PCs to compute and store (50 by default).
#'
#' @return A Seurat Object with quality control and basic data pre-process.
#'
#' @export
#'
#' @examples
#' pbmc_example <- scqc(small_RNA, min.cells = 1, min.features = 10, nfeatures = 100, npcs = 10)
#' pbmc_example
scqc <- function(obj, datatype = "RNA", min.cells = 10, min.features = 1000, percent.mt = 5, oversd = NULL, normalization.method = "LogNormalize",
                 scale.factor = 10000, selection.method = "vst", nfeatures = 2000, npcs = 50) {

  if (datatype == "RNA"){
  # check Seurat object
  info <- try(Seurat::Project(obj), silent = TRUE)

  if (grepl("Error", info) == TRUE) {
    cat("The input is not a Seurat object, next to transform gene express count matrix to Seurat object.\n")
    obj <- Seurat::CreateSeuratObject(counts = obj, min.cells = min.cells, min.features = min.features)
  }

  # Do quality control
  obj[["percent.mt"]] <- Seurat::PercentageFeatureSet(obj, pattern = "^MT-")
  obj <- subset(obj, cells = which(obj[["percent.mt"]] < percent.mt))

  if (is.null(oversd) == FALSE){
    library_size <- obj$nCount_RNA
    sd_library_size <- sd(library_size)
    mean_library_size <- mean(library_size)
    thres_library_size <- mean_library_size + oversd * sd_library_size
    obj <- subset(obj, cells = which(obj[["nCount_RNA"]] < thres_library_size))
  }

  # Normalize object
  obj <- Seurat::NormalizeData(obj, normalization.method = normalization.method, scale.factor = scale.factor)

  # Identification of highly variable features (feature selection)
  obj <- Seurat::FindVariableFeatures(obj, selection.method = selection.method, nfeatures = nfeatures)

  # Scaling the data
  all.genes <- rownames(obj)
  obj <- Seurat::ScaleData(obj, features = all.genes)

  # Perform linear dimensional reduction
  obj <- Seurat::RunPCA(obj, features = Seurat::VariableFeatures(object = obj), verbose = FALSE, npcs = npcs)

  # return Seurat object
  return(obj)
  } else {
    obj <- Seurat::CreateSeuratObject(counts = obj, min.cells = min.cells, min.features = min.features, assay = datatype)
    # return Seurat object
    return(obj)
  }
}
