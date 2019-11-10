#' Single cell quality control.
#'
#' Check whether the input is a seurat object or not and the do quality control.
#'
#' @importFrom Seurat CreateSeuratObject
#' @importFrom Seurat PercentageFeatureSet
#' @importFrom stats sd
#'
#' @param obj A seurat object or a count gene expression.
#' @param percent.mt The highest percentage of reads that map to the mitochondrial genome.
#' @param min.cells Include features detected in at least this many cells. Will subset the counts matrix as well. To reintroduce excluded features, create a new object with a lower cutoff.
#' @param min.features Include cells where at least this many features are detected.
#' @param oversd Remove cells whose library size is greater than mean + oversd * sd. Default is null, which doesn't remove cells.
#'
#' @return If the input is not Seurat object, it will assume it is a count gene expression matrix and transform it as a seurat object. And then do quality control to select cells. Otherwise, it will do quality control directly.
#' @export
#'
#' @examples
#' pbmc_example <- scqc(pbmc_small, min.cells = 1, min.features = 10)
#' pbmc_example
scqc <- function(obj, min.cells = 3, min.features = 200, percent.mt = 5, oversd = NULL) {
  # check Seurat object
  info <- try(Seurat::Project(obj), silent = TRUE)

  if (grepl("Error", info) == TRUE) {
    cat("The input is not a Seurat object, next to transform gene express count matrix to Seurat object.\n")
    obj <- Seurat::CreateSeuratObject(counts = obj, min.cells = min.cells, min.features = min.features)
  }

  # Do quality control
  obj[["percent.mt"]] <- Seurat::PercentageFeatureSet(obj, pattern = "^MT-")
  obj <- subset(obj, cells = which(obj[["percent.mt"]] < percent.mt))

  if (is.null(oversd) == TRUE){
    return(obj)
  } else{
    library_size <- obj$nCount_RNA
    sd_library_size <- sd(library_size)
    mean_library_size <- mean(library_size)
    thres_library_size <- mean_library_size + oversd * sd_library_size
    obj <- subset(obj, cells = which(obj[["nCount_RNA"]] < thres_library_size))
    return(obj)
  }
}
