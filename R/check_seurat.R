#' Check Seurat.
#'
#' Check whether the input is a seurat object or not.
#'
#' @importFrom Seurat CreateSeuratObject
#' @param obj A seurat object or a count gene expression.
#' @param min.cells Include features detected in at least this many cells. Will subset the counts matrix as well. To reintroduce excluded features, create a new object with a lower cutoff.
#' @param min.features Include cells where at least this many features are detected.
#'
#' @return If the input is a seurat object, just return it. If not, it will assume it is a count gene expression matrix and transform it as a seurat object.
#' @export
#'
#' @examples
#' pbmc_example <- check_seurat(pbmc_raw, min.cells = 1, min.features = 10)
#' pbmc_example
check_seurat <- function(obj, min.cells = 3, min.features = 200) {
  # check Seurat object
  info <- try(Seurat::Project(obj), silent = TRUE)

  if (grepl("Error", info) == TRUE) {
    cat("The input is not a Seurat object, next to transform gene express count matrix to Seurat object.")
    obj <- Seurat::CreateSeuratObject(counts = obj, min.cells = min.cells, min.features = min.features)
    return(obj)
  } else{
    return(obj)
  }
}
