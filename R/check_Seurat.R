
#' Check whether the object is Seurat object or not
#' @importFrom Seurat Project
#' @param obj The object to check whether it is Seurat object or not
#' @return An error message if the objective is not a Seurat objective and "TRUE" if it is a Seurat objective
#' @export
#' @examples
#' library(Seurat)
#' inputMatrix <- matrix(data = rnbinom(n = 1e6, size = 10, prob = .9), nrow = 5000)
#' rownames(inputMatrix) <- paste0('Gene_', seq_len(nrow(inputMatrix)))
#' colnames(inputMatrix) <- paste0('Cell_', seq_len(ncol(inputMatrix)))
#' dta <- CreateSeuratObject(counts = inputMatrix, min.cells = 3, min.features = 20)
#' dta <- check_cluster(dta)

check_Seurat <- function(obj){

  info <- try(Seurat::Project(obj), silent = TRUE)

  if (grepl("Error", info) == TRUE) {
    stop(paste0("The input should be a Seurat object!"))
  } else
    return(TRUE)
}
