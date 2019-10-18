
#' Check whether the object is Seurat object or not
#' @importFrom Seurat Project
#' @param obj The object to check whether it is Seurat object or not
#' @return An error message if the objective is not a Seurat objective and "TRUE" if it is a Seurat objective
#' @export

check_Seurat <- function(obj){

  info <- try(Seurat::Project(obj), silent = TRUE)

  if (grepl("Error", info) == TRUE) {
    stop(paste0("The input should be a Seurat object!"))
  } else
    return(TRUE)
}
