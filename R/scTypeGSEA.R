#' Gene set enrichment analysis for single cell data.
#'
#' @importFrom Seurat Project
#' @importFrom Seurat RenameIdents
#' @param obj An Seurat object.
#' @param ...
#'
#' @return An Seurat object with celltype for each cluster.
#' @export
#'
scTypeGSEA <- function(obj, ...) {

  # check Seurat object
  info <- try(Seurat::Project(obj), silent = TRUE)

  if (grepl("Error", info) == TRUE) {
    stop(paste0("The input should be a Seurat object!"))
  }

  # check Cluster
  obj <- check_cluster(obj)

  # Run DE test to get genes rank for each cluster
  cluster_list <- Test_DE_cluster(obj)

  # Do GSEA analysis
  cluster_celltype <- GSEA_analysis(obj, cluster_list = cluster_list)

  # Add cell type for Seurat object
  new.cluster.ids <- cluster_celltype
  names(new.cluster.ids) <- levels(obj)
  obj <- Seurat::RenameIdents(obj, new.cluster.ids)
  return(obj)
}
