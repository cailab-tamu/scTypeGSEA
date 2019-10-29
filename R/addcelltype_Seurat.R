#' Title Add cell type to Seurat object
#'
#' @importFrom Seurat RenameIdents
#' @param obj A Seurat object with cluster.
#' @param cluster_celltype A list, the name is cluster and the entry is cell type.
#'
#' @return A Seurat object with cell type for each cell.
#' @export
#'
#' @examples
#' pbmc_example <- check_cluster(pbmc_small, nfeatures = 100, npcs = 10,
#'                               dims = 1:10, k.param = 5, resolution = 0.75)
#' cluster_celltype <- list(Cluster0 = "NK cells", Cluster1 = "Dendritic cells")
#' pbmc_example <- addcelltype_Seurat(pbmc_example, cluster_celltype)
#' head(pbmc_example@active.ident)
#'
addcelltype_Seurat <- function(obj, cluster_celltype){

  # Add cell type for Seurat object
  new.cluster.ids <- cluster_celltype
  names(new.cluster.ids) <- levels(obj)
  obj <- Seurat::RenameIdents(obj, new.cluster.ids)
  return(obj)
}
