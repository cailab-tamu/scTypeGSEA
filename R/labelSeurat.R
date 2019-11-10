#' Add cell type to Seurat object
#'
#' @importFrom Seurat RenameIdents
#' @param obj A Seurat object with cluster.
#' @param cluster_celltype A matrix. Each row is one cell type, contains NSF and padj.
#'
#' @return A list contains two entries. The first is a Seurat object with cell type for each cell and the second is a cell matrix.
#' @export
#'
#' @examples
#' pbmc_example <- doClustering(pbmc_test, nfeatures = 100, npcs = 10,
#'                               dims = 1:10, k.param = 5, resolution = 0.75)
#' cluster_list <- getFC(pbmc_example, min.pct = 0.25, test.use = 'MAST')
#' cluster_celltype <- doGSEA(cluster_list = cluster_list, minSize = 5)
#' pbmc_example <- labelSeurat(pbmc_example, cluster_celltype)$obj
#'
labelSeurat <- function(obj, cluster_celltype){

  # Add cell type for Seurat object
  new.cluster.ids <- cluster_celltype[, 1]
  names(new.cluster.ids) <- levels(obj)
  obj <- Seurat::RenameIdents(obj, new.cluster.ids)

  # generate cell matrix
  cell_mat <- matrix(NA, nrow = ncol(obj), ncol = 3)
  rownames(cell_mat) <- colnames(obj)
  colnames(cell_mat) <- c("ClusterID", "Cell Type", "padj")
  cell_mat[, 1] <- obj@active.ident
  ncluster <- nrow(cluster_celltype)
  for (i in 1:ncluster){
    index <- which(cell_mat[, 1] == i)
    cell_mat[index, 2] <- cluster_celltype[i, 1]
    cell_mat[index, 3] <- cluster_celltype[i, 3]
  }

  return(list(obj = obj, cell_mat = cell_mat))
}
