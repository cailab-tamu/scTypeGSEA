#' Add cell type to Seurat object or other data type.
#'
#' @importFrom Seurat RenameIdents
#' @importFrom Seurat Project
#'
#' @param obj A Seurat object with cluster or any matrix where each column is a cell.
#' @param datatype Data type to label cells, which can be "sc" for single cell data and "others" for all other data type.
#' @param cluster_celltype A matrix. Each row is one cell type, contains NSF and padj.
#'
#' @return For single cell type of data, it will return a list including a Seurat object with cell type for each cell and a cell type matrix. For other data type, it will only return cell type matrix.
#'
#' @export
#'
#' @examples
#' pbmc_example <- scqc(pbmc_small, min.cells = 1, min.features = 10, nfeatures = 100, npcs = 10)
#' pbmc_example <- doClustering(pbmc_example, dims = 1:10, k.param = 5, resolution = 0.75)
#' cluster_list <- getFC(pbmc_example, min.pct = 0.25, test.use = "MAST")
#' cluster_celltype <- doGSEA(cluster_list = cluster_list, minSize = 5)
#' celltype <- labelCelltype(pbmc_example, datatype = "sc", cluster_celltype)
#'
labelCelltype <- function(obj, datatype = "sc", cluster_celltype){

  if (datatype == "sc"){
    # check Seurat object
    info <- try(Seurat::Project(obj), silent = TRUE)
    if (grepl("Error", info) == TRUE) {
      stop(cat("The input should be Seurat Object for single cell data type.\n"))
    }

    # Add cell type for Seurat object
    new.cluster.ids <- cluster_celltype[, 1]
    names(new.cluster.ids) <- levels(obj)
    obj <- Seurat::RenameIdents(obj, new.cluster.ids)

    # generate cell type matrix
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
    # return cell type matrix
    return(list(obj = obj, cell_mat = cell_mat))
  } else{
    # generate cell type matrix
    cell_mat <- matrix(NA, nrow = ncol(obj), ncol = 3)
    rownames(cell_mat) <- colnames(obj)
    colnames(cell_mat) <- c("ClusterID", "Cell Type", "padj")
    ncluster <- nrow(cluster_celltype)
    for (i in 1:ncluster){
      index <- which(cell_mat[, 1] == i)
      cell_mat[index, 1] <- paste0("Cluster ", i)
      cell_mat[index, 2] <- cluster_celltype[i, 1]
      cell_mat[index, 3] <- cluster_celltype[i, 3]
    }
    # return cell type matrix
    return(cell_mat)
  }
}
