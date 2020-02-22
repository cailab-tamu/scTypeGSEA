#' Add cell type to Seurat object.
#'
#' @importFrom Seurat RenameIdents
#' @importFrom Seurat Project
#'
#' @param obj A Seurat object with cluster.
#' @param cluster_celltype A matrix. Each row is one cell type, contains NSF and padj.
#'
#' @return It will return a list including a Seurat object with cell type for each cell and a cell type matrix.
#'
#' @export
#'
#' @examples
#' # It may take several seconds to run the example.
#' # pbmc_example <- scqc(small_RNA, min.cells = 1, min.features = 10, nfeatures = 100, npcs = 10)
#' # pbmc_example <- doClustering(pbmc_example, dims = 1:10, k.param = 5, resolution = 0.75)
#' # cluster_list <- getFC(pbmc_example, min.pct = 0.25, test.use = "MAST")
#' # cluster_celltype <- doGSEA(cluster_list = cluster_list, minSize = 5)
#' # celltype <- labelCelltype(pbmc_example, cluster_celltype)
#'
labelCelltype <- function(obj, cluster_celltype){

    # check Seurat object
    info <- try(Seurat::Project(obj), silent = TRUE)
    if (grepl("Error", info) == TRUE) {
      stop(cat("The input should be Seurat Object for single cell data type.\n"))
    }

    # check number of cluster
    if (nrow(cluster_celltype) != length(unique(obj@meta.data$seurat_clusters))){
      stop(cat("The number of cluster between 'cluster_celltype' and Seurat object is not consistent."))
    }

    # Add cell type for Seurat object
    new.cluster.ids <- as.character(cluster_celltype[, 1])
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
}
