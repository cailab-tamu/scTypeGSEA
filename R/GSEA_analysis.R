#' Title Runs preranked gene set enrichment analysis.
#'
#' The function do gene set enrichment analysis(GSEA) for each cluster with its gene ranks. It will return cell type for each cluster.
#'
#' @importFrom fgsea fgsea
#' @param obj A Seurat object with cluster.
#' @param cluster_list A ranked gene list for each cluster.
#' @param db The cell type data base to use. It should be 'PanglaoDB_list' for 'PanglaoDB' data base and 'GSEA_list' for 'GSEA' data base.
#' @param minSize Minimal size of a gene set to test. All pathways below the threshold are excluded.
#' @param maxSize Maximal size of a gene set to test. All pathways above the threshold are excluded.
#' @param nperm Number of permutations to do. Minimial possible nominal p-value is about 1/nperm
#'
#' @return Cell type for each cluster.
#' @export
#'
GSEA_analysis <- function(obj, cluster_list, db = "PanglaoDB_list", minSize = 15, maxSize = 500, nperm = 10000) {

  # number of cluster
  ncluster <- length(cluster_list)
  cluster_celltype <- rep(NA, ncluster)
  names(cluster_celltype) <- names(cluster_list)

  # decide the database to use
  if (db == db) {
    pathways <- readRDS("inst/PanglaoDB_list.rds")
  } else {
    pathways <- readRDS("inst/GSEA_list.rds")
  }

  # Do fgsea to each cluster
  for (i in 1:ncluster) {
    cat(paste0("Do GSEA for cluster"), i - 1, "\n")

    Ranks <- cluster_list[[i]]

    fgseaRes <- fgsea::fgsea(pathways = pathways, stats = Ranks, minSize = minSize, maxSize = maxSize, nperm = nperm)

    fgseaRes <- fgseaRes[order(fgseaRes[, "NES"], -fgseaRes[, "padj"], decreasing = TRUE), ]

    # decide the cell type
    if (fgseaRes[, "padj"][1] > 0.05) {
      cluster_celltype[i] <- "unidentified"
    } else {
      cluster_celltype[i] <- fgseaRes[, "pathway"][1]
    }
  }

  return(cluster_celltype)
}
