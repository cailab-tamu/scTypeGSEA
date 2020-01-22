#' Title Assign cell type.
#'
#' This is the main function of scTypeGSEA, which can do quality control, data pre-process, cluster, get fold changes, do GSEA and label the cell in one step.
#'
#' @param obj A seurat object or or any matrix where each column is a cell.
#' @param datatype Data type for your data, which can be "sc" for single cell data and "others" for all other data type.
#' @param percent.mt A decimal value between 0 and 1. Define the highest percentage of reads that map to the mitochondrial genome.
#' @param min.cells An integer value. Include features detected in at least this many cells.
#' @param min.features An integer value. Include cells where at least this many features are detected.
#' @param oversd Remove cells whose library size is greater than mean + oversd * sd. Default is null, which doesn't remove cells.
#' @param normalization.method Method for normalization. Include 'LogNormalize', 'CLR' and 'RC'.
#' @param scale.factor Sets the scale factor for cell-level normalization.
#' @param selection.method How to choose top variable features. Include 'vst', 'mean.var.plot' and 'dispersion'.
#' @param nfeatures An integer value. Define the number of features to select as top variable features.
#' @param npcs An integer value. Define total Number of PCs to compute and store (50 by default).
#' @param cluster_cell The cluster result for cells if it is already known.
#' @param dims An integer value. Define dimensions of reduction to use as input. (Do cluster for single cell data..)
#' @param k.param An integer value. Defines k for the k-nearest neighbor algorithm. (Do cluster for single cell data..)
#' @param resolution Value of the resolution parameter, use a value above (below) 1.0 if you want to obtain a larger (smaller) number of communities. (Do cluster for single cell data..)
#' @param method The agglomeration method to be used for hierarchical clustering, defalut is "complete". (Do cluster for other data type.)
#' @param min.pct A decimal value between 0 and 1. Only test genes that are detected in a minimum fraction of min.pct cells in either of the two populations. Meant to speed up the function by not testing genes that are very infrequently expressed.
#' @param logfc.threshold A decimal value between 0 and 1. Limit testing to genes which show, on average, at least X-fold difference (log-scale) between the two groups of cells. Increasing logfc.threshold speeds up the function, but can miss weaker signals.
#' @param test.use Denotes which test to use. Available options are 'wilcox', 'bimod', 'roc', 'negbinom', 'poisson', 'LR', 'MAST' and 'DESeq2'.
#' @param ncluster An integer, which is the number of cluster when your input including results from hierarchical clustering.
#' @param db The cell type data base to use. For single cell data, we provide two data base, one is 'PanglaoDB' data base (db = 'PanglaoDB_list'), the other one is 'GSEA' data base (db = 'GSEA_list').It can also be a path to the new (referential) data base that hope to be used, the file must be 'rds' format.
#' @param minSize An integer value. Minimal size of a gene set to test. All pathways below the threshold are excluded.
#' @param maxSize An integer value. Maximal size of a gene set to test. All pathways above the threshold are excluded.
#'
#' @return For single cell type of data, it will return a list including a Seurat object with cell type for each cell and a cell type matrix. For other data type, it will only return cell type matrix.
#'
#' @export
#'
#' @examples
#' pbmc_res <- assignCellType(obj = pbmc_small, min.cells = 1, min.features = 10, nfeatures = 100,
#'                            npcs = 10, dims = 1:10, k.param = 5, resolution = 0.75,
#'                            min.pct = 0.25, test.use = "MAST", minSize = 5)
assignCellType <- function(obj, datatype = "sc", min.cells = 3, min.features = 200, percent.mt = 5, oversd = NULL,
                           normalization.method = "LogNormalize", scale.factor = 10000, selection.method = "vst", nfeatures = 2000,
                           npcs = 50, cluster_cell = NULL, dims = 1:50, k.param = 20, resolution = 0.5, method = "complete",
                           min.pct = 0.1, test.use = "wilcox", ncluster = 3,
                           logfc.threshold = 0.25, db = "PanglaoDB_list",minSize = 15, maxSize = 500){
  if (datatype = "sc"){
    # quality control and data pre-process
    obj <- scqc(obj, min.cells = min.cells, min.features = min.features, percent.mt = percent.mt, oversd = oversd, normalization.method = normalization.method,
                scale.factor = scale.factor, selection.method = selection.method, nfeatures = nfeatures, npcs = npcs)

    # do clustering
    obj <- doClustering(obj, datatype = datatype, cluster_cell = cluster_cell, dims = dims, k.param = k.param, resolution = resolution,
                        method = method)

    # get fold changes
    cluster_list <- getFC(obj, datatype = datatype, min.pct = min.pct, test.use = test.use, logfc.threshold = logfc.threshold,
                          ncluster = ncluster)

    # Do GSEA
    cluster_celltype <- doGSEA(cluster_list = cluster_list, db = db, minSize = minSize, maxSize = maxSize)

    # label the seurat object
    return(labelCelltype(obj, datatype = datatype, cluster_celltype))
  } else{
    # do clustering
    obj <- doClustering(obj, datatype = datatype, cluster_list = cluster_list, dims = dims, k.param = k.param, resolution = resolution,
                        method = method)

    # get fold changes
    cluster_list <- getFC(obj, datatype = datatype, min.pct = min.pct, test.use = test.use, logfc.threshold = logfc.threshold,
                          ncluster = ncluster)

    # Do GSEA
    cluster_celltype <- doGSEA(cluster_list = cluster_list, db = db, minSize = minSize, maxSize = maxSize)

    # return cell type matrix
    return(labelCelltype(obj, datatype = datatype, cluster_celltype))

  }
}
