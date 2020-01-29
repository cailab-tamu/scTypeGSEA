#' Assign cell type.
#'
#' This is the main function of scTypeGSEA, which can do quality control, data pre-process, cluster, get fold changes, do GSEA and label the cell in one step.
#'
#' @param obj A seurat object or or any matrix where each column is a cell.
#' @param datatype Data type for your data, which can be "RNA" for scRNAseq data, "ATAC" for scATACseq data or any other data type.
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
#' @param hclustmethod The agglomeration method to be used for hierarchical clustering, defalut is "complete". (Do cluster for other data type.)
#' @param ncluster An integer, which is the number of cluster when your input including results from hierarchical clustering.
#' @param min.pct A decimal value between 0 and 1. Only test genes that are detected in a minimum fraction of min.pct cells in either of the two populations. Meant to speed up the function by not testing genes that are very infrequently expressed.
#' @param logfc.threshold A decimal value between 0 and 1. Limit testing to genes which show, on average, at least X-fold difference (log-scale) between the two groups of cells. Increasing logfc.threshold speeds up the function, but can miss weaker signals.
#' @param test.use Denotes which test to use. Available options are 'wilcox', 'bimod', 'roc', 'negbinom', 'poisson', 'LR', 'MAST' and 'DESeq2'. The defalut is "MAST" for scRNAseq data, we suggest to use 'wilcox' for other data type.
#' @param db The cell type data base to use. For single cell data, we provide two data base, one is 'PanglaoDB' data base (db = 'PanglaoDB_list'), the other one is 'GSEA' data base (db = 'GSEA_list'). It can also be a path to the new (referential) data base that hope to be used, the file must be 'rds' format.
#' @param minSize An integer value. Minimal size of a gene set to test. All pathways below the threshold are excluded.
#' @param maxSize An integer value. Maximal size of a gene set to test. All pathways above the threshold are excluded.
#' @param metadata Additional cell-level metadata to add to the Seurat object. Should be a data frame where the rows are cell names and the columns are additional metadata fields. (Only for "ATAC" data)
#' @param fragmentpath Path to tabix-indexed fragments file. (Only for "ATAC" data)
#' @param annotation.file Path to GTF annotation file. (Only for "ATAC" data)
#' @param qualitycontrol Whether to do quality control for peak matrix. If it is TRUE, it will compute QC Metrics and delete some cells. One needs both metadata data and fragments file to achieve quality control. (For "ATAC" data)
#' @param alpha Integer value 0 or 1 to decide which method to use. alpha = 0 is to use annotation file, and alpha = 1 is to use genomic fragments file. (Only for "ATAC" data)
#' @param seq.levels Which seqlevels to keep (corresponds to chromosomes usually). (Only for "ATAC" data)
#' @param include.body Include the gene body? (Only for "ATAC" data)
#' @param upstream Number of bases upstream to consider. (Only for "ATAC" data)
#' @param downstream Number of bases downstream to consider. (Only for "ATAC" data)
#' @param EnsDbobj For toSAF a GRangesList object. (Only for "ATAC" data)
#' @param chunk Number of chunks to use when processing the fragments file. Fewer chunks may enable faster processing, but will use more memory. (Only for "ATAC" data)
#' @param filter A filter describing which results to retrieve from the database. (Only for "ATAC" data)
#'
#' @return It will return the Seurat object with cell type, a cell type matrix and a cluster list.
#'
#' @export
#'
#' @examples
#' pbmc_example <- assignCellType(small_RNA, min.cells = 1, min.features = 10,
#'                                nfeatures = 100, npcs = 10,
#'                                dims = 1:10, k.param = 5, resolution = 0.75,
#'                                min.pct = 0.25, test.use = "MAST", minSize = 5)
#'
assignCellType <- function(obj, datatype = "RNA", min.cells = 3, min.features = 200, percent.mt = 5, oversd = NULL, normalization.method = "LogNormalize",
                           scale.factor = 10000, selection.method = "vst", nfeatures = 2000, npcs = 50,
                           cluster_cell = NULL, dims = 1:50, k.param = 20, resolution = 0.5, hclustmethod = "complete", ncluster = 3,
                           min.pct = 0.1, test.use = "wilcox", logfc.threshold = 0.1,
                           db = "PanglaoDB_list",minSize = 15, maxSize = 500,
                           metadata = NULL, fragmentpath = NULL, annotation.file = NULL, qualitycontrol = FALSE,
                           alpha = 0, seq.levels = c(1:22, "X", "Y"), include.body = TRUE, upstream = 2000, downstream = 0,
                           EnsDbobj = EnsDb.Hsapiens.v75, chunk = 50, filter = ~ gene_biotype == "protein_coding"){
  if (datatype == "ATAC"){
    # transform atac to gene activity matrix
    obj <- atac2rna(obj, metadata = metadata, fragmentpath = metadata, annotation.file = annotation.file, qualitycontrol = qualitycontrol, alpha = alpha,
                                seq.levels = seq.levels, include.body = include.body, upstream = upstream, downstream = downstream,
                                EnsDbobj = EnsDbobj, chunk = chunk, filter = filter)
    datatype = "RNA"
  }
  # quality control and data pre-process
  obj <- scqc(obj, datatype = datatype, min.cells = min.cells, min.features = min.features, percent.mt = percent.mt, oversd = oversd, normalization.method = normalization.method,
              scale.factor = scale.factor, selection.method = selection.method, nfeatures = nfeatures, npcs = npcs)

  # do clustering
  obj <- doClustering(obj, datatype = datatype, cluster_cell = cluster_cell, dims = dims, k.param = k.param, resolution = resolution,
                      hclustmethod = hclustmethod, ncluster = ncluster)

  # get fold changes
  cluster_list <- getFC(obj, min.pct = min.pct, test.use = test.use, logfc.threshold = logfc.threshold)

  # Do GSEA
  cluster_celltype <- doGSEA(cluster_list = cluster_list, db = db, minSize = minSize, maxSize = maxSize)

  # label the seurat object
  obj <- labelCelltype(obj, cluster_celltype)
  return(list(Seurat_obj = obj$obj, cell_mat = obj$cell_mat, cluster_list = cluster_list, cluster_celltype = cluster_celltype))
}
