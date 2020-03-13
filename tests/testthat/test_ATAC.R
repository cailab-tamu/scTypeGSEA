# annotation.file <- "~/Documents/Single cell/package example/R package/atac2rna/dataset/Homo_sapiens.GRCh37.82.gtf"
# ATAC_example <- atac2rna(small_pbmc_rna, annotation.file = annotation.file)
# test_that("scqc works", {
#   # check RNA matrix
#   expect_equal(is.null(ATAC_example@assays$RNA), FALSE)
#   # check dimension reduction
#   expect_equal(DefaultAssay(ATAC_example), "RNA")
# })
#
# ATAC_example <- scqc(ATAC_example, min.cells = 1, min.features = 10, nfeatures = 100, npcs = 10)
# test_that("scqc works", {
#   # check it becomes a Seurat object
#   expect_equal(grepl("Error", try(Seurat::Project(ATAC_example), silent = TRUE)), FALSE)
#   # check dimension reduction
#   expect_equal(is.null(ATAC_example@reductions$pca), FALSE)
# })
#
# ATAC_example <- doClustering(ATAC_example, dims = 1:10, k.param = 5, resolution = 0.75)
# test_that("doClustering works", {
#   # check there is cluster
#   expect_equal(is.null(ATAC_example@meta.data$seurat_clusters), FALSE)
# })
#
# cluster_list <- getFC(ATAC_example, min.pct = 0.25, test.use = "MAST")
# test_that("getFC works", {
#   # check every cluster has a list
#   expect_equal(length(cluster_list), length(unique(ATAC_example@meta.data$seurat_clusters)))
#   # check list is not null
#   expect_equal(is.null(cluster_list[[1]]), FALSE)
# })
#
# cluster_celltype <- doGSEA(cluster_list = cluster_list, minSize = 5)
# test_that("doGSEA works", {
#   # check each cluster has a cell type
#   expect_equal(nrow(cluster_celltype), length(unique(ATAC_example@meta.data$seurat_clusters)))
# })
#
# ATAC_example_res1 <- labelCelltype(ATAC_example, cluster_celltype)
# set.seed(1234)
# ATAC_example_res2 <- assignCellType(pbmc_small_atac, datatype = "ATAC", annotation.file = annotation.file,
#                                     min.cells = 1, min.features = 10, nfeatures = 100, npcs = 10,
#                                     dims = 1:10, k.param = 5, resolution = 0.75,
#                                     min.pct = 0.25, test.use = "MAST", minSize = 5)
# test_that("assignCellType works", {
#   expect_equal(ATAC_example_res1$cell_mat, ATAC_example_res2$cell_mat)
# })
