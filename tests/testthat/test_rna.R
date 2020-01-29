set.seed(1234)
load("data/small_RNA.rda")
load("data/PanglaoDB_list.rda")

pbmc_example <- scqc(small_RNA, min.cells = 1, min.features = 10, nfeatures = 100, npcs = 10)
test_that("scqc works", {
  # check it becomes a Seurat object
  expect_equal(grepl("Error", try(Seurat::Project(pbmc_example), silent = TRUE)), FALSE)
  # check dimension reduction
  expect_equal(is.null(pbmc_example@reductions$pca), FALSE)
})

pbmc_example <- doClustering(pbmc_example, dims = 1:10, k.param = 5, resolution = 0.75)
test_that("doClustering works", {
  # check there is cluster
  expect_equal(is.null(pbmc_example@meta.data$seurat_clusters), FALSE)
})

cluster_list <- getFC(pbmc_example, min.pct = 0.25, test.use = "MAST")
test_that("getFC works", {
  # check every cluster has a list
  expect_equal(length(cluster_list), length(unique(pbmc_example@meta.data$seurat_clusters)))
  # check list is not null
  expect_equal(is.null(cluster_list[[1]]), FALSE)
})

cluster_celltype <- doGSEA(cluster_list = cluster_list, minSize = 5)
test_that("doGSEA works", {
  # check each cluster has a cell type
  expect_equal(nrow(cluster_celltype), length(unique(pbmc_example@meta.data$seurat_clusters)))
})

pbmc_example_res1 <- labelCelltype(pbmc_example, cluster_celltype)
set.seed(1234)
pbmc_example_res2 <- assignCellType(pbmc_small, min.cells = 1, min.features = 10, nfeatures = 100, npcs = 10,
                                    dims = 1:10, k.param = 5, resolution = 0.75,
                                    min.pct = 0.25, test.use = "MAST", minSize = 5)
test_that("assignCellType works", {
  expect_equal(pbmc_example_res1$cell_mat, pbmc_example_res2$cell_mat)
})






