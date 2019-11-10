pbmc_example <- scqc(pbmc_small, min.cells = 1, min.features = 10)
pbmc_example <- doClustering(pbmc_example, nfeatures = 100, npcs = 10, dims = 1:10, k.param = 5, resolution = 0.75)
test_that("check_cluster works", {
  # check there is cluster
  expect_equal(is.null(pbmc_example@meta.data$seurat_clusters), FALSE)
})
