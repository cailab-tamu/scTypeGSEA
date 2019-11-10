pbmc_example <- doClustering(pbmc_test, nfeatures = 100, npcs = 10, dims = 1:10, k.param = 5, resolution = 0.5)
test_that("check_cluster works", {
  # check there is cluster
  expect_equal(is.null(pbmc_example@meta.data$seurat_clusters), FALSE)
})
