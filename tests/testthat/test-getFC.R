pbmc_example <- doClustering(pbmc_test, nfeatures = 100, npcs = 10, dims = 1:10, k.param = 5, resolution = 0.75)
cluster_list <- getFC(pbmc_example, min.pct = 0.25, test.use = "MAST")

test_that("Test_DE_cluster works", {
  # check every cluster has a list
  expect_equal(length(cluster_list), length(unique(pbmc_example@meta.data$seurat_clusters)))
  # check list is not null
  expect_equal(is.null(cluster_list[[1]]), FALSE)
  # check list
  expect_equal(is.list(cluster_list), TRUE)
})
