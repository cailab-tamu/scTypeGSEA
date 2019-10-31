pbmc_example <- check_cluster(pbmc_test, nfeatures = 100, npcs = 10, dims = 1:10, k.param = 5, resolution = 0.75)
cluster_list <- Test_DE_cluster(pbmc_example, min.pct = 0.25, test.use = 'MAST')
cluster_celltype <- GSEA_analysis(cluster_list = cluster_list, minSize = 5, nperm = 1000)

test_that("GSEA_analysis works", {
  # check each cluster has a cell type
  expect_equal(length(cluster_celltype), length(unique(pbmc_example@meta.data$seurat_clusters)))
  # check vector
  expect_equal(is.vector(cluster_celltype), TRUE)
})
