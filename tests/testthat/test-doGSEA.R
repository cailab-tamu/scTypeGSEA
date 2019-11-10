pbmc_example <- scqc(pbmc_small, min.cells = 1, min.features = 10)
pbmc_example <- doClustering(pbmc_example, nfeatures = 100, npcs = 10, dims = 1:10, k.param = 5, resolution = 0.75)
cluster_list <- getFC(pbmc_example, min.pct = 0.25, test.use = "MAST")
cluster_celltype <- doGSEA(cluster_list = cluster_list, minSize = 5)

test_that("GSEA_analysis works", {
  # check each cluster has a cell type
  expect_equal(nrow(cluster_celltype), length(unique(pbmc_example@meta.data$seurat_clusters)))
  # check vector
  expect_equal(is.matrix(cluster_celltype), TRUE)
})
