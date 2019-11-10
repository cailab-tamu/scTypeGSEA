pbmc_example <- scqc(pbmc_small, min.cells = 1, min.features = 10)
pbmc_example <- doClustering(pbmc_example, nfeatures = 100, npcs = 10, dims = 1:10, k.param = 5, resolution = 0.75)
cluster_list <- getFC(pbmc_example, min.pct = 0.25, test.use = "MAST")
cluster_celltype <- doGSEA(cluster_list = cluster_list, minSize = 5)
pbmc_example <- labelSeurat(pbmc_example, cluster_celltype)$obj
check_label <- match(pbmc_example@active.ident, cluster_celltype)

test_that("addcelltype_Seurat works", {
  # check substring
  expect_equal(sum(is.na(check_label)), 0)
  # check label with cluster is 1 to 1
  expect_equal(sum(table(check_label, pbmc_example@meta.data$seurat_clusters)),  ncol(pbmc_example))
})
