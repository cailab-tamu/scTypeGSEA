library(Seurat)
inputData <- Read10X('~/Documents/Single cell/package example/R package/scTypeGSEA/pbmc1k/filtered_feature_bc_matrix/')
obj <- inputData$`Antibody Capture`
obj <- scqc(obj, datatype = "Antibody Capture", min.cells = 1, min.features = 10)
obj <- doClustering(obj, datatype = "Antibody Capture")
cluster_list1 <- getFC(obj, min.pct = 0.25, test.use = "wilcox")
cluster_celltype1 <- data.frame("cell type" = paste0("type", 1:3), "NES" = c(2.324, 2.63, 2.765), "padj" = c(3.54e-9, 4.36e-5, 2.546e-8))
obj_res <- labelCelltype(obj, cluster_celltype1)
