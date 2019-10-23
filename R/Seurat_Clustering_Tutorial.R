# library(Matrix)
# library(ggplot2)
# library(dplyr)
# library(Seurat)
#
# ##Setup the Seurat Object
# pbmc.data <- Read10X(data.dir = "~/Documents/single cell/package example/R package/Seurat/cluster/data set")#32738 x 2700
# pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
# pbmc # 13714 x 2700
#
# ##Standard pre-processing workflow
#
# # The [[ operator can add columns to object metadata. This is a great place to stash QC stats
# pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
# # Visualize QC metrics as a violin plot
# VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
# # FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# # for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
# plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
# plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
# CombinePlots(plots = list(plot1, plot2))
# pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)# 13714 x 2638
#
# ##Normalizing the data
# pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
# #pbmc <- NormalizeData(pbmc)
#
# ##Identification of highly variable features (feature selection)
# pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
# # Identify the 10 most highly variable genes
# top10 <- head(VariableFeatures(pbmc), 10)
# # plot variable features with and without labels
# plot1 <- VariableFeaturePlot(pbmc)
# plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
# CombinePlots(plots = list(plot1, plot2))
#
# ##Scaling the data
# all.genes <- rownames(pbmc)
# pbmc <- ScaleData(pbmc, features = all.genes)
# #pbmc[["RNA"]]@scale.data
# #pbmc <- ScaleData(pbmc)
#
# ##Perform linear dimensional reduction
# pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc), verbose = FALSE)
# # Examine and visualize PCA results a few different ways
# print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)
# VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")
# DimPlot(pbmc, reduction = "pca")
#
# ##Determine the ‘dimensionality’ of the dataset
# # NOTE: This process can take a long time for big datasets, comment out for expediency. More
# # approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# # computation time
# pbmc <- JackStraw(pbmc, num.replicate = 100)
# pbmc <- ScoreJackStraw(pbmc, dims = 1:20)
# JackStrawPlot(pbmc, dims = 1:15)
# ElbowPlot(pbmc)
#
# ##Cluster the cells
# pbmc <- FindNeighbors(pbmc, dims = 1:10)
# pbmc <- FindClusters(pbmc, resolution = 0.5)
# # Look at cluster IDs of the first 5 cells
# head(Idents(pbmc), 5)
#
# ##Run non-linear dimensional reduction (UMAP/tSNE)
# pbmc <- RunTSNE(pbmc, dims = 1:10)
# DimPlot(pbmc, reduction = "tsne")
# pbmc <- RunUMAP(pbmc, dims = 1:10)
# DimPlot(pbmc, reduction = "umap")
#
# ##Finding differentially expressed features (cluster biomarkers)
# cluster1.markers <- FindMarkers(pbmc, ident.1 = 1, min.pct = 0.25)
# head(cluster1.markers, n = 5)
# # find all markers distinguishing cluster 5 from clusters 0 and 3
# cluster5.markers <- FindMarkers(pbmc, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25)
# head(cluster5.markers, n = 5)
# # find markers for every cluster compared to all remaining cells, report only the positive ones
# pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
# pbmc.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
#
# cluster1.markers <- FindMarkers(pbmc, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
#
# VlnPlot(pbmc, features = c("MS4A1", "CD79A"))
#
# # you can plot raw counts as well
# VlnPlot(pbmc, features = c("NKG7", "PF4"), slot = "counts", log = TRUE)
# FeaturePlot(pbmc, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP",
#                                "CD8A"))
# top10 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
# DoHeatmap(pbmc, features = top10$gene) + NoLegend()
#
# ##Assigning cell type identity to clusters
# new.cluster.ids <- c("Naive CD4 T", "Memory CD4 T", "CD14+ Mono", "B", "CD8 T", "FCGR3A+ Mono",
#                      "NK", "DC", "Platelet")
# names(new.cluster.ids) <- levels(pbmc)
# pbmc <- RenameIdents(pbmc, new.cluster.ids)
# DimPlot(pbmc, reduction = "tsne", label = TRUE, pt.size = 0.5) + NoLegend()
