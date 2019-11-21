## ---- warning=FALSE, message=FALSE, eval=FALSE--------------------------------
#  ## You may need following codes to install dependent packages.
#  # if (!requireNamespace("BiocManager", quietly = TRUE))
#  #     install.packages("BiocManager")
#  #
#  # BiocManager::install("fgsea")
#  # BiocManager::install("MAST")
#  # install.packages("Seurat")
#  
#  library(devtools)
#  install_github("cailab-tamu/scTypeGSEA")

## -----------------------------------------------------------------------------
library(scTypeGSEA)

## ---- warning=FALSE, tidy = T, tidy.opts=list(width.cutoff = 35)--------------
pbmc_res <- assignCellType(obj = pbmc_small, min.cells = 1, min.features = 10, nfeatures = 100, 
                           npcs = 10, dims = 1:10, k.param = 5, resolution = 0.75, min.pct = 0.25, 
                           test.use = "MAST", minSize = 5)

## -----------------------------------------------------------------------------
pbmc_res$obj

## -----------------------------------------------------------------------------
head(pbmc_res$cell_mat)

## ---- warning=FALSE-----------------------------------------------------------
library(scTypeGSEA)
pbmc <- scqc(pbmc_raw, min.cells = 3, min.features = 200, percent.mt = 5)
pbmc

## ---- warning=FALSE, tidy = T, tidy.opts=list(width.cutoff = 35)--------------
set.seed(47)
pbmc <- doClustering(pbmc, normalization.method = "LogNormalize", scale.factor = 10000, 
                     selection.method = "vst", nfeatures = 2000, npcs = 50, dims = 1:50, 
                     k.param = 20, resolution = 0.5)
head(pbmc@meta.data$seurat_clusters)

## ---- eval=FALSE--------------------------------------------------------------
#  ## don't run
#  pbmc <- check_cluster(pbmc, k.param = 30, resolution = 0.75, doit = TRUE)

## ---- warning=FALSE-----------------------------------------------------------
cluster_list <- getFC(pbmc, min.pct = 0.25, test.use = "wilcox")
head(cluster_list[[1]])

## ---- warning=FALSE-----------------------------------------------------------
cluster_celltype <- doGSEA(cluster_list = cluster_list, db = "PanglaoDB_list",
                                  minSize = 15, maxSize = 500)
head(cluster_celltype)

## ---- eval=FALSE--------------------------------------------------------------
#  ## Don't run
#  cluster_celltype <- GSEA_analysis(cluster_list = cluster_list, otherdb = "path/to/your/rdsfile")

## ---- warning=FALSE-----------------------------------------------------------
pbmc_res <- labelSeurat(pbmc, cluster_celltype)
pbmc <- pbmc_res$obj
head(pbmc@active.ident)

## ---- warning=FALSE-----------------------------------------------------------
head(pbmc_res$cell_mat)

## ---- warning=FALSE, message=FALSE, fig.height=5------------------------------
library(Seurat)
pbmc <- pbmc_res$obj
pbmc <- Seurat::RunTSNE(pbmc, dims = 1:30)
Seurat::DimPlot(pbmc, reduction = "tsne", label = TRUE, pt.size = 0.5) + NoLegend()

