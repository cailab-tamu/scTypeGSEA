## ---- warning=FALSE, message=FALSE, eval=FALSE--------------------------------
#  # You may need following codes to install dependent packages.
#  if (!requireNamespace("BiocManager", quietly = TRUE))
#      install.packages("BiocManager")
#  
#  BiocManager::install("fgsea")
#  BiocManager::install("MAST")
#  install.packages("Seurat")
#  
#  library(devtools)
#  install_github("cailab-tamu/scTypeGSEA")

## -----------------------------------------------------------------------------
library(scTypeGSEA)

## ---- warning=FALSE-----------------------------------------------------------
pbmc_example_res <- assignCellType(small_RNA, min.cells = 1, min.features = 10, 
                                   nfeatures = 100, npcs = 10,
                                   dims = 1:10, k.param = 5, resolution = 0.75,
                                   min.pct = 0.25, test.use = "MAST", minSize = 5)

## -----------------------------------------------------------------------------
pbmc_example_res$Seurat_obj

## -----------------------------------------------------------------------------
head(pbmc_example_res$cell_mat)

## -----------------------------------------------------------------------------
head(pbmc_example_res$cluster_celltype)

## ---- warning=FALSE-----------------------------------------------------------
library(scTypeGSEA)
pbmc <- scqc(pbmc_raw, min.cells = 3, min.features = 200, percent.mt = 5, 
             normalization.method = "LogNormalize", scale.factor = 10000, 
            selection.method = "vst", nfeatures = 2000, npcs = 50)
pbmc

## ---- warning=FALSE-----------------------------------------------------------
set.seed(47)
pbmc <- doClustering(pbmc, dims = 1:50, k.param = 20, resolution = 0.5)
head(pbmc@meta.data$seurat_clusters)

## ---- eval=FALSE--------------------------------------------------------------
#  ## don't run
#  pbmc <- doClustering(pbmc, cluster_cell = Your_cluster_list_for_cell)

## ---- warning=FALSE-----------------------------------------------------------
cluster_list <- getFC(pbmc, min.pct = 0.25, test.use = "MAST")
head(cluster_list[[1]])

## ---- warning=FALSE-----------------------------------------------------------
cluster_celltype <- doGSEA(cluster_list = cluster_list, db = "PanglaoDB_list",
                            minSize = 15, maxSize = 500)
head(cluster_celltype)

## ---- eval=FALSE--------------------------------------------------------------
#  ## Don't run
#  cluster_celltype <- doGSEA(cluster_list = cluster_list, db = "path/to/your/rdsfile")

## ---- warning=FALSE-----------------------------------------------------------
pbmc_res <- labelCelltype(pbmc, cluster_celltype)

## -----------------------------------------------------------------------------
pbmc_res$obj

## -----------------------------------------------------------------------------
head(pbmc_res$cell_mat)

## ---- warning=FALSE, message=FALSE, fig.height=5------------------------------
library(Seurat)
pbmc <- pbmc_res$obj
pbmc <- Seurat::RunTSNE(pbmc, dims = 1:30)
Seurat::DimPlot(pbmc, reduction = "tsne", label = TRUE, pt.size = 0.5) + NoLegend()

## ---- warning=FALSE, message=FALSE, eval=FALSE--------------------------------
#  if (!requireNamespace("BiocManager", quietly = TRUE))
#      install.packages("BiocManager")
#  
#  BiocManager::install("GenomeInfoDb")
#  BiocManager::install("ensembldb")
#  BiocManager::install("EnsDb.Hsapiens.v75")
#  BiocManager::install("GenomicRanges")
#  
#  library(devtools)
#  install_github("timoast/signac")

## ---- message = FALSE, warning = FALSE----------------------------------------
annotation.file <- "~/Documents/Single cell/package example/R package/atac2rna/dataset/Homo_sapiens.GRCh37.82.gtf"
ATAC_example_res <- assignCellType(small_ATAC, datatype = "ATAC", annotation.file = annotation.file,
                                    min.cells = 1, min.features = 10, nfeatures = 100, npcs = 10,
                                    dims = 1:10, k.param = 5, resolution = 0.75,
                                    min.pct = 0.25, test.use = "MAST", minSize = 5)

## -----------------------------------------------------------------------------
ATAC_example_res$Seurat_obj
head(ATAC_example_res$cell_mat)
head(ATAC_example_res$cluster_celltype)

## ---- eval = FALSE------------------------------------------------------------
#  ## don't run
#  ATAC_example <- atac2rna(small_ATAC, annotation.file = annotation.file)
#  ATAC_example <- scqc(ATAC_example, min.cells = 1, min.features = 10, nfeatures = 100, npcs = 10)
#  ATAC_example <- doClustering(ATAC_example, dims = 1:10, k.param = 5, resolution = 0.75)
#  cluster_list <- getFC(ATAC_example, min.pct = 0.25, test.use = "MAST")
#  cluster_celltype <- doGSEA(cluster_list = cluster_list, minSize = 5)
#  ATAC_example_res <- labelCelltype(ATAC_example, cluster_celltype)

## ---- eval = FALSE------------------------------------------------------------
#  obj <- scqc(dta, datatype = "Antibody Capture", min.cells = 1, min.features = 10)

## ---- eval = FALSE------------------------------------------------------------
#  obj <- doClustering(obj, datatype = "Antibody Capture")

## ---- eval = FALSE------------------------------------------------------------
#  cluster_list <- getFC(obj, min.pct = 0.25, test.use = "wilcox")

## ---- eval = FALSE------------------------------------------------------------
#  cluster_celltype <- doGSEA(cluster_list = cluster_list, db = "path/to/your/rdsfile")

## ---- eval = FALSE------------------------------------------------------------
#  obj_res <- labelCelltype(obj, cluster_celltype)

