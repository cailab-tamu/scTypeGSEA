# scTypeGSEA

This package is used to assign cell type for each cell in single-cell data. The first step is to use the “Seurat” package to do data pre-process and cell cluster. After clustering the cell, our purpose is to assign cell types for each cluster. The second step is to use differential gene expression analysis to find full changes for each cluster. Then, for each cluster, we get a list of full changes. The final step is to use the Gene set enrichment analysis (GSEA) method to decide the cell type for each cluster. GSEA will give us the normalized enrichment score (NES) and NES shows which cell type is significant. We set each cluster the cell type with the highest NES.

## Installation

```{r}
library(devtools)
install_github("cailab-tamu/scTypeGSEA")
library(scTypeGSEA)
```

## Quick example

One can do quality control, data pre-process, cluster, get full changes, do GSEA and label the cell in one step.
```{r}
pbmc_res <- assignCellType(obj = pbmc_small, min.cells = 1, min.features = 10, nfeatures = 100, npcs = 10, dims = 1:10, k.param = 5, resolution = 0.75, min.pct = 0.25, test.use = "MAST", minSize = 5)
```

It will return a list. The first entry is a Seurat object.
```{r}
pbmc_res$obj
```

The second entry is a matrix contains cluster-ID, cell type and p-value for each cell.
```{r}
head(pbmc_res$cell_mat)
```




