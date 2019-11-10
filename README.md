# scTypeGSEA

This package is used to label cell type for each cell in single-cell data. The first step is to use “Seurat” package to do data pre-process and cells cluster. After clustering the cell, our purpose is to label cell type for each cluster. The second step is to do use differential gene expression analysis to find marker genes for each cluster, that is to find which genes have higher expression in given cluster and lower expression in others. Then for each cluster we get a list of marker genes. It is also known that each cell type has its own marker genes, thus the final step is to use the Gene set enrichment analysis (GSEA) method to compare two list of genes. For each cluster, we compare its marker genes to all cell type marker genes. GSEA will give us which cell type marker genes are the most relevant to such cluster marker genes. That is the cell type for such cluster.This package aims to assemble all these steps to give cell type for each cluster.

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




