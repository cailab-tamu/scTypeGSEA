# scTypeGSEA

This package is designed to assign cell type labels for each identified cluster in single-cell data. This package uses the “Seurat” R package to do data pre-processing and cell clustering. After clustering the cells, we use differential gene expression analysis to compute the fold-change in gene expression by comparing the cluster profile against all the other identified clusters together. Then we use the Gene Set Enrichment Analysis (GSEA) technique to compute the enrichment (NES and their associated P-value) of marker genes defined for a set of cell types. GSEA analysis provides us the most statistically relevant cell type for each cluster that is finally assigned to the group.

## Installation

You can use following codes to install the package.

```{r}
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("fgsea")
BiocManager::install("MAST")
install.packages("Seurat")

library(devtools)
install_github("cailab-tamu/scTypeGSEA")
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

For more details about this pipeline, please read "vignettes".




