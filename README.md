# scTypeGSEA

This package is designed to assign cell type labels for each identified cluster in single-cell data. This package uses the “Seurat” R package to do data pre-processing and cell clustering. After clustering the cells, we use differential gene expression analysis to compute the fold-change in gene expression by comparing the cluster profile against all the other identified clusters together. Then we use the Gene Set Enrichment Analysis (GSEA) technique to compute the enrichment (NES and their associated P-value) of marker genes defined for a set of cell types. GSEA analysis provides us the most statistically relevant cell type for each cluster that is finally assigned to the group.

## Installation:

You can use following codes to install the package. To use some functions in "fgsea", your R version must be >= 3.6.0.

```{r}
## You may need following codes to install dependent packages.

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("fgsea")
BiocManager::install("MAST")

library(devtools)
install_github("cailab-tamu/scTypeGSEA")
library(scTypeGSEA)
```

## Available functions:

|Code| Function |
|:-|:-|
|scqc|Performing single-cell data quality control|
|doClustering|Performing data process and clustering|
|getFC|Using differential gene expression analysis to compute the fold-change in gene expression|
|doGSEA|Doing gene set enrichment analysis(GSEA) for each cluster with its gene ranks|
|labelSeurat|Adding cell type to Seurat object|
|assignCellType|Doing quality control, data pre-process, cluster, get full changes and do GSEA to label the cell.|

## Quick example:

Here we use a toy data set "pbmc_small" to show our main function "assignCellType". One can do quality control, data pre-process, cluster, get fold changes, do GSEA and label the cell in one step.
```{r}
pbmc_res <- assignCellType(obj = pbmc_small, min.cells = 1, min.features = 10, nfeatures = 100, npcs = 10, dims = 1:10, k.param = 5, resolution = 0.75, min.pct = 0.25, test.use = "MAST", minSize = 5)
```

It will return a list with 2 slots. The first slot is a Seurat object.
```{r}
pbmc_res$obj
```

```
An object of class Seurat 
230 features across 80 samples within 1 assay 
Active assay: RNA (230 features)
 2 dimensional reductions calculated: pca, tsne
```

The second slot is a matrix contains cluster-ID, cell type and p-value for each cell.
```{r}
head(pbmc_res$cell_mat)
```

```
               ClusterID Cell Type padj                 
ATGCCAGAACGACT "1"       "T_cells" "1.5163313773336e-08"
CATGGCCTGTGCAT "1"       "T_cells" "1.5163313773336e-08"
GAACCTGATGAACC "1"       "T_cells" "1.5163313773336e-08"
TGACTGGATTCTCA "1"       "T_cells" "1.5163313773336e-08"
AGTCAGACTGCACA "1"       "T_cells" "1.5163313773336e-08"
TCTGATACACGTGT "1"       "T_cells" "1.5163313773336e-08"
```
For more details about this pipeline, please read [vignettes](https://github.com/cailab-tamu/scTypeGSEA/blob/master/inst/doc/Example_scTypeGSEA.pdf).




