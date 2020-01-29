# scTypeGSEA

This package is designed to assign cell type labels for each identified cluster in single-cell data. This package uses the “Seurat” R package to do data pre-processing and cell clustering. After clustering the cells, we use differential gene expression analysis to compute the fold-change in gene expression by comparing the cluster profile against all the other identified clusters together. Then we use the Gene Set Enrichment Analysis (GSEA) technique to compute the enrichment (NES and their associated P-value) of marker genes defined for a set of cell types. GSEA analysis provides us with the most statistically relevant cell type for each cluster that is finally assigned to the group.

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
|atac2rna|Convert a peak matrix to a gene activity matrix|
|scqc|Performing single-cell data quality control|
|doClustering|Performing data process and clustering|
|getFC|Using differential gene expression analysis to compute the fold-change in gene expression|
|doGSEA|Doing gene set enrichment analysis(GSEA) for each cluster with its gene ranks|
|labelSeurat|Adding cell type to Seurat object|
|assignCellType|Doing quality control, data pre-process, cluster, get full changes and do GSEA to label the cell.|

## Quick example:

Here we use a toy data set "small_RNA" to show our main function "assignCellType" for single cell RNA sequence data. This function can achieve quality control, data pre-process, cluster, get fold changes, do GSEA and label the cell in one step.
```{r, tidy = TRUE, tidy.opts=list(width.cutoff = 50)}
pbmc_example_res <- assignCellType(small_RNA, min.cells = 1, min.features = 10, 
                                   nfeatures = 100, npcs = 10,
                                   dims = 1:10, k.param = 5, resolution = 0.75,
                                   min.pct = 0.25, test.use = "MAST", minSize = 5)
```

It will return a list with 4 slots. The first slot is a Seurat object.
```{r}
pbmc_example_res$Seurat_obj
```

```
An object of class Seurat 
230 features across 80 samples within 1 assay 
Active assay: RNA (230 features)
 1 dimensional reduction calculated: pca
```

The second slot is a dataframe including "cluster-ID", "cell type" and "p-value" for each cell.
```{r}
head(pbmc_example_res$cell_mat)
```

```
               ClusterID Cell Type  padj                  
ATGCCAGAACGACT "1"       "NK_cells" "7.18829936801286e-10"
CATGGCCTGTGCAT "1"       "NK_cells" "7.18829936801286e-10"
GAACCTGATGAACC "1"       "NK_cells" "7.18829936801286e-10"
TGACTGGATTCTCA "1"       "NK_cells" "7.18829936801286e-10"
AGTCAGACTGCACA "1"       "NK_cells" "7.18829936801286e-10"
TCTGATACACGTGT "1"       "NK_cells" "7.18829936801286e-10"
```

The third slot is "cluster list". For each cluster, there will be a ranked gene/feature list, which is used to do DSEA. The forth slot is dataframe including "cell type", "NES" and "padj" for each cluster.

```{r}
head(pbmc_example_res$cluster_celltype)
```

```
         cell type         NES                padj                  
Cluster0 "NK_cells"        "2.86436828618051" "7.18829936801286e-10"
Cluster1 "Dendritic_cells" "2.54600207137584" "6.306799398822e-06"  
```

Please read [vignettes](https://github.com/cailab-tamu/scTypeGSEA/blob/master/doc/Example_scTypeGSEA.pdf) for more details about how this function (pipeline) works and how to deal with other data type, for example, ATAC data. 
