# scTypeGSEA

This package is used to label cell type for each cell in single-cell data. The first step is to use “Seurat” package to do data pre-process and cells cluster. After clustering the cell, our purpose is to label cell type for each cluster. The second step is to do use differential gene expression analysis to find marker genes for each cluster, that is to find which genes have higher expression in given cluster and lower expression in others. Then for each cluster we get a list of marker genes. It is also known that each cell type has its own marker genes, thus the final step is to use the Gene set enrichment analysis (GSEA) method to compare two list of genes. For each cluster, we compare its marker genes to all cell type marker genes. GSEA will give us which cell type marker genes are the most relevant to such cluster marker genes. That is the cell type for such cluster.This package aims to assemble all these steps to give cell type for each cluster.

## Installation

```{r}
library(devtools)
install_github("cailab-tamu/scTypeGSEA")
```

## Small Example
Next to use small data set "pbmc_small" quickly go through this pipline. First is to load data set. If it is a gene expression count matrix, you can use "check_seurat" function to create a Seurat object.
```{r}
pbmc_example <- check_seurat(pbmc_raw, min.cells = 1, min.features = 10)
```
Next to do data pre-process and cluster.
```{r}
pbmc_example <- check_cluster(pbmc_test, nfeatures = 100, npcs = 10, dims = 1:10, k.param = 5, resolution = 0.5)
```

To get rank of gene list, we need to do gene differential expression analysis.
```{r}
cluster_list <- Test_DE_cluster(pbmc_example, min.pct = 0.25, test.use = "MAST")
```

After getting rank of gene list for each cluster, we can use GSEA method to find the cell type for each cluster.
```{r}
cluster_celltype <- GSEA_analysis(cluster_list = cluster_list, minSize = 5, nperm = 1000)
```

The last step is to add label cell type to our Seurat object.
```{r}
pbmc_example <- addcelltype_Seurat(pbmc_example, cluster_celltype)
```


