This package is used to label cell type for each cell in single-cell data. The first step is to use “Seurat” package to do data pre-process and cells cluster. After clustering the cell, our purpose is to label cell type for each cluster. The second step is to do use differential gene expression analysis to find marker genes for each cluster, that is to find which genes have higher expression in given cluster and lower expression in others. Then for each cluster we get a list of marker genes. It is also known that each cell type has its own marker genes, thus the final step is to use the Gene set enrichment analysis (GSEA) method to compare two list of genes. For each cluster, we compare its marker genes to all cell type marker genes. GSEA will give us which cell type marker genes are the most relevant to such cluster marker genes. That is the cell type for such cluster.This package aims to assemble all these steps to give cell type for each cluster.
