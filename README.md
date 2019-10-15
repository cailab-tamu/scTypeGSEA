This package is used to label the cell type. The input should be a Seurat object. (I may achieve the function that the input can also be a gene expression matrix with clustering label for cells later.)
Step1: For each gene in each cluster, we compared its gene expression with that of other clusters and selected the highly expressed genes. Thus, each cluster will have a list of genes.

Step2: For each list of genes, we compared it with the database. For each cell type is the database, we compute the number of intersecting genes between the highly expressed genes regarding that cell type and this list of genes.

Step3: For each list of genes, we set the cell type as the label who has the largest number of intersections. 
