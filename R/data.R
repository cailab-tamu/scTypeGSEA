#' A small example version of the PBMC dataset
#'
#' A subsetted version of 10X Genomics' 3k PBMC dataset
#'
#' @format A txt file with 230 rows and 80 columns
#' @source \url{https://github.com/satijalab/seurat/blob/master/inst/extdata/pbmc_raw.txt}
"pbmc_raw"

#' A small example version of the PBMC dataset
#'
#' A subsetted version of 10X Genomics' 3k PBMC dataset
#'
#' @format A Seurat object with the following slots filled
#' \describe{
#'   \item{assays}{
#'   \itemize{Currently only contains one assay ("RNA" - scRNA-seq expression data)
#'   \item{counts - Raw expression data}
#'   \item{data - Normalized expression data}
#'   \item{scale.data - Scaled expression data}
#'   \item{var.features - names of the current features selected as variable}
#'   \item{meta.features - Assay level metadata such as mean and variance}
#'    }}
#'   \item{meta.data}{Cell level metadata}
#'   \item{active.assay}{Current default assay}
#'   \item{active.ident}{Current default idents}
#'   \item{graphs}{Neighbor graphs computed, currently stores the SNN}
#'   \item{reductions}{Dimensional reductions: currently PCA and tSNE}
#'   \item{version}{Seurat version used to create the object}
#'   \item{commands}{Command history}
#' }
"pbmc_test"

#' PBMC raw data set.
#'
#' A version of 10X Genomics' 3k PBMC dataset
#'
#' @format A txt file with 32738 rows and 2700 columns
"pbmc_raw"

#' Cell type gene expression markers from GSEA data base
#'
#' A list contains cell types with their marker gene
#'
#' @format A list with 257 cell types
#' @source \url{http://software.broadinstitute.org/gsea/msigdb/supplementary_genesets.jsp#SCSig}
"GSEA_list"

#' Cell type gene expression markers from GSEA data base
#'
#' A list contains cell types with their marker gene
#'
#' @format A list with 178 cell types
#' @source \url{https://panglaodb.se/markers.html}
"PanglaoDB_list"
