#' Convert a peak matrix to a gene activity matrix
#'
#' This function takes in a peak matrix and an annotation file (gtf)and collapse the peak matrix to a gene activity matrix.
#' @importFrom Seurat CreateSeuratObject
#' @importFrom Seurat CreateGeneActivityMatrix
#' @importFrom Seurat CreateAssayObject
#' @importFrom Seurat DefaultAssay
#'
#' @param peaks Matrix of peak counts.
#' @param annotation.file Path to GTF annotation file.
#' @param seq.levels Which seqlevels to keep (corresponds to chromosomes usually).
#' @param include.body Include the gene body?
#' @param upstream Number of bases upstream to consider.
#' @param downstream Number of bases downstream to consider.
#'
#' @return A Seurat object with both ATAC data matrix and gene activity matrix.
#'
#' @export
#'
#' @examples
#' \dontrun{annotation.file <- "../data/Homo_sapiens.GRCh37.82.gtf"
#'          ATAC_example <- atac2rna(pbmc_small_atac, annotation.file = annotation.file)}
#'
atac2rna <- function(peaks, annotation.file = NULL, seq.levels = c(1:22, "X", "Y"),
                     include.body = TRUE, upstream = 2000, downstream = 0){
  ## create Seurat Object
  object <- Seurat::CreateSeuratObject(
    counts = peaks,
    assay = 'peaks',
    project = 'ATAC',
    min.cells = 1
  )

  if (is.null(annotation.file) == TRUE){
    stop(paste("We need 'GTF' file for 'CreateGeneActivityMatrix' function!"))
  }

  activity.matrix <- Seurat::CreateGeneActivityMatrix(peak.matrix = peaks, annotation.file = annotation.file,
                                              seq.levels = seq.levels, upstream = upstream, downstream = downstream,
                                                include.body = include.body, verbose = FALSE)

  object[['RNA']] <- Seurat::CreateAssayObject(counts = activity.matrix)
  Seurat::DefaultAssay(object) <- 'RNA'
  # return result
  return(object)
}
