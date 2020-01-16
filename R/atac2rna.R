#' Title Convert a peak matrix to a gene activity matrix
#'
#' This function provide two methods to create gene acticity matrix. One is to take in a peak matrix and an annotation file (gtf) and
#' collapse the peak matrix to a gene activity matrix. And the other one is construct a feature x cell matrix from a genomic fragments file.
#'
#' @importFrom Seurat CreateSeuratObject
#' @importFrom Seurat CreateGeneActivityMatrix
#' @importFrom Signac NucleosomeSignal
#' @importFrom Signac SetFragments
#' @importFrom Signac Extend
#' @importFrom Signac FeatureMatrix
#' @importFrom Signac GRangesToString
#' @importFrom ensembldb genes
#' @importFrom GenomeInfoDb keepStandardChromosomes
#'
#' @param peaks Matrix of peak counts.
#' @param metadata Additional cell-level metadata to add to the Seurat object. Should be a data frame where the rows are cell names and the columns are additional metadata fields.
#' @param fragmentpath Path to tabix-indexed fragments file.
#' @param annotation.file Path to GTF annotation file.
#' @param qualitycontrol Whether to do quality control for peak matrix. If it is TRUE, it will compute QC Metrics and delete some cells. One needs both metadata data and fragments file to achieve quality control.
#' @param alpha Integer value 0 or 1 to decide which method to use. alpha = 0 is to use annotation file, and alpha = 1 is to use genomic fragments file.
#' @param seq.levels Which seqlevels to keep (corresponds to chromosomes usually).
#' @param include.body Include the gene body?
#' @param upstream Number of bases upstream to consider.
#' @param downstream Number of bases downstream to consider.
#' @param EnsDbobj For toSAF a GRangesList object. .
#' @param chunk Number of chunks to use when processing the fragments file. Fewer chunks may enable faster processing, but will use more memory.
#' @param filter A filter describing which results to retrieve from the database.
#'
#' @return A signcle cell RNA sequence matrix.
#' @export
#'
atac2rna <- function(peaks, metadata = NULL, fragmentpath = NULL, annotation.file = NULL, qualitycontrol = FALSE, alpha = 0,
                     seq.levels = c(1:22, "X", "Y"), include.body = TRUE, upstream = 2000, downstream = 0,
                     EnsDbobj = EnsDb.Hsapiens.v75, chunk = 50, filter = ~ gene_biotype == "protein_coding"){
  ## do quality control
  if (qualitycontrol = TRUE){

    if (is.null(metadata) = TRUE){
      stop(paste("We need metadata to do quality control!"))
    }

    object <- Seurat::CreateSeuratObject(
      counts = peaks,
      assay = 'peaks',
      project = 'ATAC',
      min.cells = 1,
      meta.data = metadata
    )

    if (is.null(fragmentpath) = TRUE){
      stop(paste("We need fragments file to do quality control!"))
    }

    fragment.path <- fragmentpath

    object <- Signac::SetFragments(
      object = object,
      file = fragment.path
    )
    # Computing QC Metrics
    object <- Signac::NucleosomeSignal(object = object)

    object$pct_reads_in_peaks <- object$peak_region_fragments / object$total * 100
    object$blacklist_ratio <- object$blacklist_region_fragments / object$peak_region_fragments

    object$nucleosome_group <- ifelse(object$nucleosome_signal > 10, 'NS > 10', 'NS < 10')

    #filter out cells
    object <- subset(object, subset = peak_region_fragments > 1000 & peak_region_fragments < 20000 & pct_reads_in_peaks > 15 & blacklist_ratio < 0.05 & nucleosome_signal < 10)
  }

  ## easy method
  if (alpha == 0){
    if (is.null(annotation.file) == TRUE){
      stop(paste("We need 'GTF' file for 'CreateGeneActivityMatrix' function!"))
    }

    peaks <- object@assays$peaks@counts
    activity.matrix <- Seurat::CreateGeneActivityMatrix(peak.matrix = peaks, annotation.file = annotation.file,
                                                seq.levels = seq.levels, upstream = upstream, downstream = downstream,
                                                include.body = include.body, verbose = FALSE)

    # return result
    return(activity.matrix)
  } else {
    # Signac method, which needs fragments file.

    if (is.null(fragmentpath) = TRUE){
      stop(paste("We need fragments file to use 'FeatureMatrix' function."))
    }

    fragment.path <- fragmentpath

    # extract gene coordinates from Ensembl, and ensure name formatting is consistent with Seurat object
    gene.coords <- ensembldb::genes(x = EnsDbobj, filter = filter)
    seqlevelsStyle(gene.coords) <- 'UCSC'
    genebody.coords <- GenomeInfoDb::keepStandardChromosomes(gene.coords, pruning.mode = 'coarse')
    genebodyandpromoter.coords <- Signac::Extend(x = gene.coords, upstream = 2000, downstream = 0)

    # create a gene by cell matrix
    gene.activities <- Signac::FeatureMatrix(
      fragments = fragment.path,
      features = genebodyandpromoter.coords,
      cells = colnames(object),
      chunk = chunk
    )

    # convert rownames from chromsomal coordinates into gene names
    gene.key <- genebodyandpromoter.coords$gene_name
    names(gene.key) <- Signac::GRangesToString(grange = genebodyandpromoter.coords)
    rownames(gene.activities) <- gene.key[rownames(gene.activities)]

    # return result
    return(gene.activities)
  }
}
