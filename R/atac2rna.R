#' Title Constructing gene activity score matrix from scATAC-seq
#'
#' This function constructs gene activity score matrix, which is just a sum of reads in peaks around each gene, from scATAC-seq data.
#'
#' @importFrom Seurat Read10X_h5
#' @importFrom Seurat CreateSeuratObject
#' @importFrom Signac SetFragments
#' @importFrom Signac NucleosomeSignal
#' @importFrom Signac RunTFIDF
#' @importFrom Signac Extend
#' @importFrom Signac FeatureMatrix
#' @importFrom Signac GRangesToString
#' @importFrom ensembldb genes
#' @importFrom ensembldb seqlevelsStyle
#' @importFrom GenomeInfoDb keepStandardChromosomes
#'
#' @param countmatrixfile A path to ATAC seq data, each row of the matrix represents a region of the genome (a ‘peak’), which must be a ".h5" file.
#' @param metadatafile A path to meta data for above matrix, which must be a ".csv" file.
#' @param fragment.path A path to fragments file, which must be in format "tsv.gz".This represents a full list of all unique fragments across all single cells.
#' @param mincells An integer value. When creating gene activity score matrix, only including features detected in at least this many cells.
#' @param minfeatures An integer value. When creating gene activity score matrix, only including cells where at least this many features are detected.
#'
#' @return A Seurat object which contains gene activity score matrix.
#' @export
#'
atac2rna <- function(countmatrixfile, metadatafile, fragment.path, mincells = 10, minfeatures = 1000){

  ## input data and set Seurat object
  counts <- Seurat::Read10X_h5(filename = countmatrixfile)

  metadata <- read.csv(
    file = metadatafile,
    header = TRUE,
    row.names = 1)

  obj <- Seurat::CreateSeuratObject(
    counts = counts,
    assay = 'peaks',
    project = 'ATAC',
    min.cells = 1,
    meta.data = metadata
  )

  obj <- Signac::SetFragments(
    object = obj,
    file = fragment.path
  )

  ## quality control
  obj <- Signac::NucleosomeSignal(object = obj)

  obj$pct_reads_in_peaks <- obj$peak_region_fragments / obj$total * 100
  obj$blacklist_ratio <- obj$blacklist_region_fragments / obj$peak_region_fragments

  obj$nucleosome_group <- ifelse(obj$nucleosome_signal > 10, 'NS > 10', 'NS < 10')

  #filter out cells
  obj <- subset(obj, subset = peak_region_fragments > 1000 & peak_region_fragments < 20000 &
                  pct_reads_in_peaks > 15 & blacklist_ratio < 0.05 & nucleosome_signal < 10)

  ## Normalization
  obj <- Signac::RunTFIDF(obj)

  ## Create a gene activity matrix

  # extract gene coordinates from Ensembl, and ensure name formatting is consistent with Seurat object
  gene.coords <- ensembldb::genes(EnsDb.Hsapiens.v75, filter = ~ gene_biotype == "protein_coding")
  ensembldb::seqlevelsStyle(gene.coords) <- 'UCSC'
  genebody.coords <- GenomeInfoDb::keepStandardChromosomes(gene.coords, pruning.mode = 'coarse')
  genebodyandpromoter.coords <- Signac::Extend(x = gene.coords, upstream = 2000, downstream = 0)

  # create a gene by cell matrix
  gene.activities <- Signac::FeatureMatrix(
    fragments = fragment.path,
    features = genebodyandpromoter.coords,
    cells = colnames(obj),
    chunk = 10
  )

  # convert rownames from chromsomal coordinates into gene names
  gene.key <- genebodyandpromoter.coords$gene_name
  names(gene.key) <- Signac::GRangesToString(grange = genebodyandpromoter.coords)
  rownames(gene.activities) <- gene.key[rownames(gene.activities)]

  # create scRNAseq Seurat object

  obj_return <- Seurat::CreateSeuratObject(counts = gene.activities, min.cells = mincells, min.features = minfeatures)

  ## return gene activity score matrix

  return(obj_return)
}
