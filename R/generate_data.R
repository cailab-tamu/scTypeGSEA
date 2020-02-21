# library(rio)
# library(devtools)
# library(Seurat)
# library(GenomeInfoDb)
# library(GenomicRanges)
# library(readxl)
#
# #################################
# ######## PanglaoDB_list #########
# #################################
# rm(list = ls())
# dta1 <- import("~/Documents/single cell/data set/GSEA/PanglaoDB_markers_17_Oct_2019.tsv")
# celltype <- dta1$`cell type`
# celltype <- unique(celltype)
# PanglaoDB_list.names <- celltype
# PanglaoDB_list.names <- gsub(" ", "_", PanglaoDB_list.names)
#
# PanglaoDB_list <- vector("list", length(PanglaoDB_list.names))
# names(PanglaoDB_list) <- PanglaoDB_list.names
#
# for (i in 1:length(PanglaoDB_list.names)) {
#   index <- which(dta1$`cell type` == celltype[i])
#   tmp <- dta1[, 2][index]
#   PanglaoDB_list[[i]] <- dta1[, 2][index]
# }
#
# # saveRDS(PanglaoDB_list, file = "inst/extdata/PanglaoDB_list.rds")
# load("~/Documents/PanglaoDB_list.rda")
# use_data(PanglaoDB_list, internal = TRUE)
#
# ##########################
# ### generate TAIR_list ###
# ##########################
# rm(list = ls())
#
# TAIR_dta <- read_xlsx("~/Downloads/3463_marker_genes.xlsx")
# TAIR_dta <- TAIR_dta[, -6]
# TAIR_dta <- as.data.frame(TAIR_dta)
#
# celltype <- colnames(TAIR_dta)
# celltype <- unique(celltype)
#
# TAIR_list <- vector("list", length(celltype))
# names(TAIR_list) <- celltype
# TAIR_dta <- rbind(TAIR_dta, NA)
#
# for (i in 1:length(celltype)) {
#   tmp <- TAIR_dta[, i]
#   index <- which(is.na(tmp) == TRUE)
#   TAIR_list[[i]] <- as.vector(TAIR_dta[1:(index[1] - 1), i])
# }
#
# # saveRDS(TAIR_list, file = "TAIR_list.rds")
# load("~/Documents//TAIR_list.rda")
# use_data(TAIR_list, internal = TRUE)
#
#
#
# ###########################
# #### GSEA_list ############
# ###########################
# rm(list = ls())
#
# dta2 <- import("~/Documents/single cell/data set/GSEA/scsig.v1.0.metadata.xls")
# celltype <- dta2$`Gene Set Standard Name`
# celltype <- unique(celltype)
# GSEA_list.names <- celltype
# GSEA_list <- vector("list", length(GSEA_list.names))
# names(GSEA_list) <- GSEA_list.names
#
# for (i in 1:length(GSEA_list.names)) {
#   entry <- dta2[i, 12]
#   entry <- strsplit(entry, ",")
#   GSEA_list[[i]] <- entry[[1]]
# }
#
# # saveRDS(GSEA_list, file = "inst/extdata/GSEA_list.rds")
# load("~/Documents//GSEA_list.rda")
# use_data(GSEA_list, internal = TRUE)
#
#
# use_data(PanglaoDB_list, GSEA_list, TAIR_list, internal = TRUE)
#
# ##########################
# ####### pbmc_raw #########
# ##########################
#
# rm(list = ls())
# pbmc_raw <- Read10X(data.dir = "~/Documents/single cell/package example/R package/Seurat/cluster/data set")
# use_data(pbmc_raw)
#
# #########################
# ###### small_RNA ########
# #########################
#
# rm(list = ls())
# small_RNA <- pbmc_small@assays$RNA@counts
# use_data(small_RNA)
#
# ########################
# ##### small_ATAC #######
# ########################
# rm(list = ls())
#
# peaks <- Read10X_h5("~/Documents/Single cell/package example/R package/atac2rna/dataset/atac_v1_pbmc_5k_filtered_peak_bc_matrix.h5")
#
# # Match peak and its gene
# peak.df <- rownames(x = peaks)
# peak.df <- do.call(what = rbind, args = strsplit(x = gsub(peak.df, pattern = ":", replacement = "-"), split = "-"))
# peak.df <- as.data.frame(x = peak.df)
# colnames(x = peak.df) <- c("chromosome", 'start', 'end')
# peaks.gr <- GenomicRanges::makeGRangesFromDataFrame(df = peak.df)
# BiocGenerics::start(peaks.gr[BiocGenerics::start(peaks.gr) == 0, ]) <- 1
# # get annotation file, select genes
# gtf <- rtracklayer::import(con = "~/Documents/Single cell/package example/R package/atac2rna/dataset/Homo_sapiens.GRCh37.82.gtf")
# gtf <- GenomeInfoDb::keepSeqlevels(x = gtf, value = c(1:22, "X", "Y"), pruning.mode = 'coarse')
# # change seqlevelsStyle if not the same
# if (!any(GenomeInfoDb::seqlevelsStyle(x = gtf) == GenomeInfoDb::seqlevelsStyle(x = peaks.gr))) {
#   GenomeInfoDb::seqlevelsStyle(gtf) <- GenomeInfoDb::seqlevelsStyle(peaks.gr)
# }
# gtf.genes <- gtf[gtf$type == 'gene']
# # Extend definition up/downstream
# upstream <- 2000
# downstream <- 0
# include.body <- TRUE
# if (include.body) {
#   gtf.body_prom <- Signac::Extend(x = gtf.genes, upstream = upstream, downstream = downstream)
# } else {
#   gtf.body_prom <- SummarizedExperiment::promoters(x = gtf.genes, upstream = upstream, downstream = downstream)
# }
#
# gene.distances <- GenomicRanges::distanceToNearest(x = peaks.gr, subject = gtf.body_prom)
# keep.overlaps <- gene.distances[rtracklayer::mcols(x = gene.distances)$distance == 0]
# peak.ids <- peaks.gr[S4Vectors::queryHits(x = keep.overlaps)]
# gene.ids <- gtf.genes[S4Vectors::subjectHits(x = keep.overlaps)]
#
# # gene names
# gene_names <- rownames(pbmc_small)
#
# # select rows of peaks
# index_ref <- queryHits(x = keep.overlaps)
# index_row <- NULL
# for (i in 1:length(gene_names)){
#   index_id <- which(gene.ids$gene_name == gene_names[i])
#   for (j in 1:length(index_id)){
#     index_row <- c(index_row, index_ref[index_id[j]])
#   }
# }
# index_row <- na.omit(index_row)
#
# # generate small peak matrix
# small_ATAC <- peaks[index_row, ] # 685 x 4654
# set.seed(1234)
# small_ATAC <- small_ATAC[, sample(1:4654, 80)] # 685 x 80
# small_ATAC <- small_ATAC[order(rownames(small_ATAC)), ]
# use_data(small_ATAC)
