# library(rio)
# library(devtools)
#
# #### first data set
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
# use_data(PanglaoDB_list)
#
# rm(list = ls())
#
# #### second data set
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
# use_data(GSEA_list)
