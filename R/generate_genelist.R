# library(rio)
#
# #### first data set
# dta1 <- import("~/Documents/single cell/data set/GSEA/PanglaoDB_markers_17_Oct_2019.tsv")
# celltype <- dta1$`cell type`
# celltype <- unique(celltype)
# mylist.names <- celltype
# mylist <- vector("list", length(mylist.names))
# names(mylist) <- mylist.names
#
# for (i in 1:length(mylist.names)){
#   index <- which(dta1$`cell type` == celltype[i])
#   mylist[[i]] <- dta1$`official gene symbol`[index]
# }
#
# saveRDS(mylist, file = "~/Documents/single cell/data set/GSEA/PanglaoDB_list.rds")
#
# rm(list = ls())
#
# #### second data set
# dta2_1 <- import("~/Documents/single cell/data set/GSEA/scsig.v1.0.metadata.xls")
# celltype <- dta2_1$`Gene Set Standard Name`
# celltype <- unique(celltype)
# mylist.names <- celltype
# mylist <- vector("list", length(mylist.names))
# names(mylist) <- mylist.names
#
# for (i in 1:length(mylist.names)){
#   entry <- dta2_1[i, 12]
#   entry <- strsplit(entry, ",")
#   mylist[[i]] <- entry
# }
#
# saveRDS(mylist, file = "~/Documents/single cell/data set/GSEA/scsig_v1.0_metadata_xls.rds")
#
# rm(list = ls())
#
# #### the third data set
# dta2_2 <- import("~/Documents/single cell/data set/GSEA/scsig.v1.0.metadata.txt", header = FALSE) ## same data set as dta2_1
#
# rm(list = ls())
