GSEA_analysis <- function(obj, cluster_list, db){

  ncluster <- length(cluster_list)
  cluster_celltype <- rep(NA, ncluster)
  names(cluster_celltype) <- names(cluster_list)

  if (db == "PanglaoDB_list"){
    data_base <- readRDS("inst/PanglaoDB_list.rds")
  } else{
    dta_base <- readRDS("inst/GSEA_list.rds")
  }

  for (i in 1:ncluster){
    Ranks <- cluster_list[[i]]

    fgseaRes <- fgsea(pathways = data_base,
                      stats = Ranks,
                      minSize = 15,
                      maxSize = 500,
                      nperm = 10000)

    fgseaRes <- fgseaRes[order(fgseaRes[, "NES"], -fgseaRes[, "padj"], decreasing = TRUE), ]

    if (fgseaRes[, "padj"][1] > 0.05){
      cluster_celltype[i] <- "unidentified"
    } else{
      cluster_celltype[i] <- fgseaRes[, "pathway"][1]
    }
  }

  return(cluster_list)
}
