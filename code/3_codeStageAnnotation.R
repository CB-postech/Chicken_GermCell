library(SingleCellExperiment)
library(scater)
library(scran)
library(Seurat)
library(org.Gg.eg.db)

Cluster2CellType <- function(SeuratObj, ClusterSlot = "seurat_clusters", # name of slot as input
                             CellTypeSlot = "CellType", # name of cell type slot to be returned
                             QueryList = list(), # list of (name of cell type) to (name of clusters)
                             ReturnPlot = T){

  celltype_vec <- c(rep(NA, ncol(SeuratObj)))
  names(celltype_vec) <- colnames(SeuratObj)
  
  for(i in names(QueryList)){
    for(j in QueryList[[i]]){
      celltype_vec[SeuratObj[[ClusterSlot]] == j] <- i
    }
  }
  
  SeuratObj[[CellTypeSlot]] <- factor(celltype_vec, levels = names(QueryList))
  
  if(ncol(SeuratObj) == sum(table(celltype_vec))){
    print(table(celltype_vec))
    cat(paste0("Total #Cells: ", sum(table(celltype_vec)), "\n"))
  }else{
    stop(paste0("Total #Cells Differ: ", ncol(SeuratObj), ", ", sum(table(celltype_vec)), "\n"))
  }
  
  remove(celltype_vec)
  
  if(ReturnPlot == T){
    p <- DimPlot(SeuratObj, group.by = CellTypeSlot, label = T, label.size = 7, pt.size = 1) +
      theme(axis.line=element_blank(),
            axis.text.x=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks=element_blank(),
            axis.title.x=element_blank(),
            axis.title.y=element_blank(),
            panel.background=element_blank(),
            panel.border=element_blank(),
            panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            plot.background=element_blank()
      )
    print(p)
  }
  return(SeuratObj)
}


seurat_female <- readRDS("seurat_female_QC_CCcorrected.rds")
seurat_male <- readRDS("seurat_male_QC_CCcorrected.rds")

#load kegg
mapped <- mappedkeys(org.Gg.egPATH2EG)
L <- as.list(org.Gg.egPATH2EG[mapped])
Kegg_ID <- names(L)

library(biomaRt)
date = Sys.Date()
ensembl <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset="ggallus_gene_ensembl", host="www.ensembl.org")
ensemblGenes <- getBM(attributes=c('ensembl_gene_id', 'external_gene_name',  'entrezgene_id', 'gene_biotype'), mart=ensembl)
ensemblGenes[ensemblGenes$external_gene_name == "", "external_gene_name"] <- ensemblGenes[ensemblGenes$external_gene_name == "", "ensembl_gene_id"]


#assign stages: female
seurat_female@assays$RNA@var.features <- list(hvg_750 = seurat_female@assays$RNA@var.features)
for(j in c(1000, 2000, 3000, 4000, 5000)){ #grouping clusters using score of KEGG path genes from N hvgs
  sce <- as.SingleCellExperiment(seurat_female)
  
  #HVG selection
  dec <- modelGeneVar(sce)
  top.hvgs <- getTopHVGs(dec, n = j)
  
  seurat_female@assays$RNA@var.features[[paste0("hvg_", j)]] <- top.hvgs
  
  features <- top.hvgs
  
  #KEGG path/genes
  L_sym_fem <- list()
  for(i in names(L)){
    tmp_df <- ensemblGenes[ensemblGenes$entrezgene_id %in% L[[as.character(i)]],]
    L_sym_fem[[as.character(i)]] <- tmp_df$external_gene_name
    L_sym_fem[[as.character(i)]] <- L_sym_fem[[as.character(i)]][L_sym_fem[[as.character(i)]] %in% features]
  }
  
  library(KEGGREST)
  for(i in c(1:length(L_sym_fem))){
    names(L_sym_fem)[i] <- paste0("gga", names(L_sym_fem)[i])
    try(names(L_sym_fem)[i] <- sapply(keggGet(names(L_sym_fem)[i]), "[[", "NAME"))
  }
  
  #scoring using KEGG path genes
  tmp_seurat <- ScaleData(seurat_female, rownames(seurat_female))
  tmp_seurat@active.ident <- tmp_seurat$seurat_clusters
  avgExprs_fem <- AverageExpression(tmp_seurat, features = rownames(tmp_seurat), slot = "scale.data")
  
  avg_act <- as.data.frame(colMeans(avgExprs_fem$RNA[L_sym_fem[[as.character(names(L_sym_fem)[1])]], ]))
  colnames(avg_act) = names(L_sym_fem)[1]
  for(i in c(2:length(L_sym_fem))){
    tmp <- as.data.frame(colMeans(avgExprs_fem$RNA[L_sym_fem[[as.character(names(L_sym_fem)[i])]], ]))
    if(length(tmp) != 0){
      avg_act <- cbind(avg_act, tmp)
      
      colnames(avg_act)[i] <- names(L_sym_fem)[i]
    } else{
      print(paste0("# Genes for KEGG # ", names(L_sym_fem)[i], " is 0."))
    }
  }
  
  avg_act <- t(avg_act)
  avg_act[avg_act > 1] <- 1
  avg_act[avg_act < -1] <- -1
  
  
  library(RColorBrewer)
  palette_length = 100
  my_color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(palette_length)
  
  
  avg_act <- avg_act[, as.character(c(14, 2, 12, 6, 0, 1, 11, 15, 5, 10, 9, 8, 7, 13, 3, 4))] #reordering clusters
  avg_act <- avg_act[rowSums(is.nan(avg_act)) == 0,]

  
  my_breaks <- c(seq(min(avg_act), 0,
                     length.out=ceiling(palette_length/2) + 1),
                 seq(max(avg_act)/palette_length,
                     max(avg_act),
                     length.out=floor(palette_length/2)))
  
  pheatmap::pheatmap(avg_act, show_rownames = T, cluster_cols = T,
                     color=my_color,
                     breaks = my_breaks,
                     #treeheight_col = 0,
                     border_color = NA, width = 15, height = 15) # Figure S1E
}




#Male
seurat_male@assays$RNA@var.features <- list(hvg_434 = seurat_male@assays$RNA@var.features)
for(j in c(1000, 2000, 3000, 4000, 5000)){
  sce <- as.SingleCellExperiment(seurat_male)
  
  #HVG selection
  dec <- modelGeneVar(sce)
  top.hvgs <- getTopHVGs(dec, n = j)
  
  seurat_male@assays$RNA@var.features[[paste0("hvg_", j)]] <- top.hvgs
  
  features <- top.hvgs
  L_sym_m <- list()
  for(i in names(L)){
    tmp_df <- ensemblGenes[ensemblGenes$entrezgene_id %in% L[[as.character(i)]],]
    L_sym_m[[as.character(i)]] <- tmp_df$external_gene_name
    L_sym_m[[as.character(i)]] <- L_sym_m[[as.character(i)]][L_sym_m[[as.character(i)]] %in% features]
  }
  
  library(KEGGREST)
  for(i in c(1:length(L_sym_m))){
    names(L_sym_m)[i] <- paste0("gga", names(L_sym_m)[i])
    try(names(L_sym_m)[i] <- sapply(keggGet(names(L_sym_m)[i]), "[[", "NAME"))
  }
  
  tmp_seurat <- ScaleData(seurat_male, rownames(seurat_male))
  tmp_seurat@active.ident <- tmp_seurat$seurat_clusters
  avgExprs_m <- AverageExpression(tmp_seurat, features = rownames(tmp_seurat), slot = "scale.data")
  
  avg_act_m <- as.data.frame(colMeans(avgExprs_m$RNA[L_sym_m[[as.character(names(L_sym_m)[1])]], ]))
  colnames(avg_act_m) = names(L_sym_m)[1]
  for(i in c(2:length(L_sym_m))){
    tmp <- as.data.frame(colMeans(avgExprs_m$RNA[L_sym_m[[as.character(names(L_sym_m)[i])]], ]))
    if(length(tmp) != 0){
      avg_act_m <- cbind(avg_act_m, tmp)
      
      colnames(avg_act_m)[i] <- names(L_sym_m)[i]
    } else{
      print(paste0("# Genes for KEGG # ", names(L_sym_m)[i], " is 0."))
    }
  }
  
  
  avg_act_m[avg_act_m > 1] <- 1
  avg_act_m[avg_act_m < -1] <- -1
  
  avg_act_m <- t(avg_act_m)
  avg_act_m <- avg_act_m[, as.character(c(11, 3, 6, 7, 9, 10, 0, 8, 2, 1, 5, 4))] #reordering clusters
  avg_act_m <- avg_act_m[rowSums(is.nan(avg_act_m)) == 0,]
  
  library(RColorBrewer)
  palette_length = 100
  my_color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(palette_length)
  
  my_breaks <- c(seq(min(avg_act_m), 0,
                     length.out=ceiling(palette_length/2) + 1),
                 seq(max(avg_act_m)/palette_length,
                     max(avg_act_m),
                     length.out=floor(palette_length/2)))
  
  pheatmap::pheatmap(avg_act_m, show_rownames = T, cluster_cols = T,
                     color=my_color,
                     breaks = my_breaks,
                     #treeheight_col = 0,
                     border_color = NA, width = 15, height = 15) # Figure S1E
}




query_female <- list("S1" = c(2, 6, 12, 14), "S2" = c(0, 1, 11, 15), "S3" = c(5, 9 ,10), "S4" = c(3, 7, 8, 13), "S5" = c(4))
seurat_female <- Cluster2CellType(SeuratObj = seurat_female, QueryList = query_female, ReturnPlot = T)

query_male <- list("S1" = c(11, 3), "S2" = c(6, 7, 9, 0, 8, 10), "S3" = c(1, 2, 5), "S4" = c(4))
seurat_male <- Cluster2CellType(SeuratObj = seurat_male, QueryList = query_male, ReturnPlot = T)

# saveRDS(seurat_female, "seurat_female_CellTypeAnnotated.rds")
# saveRDS(seurat_male, "seurat_male_CellTypeAnnotated.rds")
