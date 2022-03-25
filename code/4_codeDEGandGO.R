library(SingleCellExperiment)
library(scater)
library(scran)
library(Seurat)

setwd("F:/Device_D/SNU_Chicken/data/")
seurat_male <- readRDS("seurat_male_QC_CCcorrected.rds")
seurat_female <- readRDS("seurat_female_QC_CCcorrected.rds")

seurat_male@active.ident <- seurat_male$CellType
seurat_female@active.ident <- seurat_female$CellType

markers_male <- FindAllMarkers(seurat_male)
markers_male <- markers_male[markers_male$p_val_adj <= 0.05, ]
write.csv(markers_male, "E2.5/OutlierSplitted/DEG/Male_CellTypeDEGs_Final.csv")

markers_female <- FindAllMarkers(seurat_female)
markers_female <- markers_female[markers_female$p_val_adj <= 0.05, ]
write.csv(markers_female, "E2.5/OutlierSplitted/DEG/Female_CellTypeDEGs_Final.csv")



library(topGO)
library(org.Gg.eg.db)

markers <- read.csv("E2.5/OutlierSplitted/DEG/Male_CellTypeDEGs_Final.csv", row.names = 1)

backgroundGenes = rownames(seurat_male)
#Upreg
for(j in unique(seurat_male$CellType)){
  allmarkers <- markers[markers$cluster == j, ]
  rownames(allmarkers) <- allmarkers[,"gene"]
  
  targetGenes = subset(allmarkers, allmarkers[,"avg_logFC"] > 0.25 & allmarkers[,"p_val_adj"] < 0.05)$gene
  allGene = factor(as.integer(backgroundGenes %in% targetGenes))
  names(allGene) = backgroundGenes
  library(plyr)
  onts= c("BP")
  tab = as.list(onts)
  names(tab) = onts
  for (i in 1:length(onts)) {
    tgd = new("topGOdata", ontology = onts[i], allGenes=allGene, nodeSize=10,
              annot=annFUN.org, mapping="org.Gg.eg.db", ID = "symbol")
    resultTopGO.elim = runTest(tgd, algorithm="elim", statistic = "Fisher")
    tab[[i]] = GenTable(tgd, Fisher.elim = resultTopGO.elim,
                        orderBy="Fisher.classic", topNodes=200, numChar = 1000L)
  }
  
  topGOResults = rbind.fill(tab)
  write.csv(topGOResults, file=paste("E2.5/OutlierSplitted/GO/Male/GO_Upreg_Male_", j, ".csv", sep=""))
}


#Downreg
for(j in unique(seurat_male$CellType)){
  allmarkers <- markers[markers$cluster == j, ]
  rownames(allmarkers) <- allmarkers[,"gene"]
  
  targetGenes = subset(allmarkers, allmarkers[,"avg_logFC"] < -0.25 & allmarkers[,"p_val_adj"] < 0.05)$gene
  allGene = factor(as.integer(backgroundGenes %in% targetGenes))
  names(allGene) = backgroundGenes
  library(plyr)
  onts= c("BP")
  tab = as.list(onts)
  names(tab) = onts
  for (i in 1:length(onts)) {
    tgd = new("topGOdata", ontology = onts[i], allGenes=allGene, nodeSize=10,
              annot=annFUN.org, mapping="org.Gg.eg.db", ID = "symbol")
    resultTopGO.elim = runTest(tgd, algorithm="elim", statistic = "Fisher")
    tab[[i]] = GenTable(tgd, Fisher.elim = resultTopGO.elim,
                        orderBy="Fisher.classic", topNodes=200, numChar = 1000L)
  }
  
  topGOResults = rbind.fill(tab)
  write.csv(topGOResults, file=paste("E2.5/OutlierSplitted/GO/Male/GO_Downreg_Male_", j, ".csv", sep=""))
}


#Female
markers <- read.csv("E2.5/OutlierSplitted/DEG/Female_CellTypeDEGs_Final.csv", row.names = 1)

backgroundGenes = rownames(seurat_female)
#Upreg
for(j in unique(seurat_female$CellType)){
  allmarkers <- markers[markers$cluster == j, ]
  rownames(allmarkers) <- allmarkers[,"gene"]
  
  targetGenes = subset(allmarkers, allmarkers[,"avg_logFC"] > 0.25 & allmarkers[,"p_val_adj"] < 0.05)$gene
  allGene = factor(as.integer(backgroundGenes %in% targetGenes))
  names(allGene) = backgroundGenes
  library(plyr)
  onts= c("BP")
  tab = as.list(onts)
  names(tab) = onts
  for (i in 1:length(onts)) {
    tgd = new("topGOdata", ontology = onts[i], allGenes=allGene, nodeSize=10,
              annot=annFUN.org, mapping="org.Gg.eg.db", ID = "symbol")
    resultTopGO.elim = runTest(tgd, algorithm="elim", statistic = "Fisher")
    tab[[i]] = GenTable(tgd, Fisher.elim = resultTopGO.elim,
                        orderBy="Fisher.classic", topNodes=200, numChar = 1000L)
  }
  
  topGOResults = rbind.fill(tab)
  write.csv(topGOResults, file=paste("E2.5/OutlierSplitted/GO/Female/GO_Upreg_Female_", j, ".csv", sep=""))
}


#Downreg
for(j in unique(seurat_female$CellType)){
  allmarkers <- markers[markers$cluster == j, ]
  rownames(allmarkers) <- allmarkers[,"gene"]
  
  targetGenes = subset(allmarkers, allmarkers[,"avg_logFC"] < -0.25 & allmarkers[,"p_val_adj"] < 0.05)$gene
  allGene = factor(as.integer(backgroundGenes %in% targetGenes))
  names(allGene) = backgroundGenes
  library(plyr)
  onts= c("BP")
  tab = as.list(onts)
  names(tab) = onts
  for (i in 1:length(onts)) {
    tgd = new("topGOdata", ontology = onts[i], allGenes=allGene, nodeSize=10,
              annot=annFUN.org, mapping="org.Gg.eg.db", ID = "symbol")
    resultTopGO.elim = runTest(tgd, algorithm="elim", statistic = "Fisher")
    tab[[i]] = GenTable(tgd, Fisher.elim = resultTopGO.elim,
                        orderBy="Fisher.classic", topNodes=200, numChar = 1000L)
  }
  
  topGOResults = rbind.fill(tab)
  write.csv(topGOResults, file=paste("E2.5/OutlierSplitted/GO/Female/GO_Downreg_Female_", j, ".csv", sep=""))
}

