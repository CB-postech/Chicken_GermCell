library(SingleCellExperiment)
library(scater)
library(scran)
library(Seurat)
library(dplyr)
library(ggplot2)
library(zoo)
library(pheatmap)
library(viridis)
library(RColorBrewer)

smoothing_function <- function(seurat, expMat, genes, path,
                               z_cutoff = 3, windowsize=150){
  
  metadata = seurat[,seurat[[path]] == 1]@meta.data
  cells = rownames(metadata)
  metadata = metadata[order(metadata$pseudotime),]
  cells = as.vector(rownames(metadata))
  print(length(cells))
  
  mat = as.matrix(expMat)
  mat = mat[genes,cells]
  zmat = t(scale(t(mat)))
  TSmat = zmat
  TSmat[is.na(TSmat)] = 0
  TSmat[TSmat>z_cutoff] = z_cutoff
  TSmat[TSmat<-(z_cutoff)] =-(z_cutoff)
  
  for(i in 1:length(rownames(TSmat))){
    y = TSmat[rownames(TSmat)[i],]
    Y_hat = rollapply(y, windowsize, mean, align = 'center', fill = 'extend')
    TSmat[rownames(TSmat)[i],] = Y_hat
  }
  return(TSmat)
}

make_dynamic_heatmap <- function(seurat, expMat, genes, paths, pseudotime, repel = F, repelGenes = c(), windowsize){
  
  switch(as.character(length(paths)),
         "1" = {heatmap_matrix <- smoothing_function(seurat, expMat, genes, paths[1], windowsize = windowsize)},
         "2" = {mat1 = smoothing_function(seurat, expMat, genes, paths[1], windowsize = windowsize)
         mat2 = smoothing_function(seurat, expMat, genes, paths[2], windowsize = windowsize)
         heatmap_matrix <- cbind(mat1[,rev(colnames(mat1))], mat2)}
  )
  
  
  annotation_col <- data.frame(row.names = colnames(heatmap_matrix),
                               "pseudotime" = round(pseudotime[colnames(heatmap_matrix)], digit=2))
  annotation_colors = list("pseudotime" = colorRampPalette(brewer.pal(11, "YlGnBu"))(100))
  heatmap_matrix[is.na(heatmap_matrix)] = 0
  
  
  heatmap_matrix[heatmap_matrix > min(abs(max(heatmap_matrix)), abs(min(heatmap_matrix)))] <- min(abs(max(heatmap_matrix)), abs(min(heatmap_matrix)))
  heatmap_matrix[heatmap_matrix < -min(abs(max(heatmap_matrix)), abs(min(heatmap_matrix)))] <- -min(abs(max(heatmap_matrix)), abs(min(heatmap_matrix)))
  
  my_color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100)
  
  if(repel == F){
    p <- pheatmap(heatmap_matrix, cluster_cols=F, cluster_rows=T, show_colnames=F,
                  clustering_method = "ward.D2",
                  annotation_col = annotation_col,
                  annotation_colors = annotation_colors,
                  fontsize_row = 7, show_rownames = F,
                  treeheight_row = 0, cutree_row = 5,
                  color = my_color
    )
  }else{
    p <- pheatmap(heatmap_matrix, cluster_cols=F, cluster_rows=T, show_colnames=F,
                  clustering_method = "ward.D2",
                  annotation_col = annotation_col,
                  annotation_colors = annotation_colors,
                  fontsize_row = 7,
                  treeheight_row = 0,
                  color = my_color
                  
    )
    p <- add.flag(p, repelGenes, repel.degree = 1)
  }
  
  print(p)
  return(p)
}


#ortholog-extracted objects
gg_male <- readRDS("seurat_Gg_male.rds")
hs_male <- readRDS("seurat_Hs_male.rds")

hs_male <- ScaleData(hs_male, features = rownames(hs_male))

#Coerce all cells to path2 to fit the function
hs_male$Path2 <- rep(1, ncol(hs_male))
names(hs_male$Path2) <- colnames(hs_male)
gg_male$Path2 <- rep(1, ncol(gg_male))
names(gg_male$Path2) <- colnames(gg_male)

genes4use <- intersect(VariableFeatures(gg_male), VariableFeatures(hs_male))


# run
pseudotime <- gg_male$pseudotime
names(pseudotime) <- colnames(gg_male)

p_gg <- make_dynamic_heatmap(seurat = gg_male, expMat = gg_male@assays$RNA@scale.data[genes4use, ],
                         genes = genes4use, paths = c("Path2"), pseudotime = pseudotime, windowsize = 50, repel = F,
)



pseudotime <- hs_male$pseudotime
names(pseudotime) <- colnames(hs_male)

p_hs <- make_dynamic_heatmap(seurat = hs_male, expMat = hs_male@assays$RNA@scale.data[genes4use, ],
                            genes = genes4use, paths = c("Path2"), pseudotime = pseudotime, windowsize = 15, repel = F,
)

gene_clusters_gg <- cutree(p_gg$tree_row, k = 5)
gene_clusters_hs <- cutree(p_hs$tree_row, k = 5)

gc_df <- as.data.frame(gene_clusters_gg)
gc_df$gene_clusters_hs <- gene_clusters_hs[rownames(gc_df)]

#Construct sankey network
sankey_df <- data.frame("HM1" = rep(0, 5), "HM2" = rep(0, 5), "HM3" = rep(0, 5), "HM4" = rep(0, 5), "HM5" = rep(0, 5)) # rows: chickens', cols: humans'

for(i in c(1:5)){
  for(j in c(1:5)){
    for(k in c(1:nrow(gc_df))){
      if(as.vector(gc_df[k, ])[1] == i & as.vector(gc_df[k, ])[2] == j){
        sankey_df[i, j] = sankey_df[i, j] + 1
      }
    }
  }
}
rownames(sankey_df) <- paste0("CM", c(1:5))

library(tidyverse)
library(networkD3)

# Reshape the matrix data
sankey_df <- sankey_df %>%
  rownames_to_column %>%
  gather(key = 'key', value = 'value', -rowname) %>%
  filter(value > 0)
colnames(sankey_df) <- c("source", "target", "value")
sankey_df$target <- paste(sankey_df$target, " ", sep="")

# Create a node data frame
nodes <- data.frame(name=c(as.character(sankey_df$source),
                           as.character(sankey_df$target)) %>%
                      unique()
)

sankey_df$IDsource=match(sankey_df$source, nodes$name)-1 
sankey_df$IDtarget=match(sankey_df$target, nodes$name)-1

# Make the Network
# set "iterations=0" to avoid automatic assignment of the box order
sankeyNetwork(Links = sankey_df, Nodes = nodes,
              Source = "IDsource", Target = "IDtarget",
              Value = "value", NodeID = "name", 
              sinksRight=FALSE,
              #colourScale=ColourScal, 
              nodeWidth=20, fontSize=0, nodePadding=20,
              iterations=0
)


gene_clusters_gg[p_gg$tree_row$order]
gene_clusters_hs[p_hs$tree_row$order]

gc_df_tmp <- gc_df

gc_df$gene_clusters_gg[gc_df$gene_clusters_gg == 5] <- "A"
gc_df$gene_clusters_gg[gc_df$gene_clusters_gg == 3] <- "B"
gc_df$gene_clusters_gg[gc_df$gene_clusters_gg == 1] <- "C"
gc_df$gene_clusters_gg[gc_df$gene_clusters_gg == 4] <- "D"
gc_df$gene_clusters_gg[gc_df$gene_clusters_gg == 2] <- "E"
gc_df$gene_clusters_gg[gc_df$gene_clusters_gg == "A"] <- "M1"
gc_df$gene_clusters_gg[gc_df$gene_clusters_gg == "B"] <- "M2"
gc_df$gene_clusters_gg[gc_df$gene_clusters_gg == "C"] <- "M3"
gc_df$gene_clusters_gg[gc_df$gene_clusters_gg == "D"] <- "M4"
gc_df$gene_clusters_gg[gc_df$gene_clusters_gg == "E"] <- "M5"

gc_df$gene_clusters_hs[gc_df$gene_clusters_hs == 5] <- "a"
gc_df$gene_clusters_hs[gc_df$gene_clusters_hs == 3] <- "b"
gc_df$gene_clusters_hs[gc_df$gene_clusters_hs == 2] <- "c"
gc_df$gene_clusters_hs[gc_df$gene_clusters_hs == 1] <- "d"
gc_df$gene_clusters_hs[gc_df$gene_clusters_hs == 4] <- "e"
gc_df$gene_clusters_hs[gc_df$gene_clusters_hs == "a"] <- "m1"
gc_df$gene_clusters_hs[gc_df$gene_clusters_hs == "b"] <- "m2"
gc_df$gene_clusters_hs[gc_df$gene_clusters_hs == "c"] <- "m3"
gc_df$gene_clusters_hs[gc_df$gene_clusters_hs == "d"] <- "m4"
gc_df$gene_clusters_hs[gc_df$gene_clusters_hs == "e"] <- "m5"

gc_df <- gc_df[order(gc_df$gene_clusters_hs), ]
gc_df <- gc_df[order(gc_df$gene_clusters_gg), ]

# write.csv(gc_df, "Human_Chicken_GeneClusterInfo_Male.csv")