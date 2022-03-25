library(SingleCellExperiment)
library(scater)
library(scran)
library(Seurat)


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



hs_male <- readRDS("SeuratObj_HumanFGC_Male_logNorm.rds")
hs_female <- readRDS("SeuratObj_HumanFGC_Female_logNorm.rds")


plot(hs_male@reductions$pca@stdev)
PCA=20

set.seed(10)
hs_male <- FindNeighbors(hs_male, reduction = "pca", dims = 1:PCA)
hs_male <- FindClusters(hs_male, resolution = 0.8)
hs_male <- RunUMAP(hs_male, reduction = "pca", dims = 1:PCA)
hs_male <- RunTSNE(hs_male, reduction = "pca", dims = 1:PCA, check_duplicates = FALSE)

tmp <- CellCycleScoring(hs_male, s.features = cc.genes.updated.2019$s.genes, g2m.features = cc.genes.updated.2019$g2m.genes, set.ident = T)
tmp$CellCycle <- factor(as.vector(tmp@active.ident), levels = c("G1", "G2M", "S"))

markers <- list("Male-S1" = c("POU5F1", "NANOG", "SALL4", "TFAP2C", "DPPA3"),
                "Male-S2" = c("SOX13", "DAZL"),
                "Male-S3" = c("DDX4", "TDRD6", "PIWIL2", "SIX1", "ZGLP1", "SYCP1", "TEX12", "SPO11")
)
markers_vec <- unlist(markers)
markers_vec <- markers_vec[markers_vec %in% rownames(hs_male@assays$RNA@data)]
length(markers_vec)

hs_male@active.ident <- hs_male$seurat_clusters
avgExprs <- AverageExpression(hs_male, features = markers_vec, assays = "RNA", slot = "data")
df <- data.frame(Markers = c(rep(names(markers)[1], length(markers[[1]])),
                             rep(names(markers)[2], length(markers[[2]])),
                             rep(names(markers)[3], length(markers[[3]]))
))
rownames(df) <- unlist(markers)

cols <- gg_color_hue(3)
cour<-list(
  Markers=cols
)

scaledExprs <- t(scale(t(avgExprs$RNA)))
scaledExprs[scaledExprs > -min(scaledExprs)] <- -min(scaledExprs)

library(RColorBrewer)
palette_length = 100
my_color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(palette_length)

my_breaks <- c(seq(min(scaledExprs), 0,
                   length.out=ceiling(palette_length/2) + 1),
               seq(max(scaledExprs)/palette_length,
                   max(scaledExprs),
                   length.out=floor(palette_length/2)))

p <- pheatmap::pheatmap(
  scaledExprs,
  cluster_cols = T, cluster_rows = F,
  annotation_row = df, annotation_colors = list(cour$Markers),
  breaks = my_breaks, color=my_color,
  labels_row = as.expression(lapply(rownames(scaledExprs), function(a) bquote(italic(.(a))))),
  border_color = F,
  treeheight_col = 0,
  angle_col =315
)


plot(hs_female@reductions$pca@stdev)
PCA=20

set.seed(10)
hs_female <- FindNeighbors(hs_female, reduction = "pca", dims = 1:PCA)
hs_female <- FindClusters(hs_female, resolution = 0.8)
hs_female <- RunUMAP(hs_female, reduction = "pca", dims = 1:PCA)
hs_female <- RunTSNE(hs_female, reduction = "pca", dims = 1:PCA, check_duplicates = FALSE)

Seurat_GroupByPlot(hs_female, "right", 5, T, "condition", "umap")


markers <- list("Female-S1" = c("POU5F1", "NANOG", "SALL4", "TFAP2C"),
                "Female-S2" = c("STRA8", "ZGLP1"),
                "Female-S3" = c("SYCP1", "SYCP3", "TEX12", "SPO11"),
                "Female-S4" = c("DPPA3", "ZP3", "NOBOX", "WEE2")
)
markers_vec <- unlist(markers)
markers_vec <- markers_vec[markers_vec %in% rownames(hs_female@assays$RNA@data)]
length(markers_vec)

hs_female@active.ident <- hs_female$seurat_clusters
avgExprs <- AverageExpression(hs_female, features = markers_vec, assays = "RNA", slot = "data")
df <- data.frame(Markers = c(rep(names(markers)[1], length(markers[[1]])),
                             rep(names(markers)[2], length(markers[[2]])),
                             rep(names(markers)[3], length(markers[[3]])),
                             rep(names(markers)[4], length(markers[[4]]))
                             ))
rownames(df) <- unlist(markers)

cols <- gg_color_hue(4)
cour<-list(
  Markers=cols
)

scaledExprs <- t(scale(t(avgExprs$RNA)))
scaledExprs[scaledExprs > -min(scaledExprs)] <- -min(scaledExprs)

library(RColorBrewer)
palette_length = 100
my_color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(palette_length)

my_breaks <- c(seq(min(scaledExprs), 0,
                   length.out=ceiling(palette_length/2) + 1),
               seq(max(scaledExprs)/palette_length,
                   max(scaledExprs),
                   length.out=floor(palette_length/2)))

p <- pheatmap::pheatmap(
  scaledExprs,
  cluster_cols = T, cluster_rows = F,
  annotation_row = df, annotation_colors = list(cour$Markers),
  breaks = my_breaks, color=my_color,
  labels_row = as.expression(lapply(rownames(scaledExprs), function(a) bquote(italic(.(a))))),
  border_color = F,
  treeheight_col = 0,
  angle_col =315
)


#C7 ==> RBFOX3+/Axon+ ==> Neuron
hs_female <- hs_female[, hs_female$seurat_clusters != 7]
hs_female <- hs_female[rowSums(hs_female@assays$RNA@counts) != 0, ]

sce <- as.SingleCellExperiment(hs_female)
clusters <- quickCluster(sce, method = "igraph")
sce <- computeSumFactors(sce, clusters = clusters)
sce <- logNormCounts(sce)

dec <- modelGeneVar(sce)
top.hvgs <- getTopHVGs(dec, n = 1000)

hs_female <- as.Seurat(sce)
hs_female@reductions$PCA <- NULL
hs_female@reductions$UMAP <- NULL
hs_female@reductions$TSNE <- NULL
VariableFeatures(hs_female) <- top.hvgs
hs_female <- ScaleData(hs_female, features = top.hvgs)
hs_female <- RunPCA(hs_female, npcs = 50, weight.by.var = F)
plot(hs_female@reductions$pca@stdev)
nPC = 20

set.seed(10)
hs_female <- FindNeighbors(hs_female, reduction = "pca", dims = 1:nPC)
hs_female <- FindClusters(hs_female, resolution = 0.8)
hs_female <- RunUMAP(hs_female, reduction = "pca", dims = 1:nPC)

markers <- list("Female-S1" = c("POU5F1", "NANOG", "SALL4", "TFAP2C"),
                "Female-S2" = c("STRA8", "ZGLP1"),
                "Female-S3" = c("SYCP1", "SYCP3", "TEX12", "SPO11"),
                "Female-S4" = c("DPPA3", "ZP3", "NOBOX", "WEE2")
)
markers_vec <- unlist(markers)
markers_vec <- markers_vec[markers_vec %in% rownames(hs_female@assays$RNA@data)]
length(markers_vec)

avgExprs <- AverageExpression(hs_female, features = markers_vec, assays = "RNA", slot = "data")
df <- data.frame(Markers = c(rep(names(markers)[1], length(markers[[1]])),
                             rep(names(markers)[2], length(markers[[2]])),
                             rep(names(markers)[3], length(markers[[3]])),
                             rep(names(markers)[4], length(markers[[4]]))
))
rownames(df) <- unlist(markers)

cols <- gg_color_hue(4)
cour<-list(
  Markers=cols
)

scaledExprs <- t(scale(t(avgExprs$RNA)))
scaledExprs[scaledExprs > -min(scaledExprs)] <- -min(scaledExprs)

library(RColorBrewer)
palette_length = 100
my_color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(palette_length)

my_breaks <- c(seq(min(scaledExprs), 0,
                   length.out=ceiling(palette_length/2) + 1),
               seq(max(scaledExprs)/palette_length,
                   max(scaledExprs),
                   length.out=floor(palette_length/2)))

pheatmap::pheatmap(
  scaledExprs,
  cluster_cols = T, cluster_rows = F,
  annotation_row = df, annotation_colors = list(cour$Markers),
  breaks = my_breaks, color=my_color,
  labels_row = as.expression(lapply(rownames(scaledExprs), function(a) bquote(italic(.(a))))),
  angle_col =315
)


query_list <- list("M_FGC1" = c(6), "M_FGC2" = c(1, 2, 4, 7), "M_FGC3" = c(0, 3, 5))
hs_male <- Cluster2CellType(hs_male, QueryList = query_list)


query_list <- list("F_FGC1" = c(0, 1, 3, 5, 7), "F_FGC2" = c(4), "F_FGC3" = c(6), "F_FGC4" = c(2))
hs_female <- Cluster2CellType(hs_female, QueryList = query_list)

# saveRDS(hs_male, "SeuratObj_HumanFGC_Male.rds")
# saveRDS(hs_female, "SeuratObj_HumanFGC_Female.rds")

