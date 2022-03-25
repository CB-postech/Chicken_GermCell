library(SingleCellExperiment)
library(scater)
library(scran)
library(Seurat)

gg_male <- readRDS("seurat_male_QC_CCcorrected.rds")
hs_male <-  readRDS("SeuratObj_HumanFGC_Male.rds")

# Extract orthologous gene from scRNAseq data
ortho <- read.csv("../Multi-species/data/human_chicken_ortholog.txt", sep = "\t")
for(i in c(1:4)){ortho[[i]] <- as.vector(ortho[[i]])}
head(ortho)

ortho[ortho$Chicken.gene.name == "", "Chicken.gene.name"] <- ortho[ortho$Chicken.gene.name == "", "Chicken.gene.stable.ID"]

ortho <- ortho[ortho$Chicken.gene.name %in% rownames(gg_male),]
ortho <- ortho[ortho$Gene.name %in% rownames(hs_male),]
dim(ortho)

ortho <- ortho[!duplicated(ortho$Chicken.gene.name) & !duplicated(ortho$Gene.name),]
dim(ortho)

VariableFeatures(gg_male) <- gg_male@assays$RNA@var.features$hvg_434

gg_male@assays$RNA@counts <- gg_male@assays$RNA@counts[rownames(gg_male) %in% ortho$Chicken.gene.name,]
gg_male@assays$RNA@data <- gg_male@assays$RNA@data[rownames(gg_male) %in% ortho$Chicken.gene.name,]
gg_male@assays$RNA@scale.data <- gg_male@assays$RNA@scale.data[rownames(gg_male@assays$RNA@scale.data) %in% ortho$Chicken.gene.name,]
gg_male@assays$RNA@var.features <- gg_male@assays$RNA@var.features[gg_male@assays$RNA@var.features %in% ortho$Chicken.gene.name]


rownames(ortho) <- ortho$Chicken.gene.name

rownames(gg_male@assays$RNA@counts) <- ortho[rownames(gg_male),  "Gene.name"]
rownames(gg_male@assays$RNA@data) <- ortho[rownames(gg_male),  "Gene.name"]
rownames(gg_male@assays$RNA@scale.data) <- ortho[rownames(gg_male@assays$RNA@scale.data),  "Gene.name"]
gg_male@assays$RNA@var.features <- ortho[gg_male@assays$RNA@var.features,  "Gene.name"]

hs_male <- hs_male[ortho$Gene.name,]

#mapping using coordinates of k-NN
library(SingleCellExperiment)
library(KernelKnn)

sce <- SingleCellExperiment(list(counts = gg_male@assays$RNA@counts))
logcounts(sce) <- gg_male@assays$RNA@data

#Select Highly variable genes
hvg <- VariableFeatures(gg_male)
hs_male <- ScaleData(hs_male, features = hvg)

hvg_mat <- as.matrix(gg_male@assays$RNA@scale.data)[hvg,]
hvg_test_mat <- as.matrix(hs_male@assays$RNA@scale.data[hvg,])

indexN <- KernelKnn::knn.index.dist(t(hvg_mat),t(hvg_test_mat),k=5,threads = 4, method = "pearson_correlation")
iN2 <- indexN$test_knn_idx
rownames(iN2) <- colnames(hvg_test_mat)


iN3 <- apply(iN2, 2, function(x) colnames(hvg_mat)[x])
idxx <- apply(iN3, 2, function(x) gg_male@reductions$umap@cell.embeddings[x,1])
idxy <- apply(iN3, 2, function(x) gg_male@reductions$umap@cell.embeddings[x,2])

prjx <- rowMeans(idxx)
prjy <- rowMeans(idxy)
names(prjx) <- rownames(iN2)
names(prjy) <- rownames(iN2)


#Plot
require(scales)

# Create vector with levels of object@ident
identities <- levels(gg_male$CellType)

# Create vector of default ggplot2 colors
my_color_palette <- hue_pal()(length(identities))       

df_test=data.frame(x=rbind(as.matrix(gg_male@reductions$umap@cell.embeddings[,1]),as.matrix(prjx)), 
                   y=rbind(as.matrix(gg_male@reductions$umap@cell.embeddings[,2]),as.matrix(prjy)), 
                   expression = c(rep("gg",length(colnames(gg_male))),paste0("hs",hs_male$CellType)))
Index = c(rep("gg",length(colnames(gg_male))),paste0("hs",hs_male$CellType))
Index = factor(Index, levels=c("gg", paste0("hs",levels(hs_male$CellType))))

for(i in c(1:length(levels(hs_male$CellType)))){
  col_vec <- RColorBrewer::brewer.pal(9, "YlGnBu")[c(4, 6, 8)]
  col_vec[col_vec != col_vec[i]] <- "grey79"
  
  p <- ggplot(df_test,
              aes(x=x, y=y, color = Index)) +
    geom_point(size=1) +
    scale_color_manual(values = c("grey79", col_vec))+
    guides(colour = guide_legend(override.aes = list(size=10)))+
    theme_bw() +
    theme(text = element_text(size = 20),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          axis.line = element_line(size = 1),
          axis.ticks = element_line(size = 1),
          legend.title = element_blank(),
          legend.key = element_blank())
  
  # ggsave(plot = p, width = 7.5, file = paste0("gg_hs_male_kNNaverage_", levels(gg_male$CellType)[i], ".pdf"))
}
col_vec <- RColorBrewer::brewer.pal(9, "YlGnBu")[c(4, 6, 8)]

p <- ggplot(df_test,
            aes(x=x, y=y, color = Index)) +
  geom_point(size=1) +
  scale_color_manual(values = c("grey79", col_vec))+
  guides(colour = guide_legend(override.aes = list(size=10)))+
  theme_bw() +
  theme(text = element_text(size = 20),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(size = 1),
        axis.ticks = element_line(size = 1),
        legend.title = element_blank(),
        legend.key = element_blank())

# ggsave(plot = p, width = 7.5, file = paste0("gg_hs_male_kNNaverage.pdf"))


p <- ggplot(df_test,
            aes(x=x, y=y, color = Index)) +
  geom_point(size=1) +
  scale_color_manual(values = c("grey79", RColorBrewer::brewer.pal(9, "YlGnBu")[c(4, 6, 8)]))+
  guides(colour = guide_legend(override.aes = list(size=10)))+
  theme_bw() +
  theme(text = element_text(size = 20),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(size = 1),
        axis.ticks = element_line(size = 1),
        legend.title = element_blank(),
        legend.key = element_blank())

# ggsave(plot = p, width = 7.5, file = paste0("gg_hs_male_kNNaverage_", levels(gg_male$CellType)[i], ".pdf"))



#assign nearest stage of the other species
rownames(iN3) <- colnames(hs_male)
mapped_vec <- c()
for(i in c(1:nrow(iN3))){
  tmp_gg <- as.vector(gg_male$CellType[iN3[i, ]])
  tmp_count <- table(tmp_gg)[order(table(tmp_gg), decreasing = T)]
  if(tmp_count[1] >= 3){
    mapped_vec <- c(mapped_vec, names(tmp_count)[1])
  }else{
    mapped_vec <- c(mapped_vec, NA)
  }
}

names(mapped_vec) <- rownames(iN3)
hs_male$gg_mapped <- factor(mapped_vec, levels = levels(gg_male$CellType))

df_test=data.frame(x=rbind(as.matrix(gg_male@reductions$umap@cell.embeddings[,1]),as.matrix(prjx)), 
                   y=rbind(as.matrix(gg_male@reductions$umap@cell.embeddings[,2]),as.matrix(prjy)), 
                   expression = c(rep("gg",length(colnames(gg_male))),paste0("hs",hs_male$gg_mapped)))
Index = c(rep("gg",length(colnames(gg_male))),paste0("hs",hs_male$gg_mapped))
Index = factor(Index, levels=c("gg", paste0("hs",levels(hs_male$gg_mapped))))

for(i in c(1:(length(levels(gg_male$CellType))+1))){
  if(i <= 4){
    col_vec <- c("#F297E5", "#D09EF5", "#AB87ED", "#68369A")
    col_vec[col_vec != col_vec[i]] <- "grey79"
    
    p <- ggplot(df_test,
                aes(x=x, y=y, color = Index)) +
      geom_point(size=1) +
      scale_color_manual(values = c("grey79", col_vec), na.value = "grey79")+
      guides(colour = guide_legend(override.aes = list(size=10)))+
      theme_bw() +
      theme(text = element_text(size = 20),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            axis.line = element_line(size = 1),
            axis.ticks = element_line(size = 1),
            legend.title = element_blank(),
            legend.key = element_blank())
    
    # ggsave(plot = p, width = 7.5, file = paste0("gg_hs_male_kNNaverage_", levels(gg_male$CellType)[i], "_mappingResult.pdf"))
  }else{
    col_vec <- rep("grey79", length(levels(gg_male$CellType)))
    
    p <- ggplot(df_test,
                aes(x=x, y=y, color = Index)) +
      geom_point(size=1) +
      scale_color_manual(values = c("grey79", col_vec), na.value = "#000000")+
      guides(colour = guide_legend(override.aes = list(size=10)))+
      theme_bw() +
      theme(text = element_text(size = 20),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            axis.line = element_line(size = 1),
            axis.ticks = element_line(size = 1),
            legend.title = element_blank(),
            legend.key = element_blank())
    
    # ggsave(plot = p, width = 7.5, file = paste0("gg_hs_male_kNNaverage_NAs_mappingResult.pdf"))
  }
}

p <- ggplot(df_test,
            aes(x=x, y=y, color = Index)) +
  geom_point(size=1) +
  scale_color_manual(values = c("grey79", c("#F297E5", "#D09EF5", "#AB87ED", "#68369A")), na.value = "#000000")+
  guides(colour = guide_legend(override.aes = list(size=10)))+
  theme_bw() +
  theme(text = element_text(size = 20),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(size = 1),
        axis.ticks = element_line(size = 1),
        legend.title = element_blank(),
        legend.key = element_blank())

# ggsave(plot = p, width = 7.5, file = paste0("gg_hs_male_kNNaverage_mappingResult.pdf"))


df_bar <- data.frame(hs = hs_male$CellType, gg = hs_male$gg_mapped)
table(is.na(hs_male$gg_mapped)) # NA = 1156

df_barplt <- rbind(data.frame("Hs_stage" = rep("Hs_Male_S1", 4), "Gg_mapped" = levels(gg_male$CellType), "Freq" = c(table(df_bar[df_bar$hs == "M_FGC1", "gg"]))),
                   data.frame("Hs_stage" = rep("Hs_Male_S2", 4), "Gg_mapped" = levels(gg_male$CellType), "Freq" = c(table(df_bar[df_bar$hs == "M_FGC2", "gg"])))
)
df_barplt <- rbind(df_barplt,
                   data.frame("Hs_stage" = rep("Hs_Male_S3", 4), "Gg_mapped" = levels(gg_male$CellType), "Freq" = c(table(df_bar[df_bar$hs == "M_FGC3", "gg"])))
)

rownames(df_barplt) <- c(1:nrow(df_barplt))


# Plot for frequency
p_freq <- ggplot(df_barplt, aes(fill=df_barplt$Gg_mapped, y=df_barplt$Freq, x=df_barplt$Hs_stage)) + 
  geom_text(aes(label=Freq), vjust = -0.3,
            position=position_dodge(0.9), size=4) +
  geom_bar(position="dodge", stat="identity") +
  scale_fill_manual(values = c("#F297E5", "#D09EF5", "#AB87ED", "#68369A")) +
  theme_bw() +
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank()
  )+
  labs(fill = "Gg_mapped")

#Add Proportion
prop_vec <- c()
for(i in c(1:length(levels(hs_male$CellType)))){
  subset_prop <- (((i-1)*length(unique(gg_male$CellType))+1):(i*(length(unique(gg_male$CellType))))) #(1,2,3,4,5), (6, 7, 8, 9, 10)...
  prop_vec <- c(prop_vec, df_barplt[subset_prop, "Freq"]/sum(df_barplt[subset_prop, "Freq"]))
  print(df_barplt[subset_prop,])
}
prop_vec <- paste(round(prop_vec, 3)*100, "%")


df_barplt$Prop <- prop_vec

p_prop <- ggplot(df_barplt, aes(fill=df_barplt$Gg_mapped, y=df_barplt$Freq, x=df_barplt$Hs_stage)) + 
  geom_bar(position="fill", stat="identity") +
  geom_text(aes(label=df_barplt$Prop),position=position_fill(vjust=0.5), size=3) +
  scale_fill_manual(values = c("#F297E5", "#D09EF5", "#AB87ED", "#68369A")) +
  theme_bw() +
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank()
  )+
  labs(fill = "Gg_mapped")

p <- cowplot::plot_grid(p_freq, p_prop, ncol = 1)

df <- data.frame(pt_bin = hs_male$gg_mapped, stage = hs_male$CellType)
ggplot(df) +
  aes(x = stage, fill = factor(pt_bin)) +
  geom_bar(position = "fill")