library(data.table)
library(Seurat)
library(SingleCellExperiment)
library(scater)
library(scran)

setwd("raw/normalized/Human/GSE86146_RAW/")

count_files <- list.files()
female_count <- list()
male_count <- list()

grep("_F_", count_files, value = T)[-length(grep("_F_", count_files))]

for(i in grep("_F_", count_files, value = T)[-length(grep("_F_", count_files))]){
  female_count[[i]] <- as.data.frame(fread(i, sep = "\t"))
  rownames(female_count[[i]]) <- female_count[[i]]$Gene
  female_count[[i]] <- female_count[[i]][, -1]
}

grep("_M_", count_files, value = T)[-1]


for(i in grep("_M_", count_files, value = T)[-1]){
  male_count[[i]] <- as.data.frame(fread(i, sep = "\t"))
  rownames(male_count[[i]]) <- male_count[[i]]$Gene
  male_count[[i]] <- male_count[[i]][, -1]
}



m4_f11 <- as.data.frame(fread("GSM2306023_M_4W_embryo1_and_F_11W_embryo1_gene_expression.txt",
                              sep = "\t"))
rownames(m4_f11) <- m4_f11$Gene
m4_f11 <- m4_f11[, -1]


length(grep("M_4W", colnames(m4_f11)))


print(dim(m4_f11))
for(i in names(female_count)){
  print(dim(female_count[[i]]))
}
for(i in names(male_count)){
  print(dim(male_count[[i]]))
}


male_count[["GSM2306023_M_4W"]] <- m4_f11[, grep("M_4W", colnames(m4_f11), value = T)]
female_count[["GSM2306023_F_11W"]] <- m4_f11[, grep("F_11W", colnames(m4_f11), value = T)]

#Check the order of gene name
for(i in 1:(length(names(female_count))-1)){
  try(print(identical(rownames(female_count[[i]]), rownames(female_count[[i+1]]))))
}
for(i in 1:(length(names(male_count))-1)){
  try(print(identical(rownames(male_count[[i]]), rownames(male_count[[i+1]]))))
}

#check the number of cells
n = 0
for(i in 1:(length(names(female_count)))){
  n = n + ncol(female_count[[i]])
}
print(paste0("number of female cell: ", n))

n = 0
for(i in 1:(length(names(male_count)))){
  n = n + ncol(male_count[[i]])
}
print(paste0("number of male cell: ", n))

#Merge Expression Matrix-Female
female_expMat <- data.frame(row.names = rownames(female_count[[1]]))
for(i in 1:(length(names(female_count)))){
  female_expMat <- cbind(female_expMat, female_count[[i]])
}
print(dim(female_expMat))

#Merge Expression Matrix-Male
male_expMat <- data.frame(row.names = rownames(male_count[[1]]))
for(i in 1:(length(names(male_count)))){
  male_expMat <- cbind(male_expMat, male_count[[i]])
}
print(dim(male_expMat))

#toSeuratObj
################################Female############################
seurat_female <- CreateSeuratObject(counts = female_expMat, project = "HumanFGC_Female",
                                    min.cells = 0, min.features = 0)
sce <- as.SingleCellExperiment(seurat_female)
sce@assays@data$logcounts <- log2(sce@assays@data$counts + 1)
sce@assays@data$logcounts <- as(sce@assays@data$logcounts, "dgCMatrix")
dec <- modelGeneVar(sce)
plot(dec$mean, dec$total, xlab="Mean log-expression", ylab="Variance")
curve(metadata(dec)$trend(x), col="blue", add=TRUE)

top.hvgs <- getTopHVGs(dec, fdr.threshold = 0.05)
length(top.hvgs)

seurat_female <- as.Seurat(sce)

seurat_female@assays$RNA@var.features <- top.hvgs
seurat_female <- ScaleData(seurat_female, features = top.hvgs)

seurat_female <- RunPCA(seurat_female, pcs.compute = 50, weight.by.var = F)

plot(seurat_female@reductions$pca@stdev)
PCA = 10

set.seed(10)
seurat_female <- FindNeighbors(seurat_female, reduction = "pca", dim = 1:PCA)
seurat_female <- FindClusters(seurat_female, resolution = 0.8)
seurat_female <- RunUMAP(seurat_female, reduction = "pca", dims = 1:PCA)


#Add batch infos into seuratObj
batch_vec <- c()
for(i in colnames(seurat_female)){
  tmp <- paste0(sapply(strsplit(i, "_"), `[`, 1), "_", sapply(strsplit(i, "_"), `[`, 2), "_",
                sapply(strsplit(i, "_"), `[`, 3))
  batch_vec <- c(batch_vec, tmp)
}

names(batch_vec) <- colnames(seurat_female)
seurat_female$batch <- factor(batch_vec, levels = c('F_5W_embryo1', 'F_5W_embryo2', 'F_7W_embryo1',
                                                    'F_8W_embryo1', 'F_10W_embryo1', 'F_11W_embryo1',
                                                    'F_12W_embryo1', 'F_14W_embryo1', 'F_18W_embryo1',
                                                    'F_18W_embryo2', 'F_20W_embryo1', 'F_20W_embryo2',
                                                    'F_23W_embryo1', 'F_23W_embryo2', 'F_24W_embryo1',
                                                    'F_24W_embryo2', 'F_26W_embryo1'
))

#Add condition infos into seuratObj
cond_vec <- c()
for(i in colnames(seurat_female)){
  tmp <- paste0(sapply(strsplit(i, "_"), `[`, 1), "_", sapply(strsplit(i, "_"), `[`, 2))
  cond_vec <- c(cond_vec, tmp)
}
names(cond_vec) <- colnames(seurat_female)
seurat_female$condition <- factor(cond_vec, levels = c('F_5W', 'F_7W',
                                                       'F_8W', 'F_10W', 'F_11W',
                                                       'F_12W', 'F_14W', 'F_18W',
                                                       'F_20W',
                                                       'F_23W', 'F_24W', 'F_26W'
))


DimPlot(seurat_female) +
  theme(axis.ticks = element_blank(),
        axis.line = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title = element_blank())


for(i in c("DAZL", "DDX4", "DND1", "NANOG")){
  p <- FeaturePlot(seurat_female, i, pt.size = 1, order = T) +
          theme(axis.ticks = element_blank(),
                axis.line = element_blank(),
                axis.text.x = element_blank(),
                axis.text.y = element_blank(),
                axis.title = element_blank()) +
          scale_colour_gradientn(colours = rev(c("#300000", "red", "#eeeeee")))
  print(p)
}


#subset fgc
fgc_female <- seurat_female[, seurat_female$seurat_clusters %in% c(0, 3, 5, 7, 9)]
fgc_female <- fgc_female[rowSums(fgc_female@assays$RNA@counts) != 0,]

sce <- as.SingleCellExperiment(fgc_female)

dec <- modelGeneVar(sce)
plot(dec$mean, dec$total, xlab="Mean log-expression", ylab="Variance")
curve(metadata(dec)$trend(x), col="blue", add=TRUE)

top.hvgs <- getTopHVGs(dec, n = 1000)
length(top.hvgs)

fgc_female <- as.Seurat(sce)
fgc_female@assays$RNA@var.features <- top.hvgs
fgc_female <- ScaleData(fgc_female, features = top.hvgs)
fgc_female <- RunPCA(fgc_female, pcs.compute = 50, weight.by.var = F)
plot(fgc_female@reductions$pca@stdev)
PCA = 10

fgc_female@reductions$PCA <- NULL
fgc_female@reductions$UMAP <- NULL

set.seed(10)
fgc_female <- FindNeighbors(fgc_female, reduction = "pca", dim = 1:PCA)
fgc_female <- FindClusters(fgc_female, resolution = 0.8)
fgc_female <- RunUMAP(fgc_female, reduction = "pca", dims = 1:PCA)
DimPlot(fgc_female) +
  theme(axis.ticks = element_blank(),
        axis.line = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title = element_blank())



#########################Male#####################################
#toSeuratObj
seurat_male <- CreateSeuratObject(counts = male_expMat, project = "HumanFGC_Male",
                                  min.cells = 0, min.features = 0)
sce <- as.SingleCellExperiment(seurat_male)
sce@assays@data$logcounts <- log2(sce@assays@data$counts + 1)
sce@assays@data$logcounts <- as(sce@assays@data$logcounts, "dgCMatrix")

dec <- modelGeneVar(sce)
plot(dec$mean, dec$total, xlab="Mean log-expression", ylab="Variance")
curve(metadata(dec)$trend(x), col="blue", add=TRUE)

top.hvgs <- getTopHVGs(dec, n = 1000)
length(top.hvgs)

seurat_male <- as.Seurat(sce)

seurat_male@assays$RNA@var.features <- top.hvgs
seurat_male <- ScaleData(seurat_male, features = top.hvgs)

seurat_male <- RunPCA(seurat_male, pcs.compute = 50, weight.by.var = F)
plot(seurat_male@reductions$pca@stdev)
PCA = 10

set.seed(10)
seurat_male <- FindNeighbors(seurat_male, reduction = "pca", dim = 1:PCA)
seurat_male <- FindClusters(seurat_male, resolution = 0.8)
seurat_male <- RunUMAP(seurat_male, reduction = "pca", dims = 1:PCA)

#Add batch infos into seuratObj
batch_vec <- c()
for(i in colnames(seurat_male)){
  tmp <- paste0(sapply(strsplit(i, "_"), `[`, 1), "_", sapply(strsplit(i, "_"), `[`, 2), "_",
                sapply(strsplit(i, "_"), `[`, 3))
  batch_vec <- c(batch_vec, tmp)
}

names(batch_vec) <- colnames(seurat_male)
seurat_male$batch <- factor(batch_vec, levels = c('M_4W_embryo1', 'M_9W_embryo1', 'M_10W_embryo1', 'M_10W_embryo2',
                                                  'M_12W_embryo1', 'M_19W_embryo1', 'M_19W_embryo2', 'M_20W_embryo1',
                                                  'M_21W_embryo1', 'M_21W_embryo2', 'M_21W_embryo3', 'M_25W_embryo1'
))

DimPlot(seurat_male, group.by = "batch") +
  theme(axis.ticks = element_blank(),
        axis.line = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title = element_blank())

#Add condition infos into seuratObj
cond_vec <- c()
for(i in colnames(seurat_male)){
  tmp <- paste0(sapply(strsplit(i, "_"), `[`, 1), "_", sapply(strsplit(i, "_"), `[`, 2))
  cond_vec <- c(cond_vec, tmp)
}
names(cond_vec) <- colnames(seurat_male)
seurat_male$condition <- factor(cond_vec, levels = c('M_4W', 'M_9W', 'M_10W',
                                                     'M_12W', 'M_19W', 'M_20W',
                                                     'M_21W', 'M_25W'
))



DimPlot(seurat_male) +
  theme(axis.ticks = element_blank(),
        axis.line = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title = element_blank())

for(i in c("DAZL", "DDX4", "DND1", "NANOG")){
  p <- FeaturePlot(seurat_male, i, pt.size = 1, order = T) +
          theme(axis.ticks = element_blank(),
                axis.line = element_blank(),
                axis.text.x = element_blank(),
                axis.text.y = element_blank(),
                axis.title = element_blank()) +
          scale_colour_gradientn(colours = rev(c("#300000", "red", "#eeeeee")))
  print(p)
}


#subset fgc
fgc_male <- seurat_male[, seurat_male$seurat_clusters %in% c(0, 2, 3, 5, 8, 9, 10)]
fgc_male <- fgc_male[rowSums(fgc_male@assays$RNA@counts) != 0,]

sce <- as.SingleCellExperiment(fgc_male)

dec <- modelGeneVar(sce)
plot(dec$mean, dec$total, xlab="Mean log-expression", ylab="Variance")
curve(metadata(dec)$trend(x), col="blue", add=TRUE)

top.hvgs <- getTopHVGs(dec, n = 1000)


fgc_male <- as.Seurat(sce)
fgc_male@assays$RNA@var.features <- top.hvgs
fgc_male <- ScaleData(fgc_male, features = top.hvgs)
fgc_male <- RunPCA(fgc_male, pcs.compute = 50, weight.by.var = F)
plot(fgc_male@reductions$pca@stdev)
PCA = 7

fgc_male@reductions$PCA <- NULL
fgc_male@reductions$UMAP <- NULL

set.seed(10)
fgc_male <- FindNeighbors(fgc_male, reduction = "pca", dim = 1:PCA)
fgc_male <- FindClusters(fgc_male, resolution = 0.8)
fgc_male <- RunUMAP(fgc_male, reduction = "pca", dims = 1:PCA)

DimPlot(fgc_male) +
  theme(axis.ticks = element_blank(),
        axis.line = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title = element_blank())


saveRDS(fgc_female, "SeuratObj_HumanFGC_Female_logNorm.rds")
saveRDS(fgc_male, "SeuratObj_HumanFGC_Male_logNorm.rds")
