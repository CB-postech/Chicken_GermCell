library(SingleCellExperiment)
library(scater)
library(scran)
library(Seurat)

seurat <- readRDS("seurat_merged_QC_CCcalculated.rds")


##############Extract Male Samples##################
seurat <- seurat[, seurat$sex == "Male"]

sce_male <- as.SingleCellExperiment(seurat)
sce_male <- sce_male[rowSums(sce_male@assays@data$counts) != 0,]

#Normalization
clusters <- quickCluster(sce_male, method="igraph")
sce_male <- computeSumFactors(sce_male, clusters=clusters)
sce_male <- logNormCounts(sce_male)

# saveRDS(sce_male, "sce_male_merged_QC_normalized.rds")

#HVG selection
dec <- modelGeneVar(sce_male)
plot(dec$mean, dec$total, xlab="Mean log-expression", ylab="Variance")
curve(metadata(dec)$trend(x), col="blue", add=TRUE)
top.hvgs <- getTopHVGs(dec, fdr.threshold = 0.05) # 434 genes
length(top.hvgs)

seurat <- as.Seurat(sce_male)
seurat@assays$RNA@var.features <- top.hvgs

seurat@reductions$PCA <- NULL
seurat@reductions$UMAP <- NULL

seurat <- ScaleData(seurat, features = top.hvgs)
#dimensionality reduction
seurat <- RunPCA(seurat, npcs = 50, weight.by.var = F)
plot(seurat@reductions$pca@stdev)
PCA = 10

set.seed(10)
seurat <- FindNeighbors(seurat, reduction = "pca", dims = 1:PCA)
seurat <- FindClusters(seurat, resolution = 0.6)
seurat <- RunUMAP(seurat, reduction = "pca", dims = 1:PCA)

#CCcorrection
#ortholog definition for cell cycle phase genes
ortho <- read.table("F:/Device_D/SNU_Chicken/data/ortholog/human_chicken_hcop_fifteen_column.txt", header = T, sep = "\t") #ensembl orthology
load("F:/Device_D/SNU_Chicken/data/ensemblGenes 20200227 .RData") #ensembl biomart

hist(nchar(grep("^LOC", ortho$chicken_symbol, value = T)))
ortho$chicken_symbol[grep("^LOC", ortho$chicken_symbol)] <- ortho$chicken_ensembl_gene[grep("^LOC", ortho$chicken_symbol)]

ortho_hum4chic <- ortho[ortho$chicken_symbol %in% rownames(seurat), ] #Ortholog of human genes against chicken, for detected features
ortho_hum4chic <- ortho_hum4chic[!duplicated(ortho_hum4chic$human_symbol) & !duplicated(ortho_hum4chic$chicken_symbol),]
rownames(ortho_hum4chic) <- ortho_hum4chic$human_symbol

#CellCycleScoring
s.genes <- ortho_hum4chic[cc.genes.updated.2019$s.genes[cc.genes.updated.2019$s.genes %in% ortho_hum4chic$human_symbol], ]$chicken_symbol
g2m.genes <- ortho_hum4chic[cc.genes.updated.2019$g2m.genes[cc.genes.updated.2019$g2m.genes %in% ortho_hum4chic$human_symbol], ]$chicken_symbol

seurat <- CellCycleScoring(seurat, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

seurat@active.ident <- factor(as.vector(seurat@active.ident), levels = c("G1", "G2M", "S"))
seurat$CellCycle <- factor(seurat@active.ident, levels = c("G1", "G2M", "S"))
seurat@active.ident <- seurat$seurat_clusters

#CellCycle Correction
seurat$CC.Difference <- seurat$S.Score - seurat$G2M.Score
seurat <- ScaleData(seurat, vars.to.regress = "CC.Difference", features = rownames(seurat))

seurat@reductions$PCA <- NULL
seurat@reductions$UMAP <- NULL

#dimensionality reduction
seurat <- RunPCA(seurat , npcs = 50, weight.by.var = F, features = VariableFeatures(seurat))
plot(seurat@reductions$pca@stdev)
PCA = 10

set.seed(10)
seurat <- FindNeighbors(seurat, reduction = "pca", dims = 1:PCA)
seurat <- FindClusters(seurat, resolution = 0.6)
seurat <- RunUMAP(seurat, reduction = "pca", dims = 1:PCA)

seurat <- CellCycleScoring(seurat, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
seurat@active.ident <- factor(as.vector(seurat@active.ident), levels = c("G1", "G2M", "S"))
seurat$CellCycle <- seurat@active.ident
seurat@active.ident <- seurat$seurat_clusters

# saveRDS(seurat, "seurat_male_QC_CCcorrected.rds")

