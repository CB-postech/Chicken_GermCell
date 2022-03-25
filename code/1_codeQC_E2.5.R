library(SingleCellExperiment)
library(scater)
library(DropletUtils)
library(scran)
library(Seurat)

raw_path <- "F:/Device_D/SNU_Chicken/data/Raw/B7/"
data_path <- "F:/Device_D/SNU_Chicken/data/"
setwd(data_path)

sce <- read10xCounts(raw_path)

br.out <- barcodeRanks(counts(sce))


plot(br.out$rank, br.out$total, log="xy", xlab="Rank", ylab="Total")
o <- order(br.out$rank)
lines(br.out$rank[o], br.out$fitted[o], col="red")

set.seed(100)
e.out <- emptyDrops(counts(sce))  ## Cells that have UMI counts lower than 100 are empty cells.
is.cell <- e.out$FDR <= 0.05
ncell <- sum(is.cell, na.rm=TRUE)

abline(h=br.out[br.out$rank == sum(is.cell, na.rm=TRUE),]$total, col="red", lty=2)
abline(h=metadata(br.out)$knee, col="dodgerblue", lty=2)
abline(h=metadata(br.out)$inflection, col="forestgreen", lty=2)
legend("bottomleft", lty=2, col=c("dodgerblue", "forestgreen", "red"), 
       legend=c("knee", "inflection", "FDR 0.05"))


sce <- sce[,which(e.out$FDR <= 0.05)]

library(biomaRt)
date = Sys.Date()
ensembl <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset="ggallus_gene_ensembl", host="www.ensembl.org")
ensemblGenes <- getBM(attributes=c('ensembl_gene_id', 'external_gene_name',  'chromosome_name', 'gene_biotype'), mart=ensembl)
rownames(ensemblGenes) <- ensemblGenes[,1]

mtGenes = ensemblGenes[ensemblGenes[,3]=="MT",]

is.mito = rownames(sce) %in% mtGenes[,1]
length(is.mito[is.mito == TRUE])

per.cell = perCellQCMetrics(sce, subsets=list(Mito=is.mito))
colData(sce) = cbind(colData(sce), per.cell)
sce$log10_total_counts <- log10(sce$sum + 1)
sce$log10_total_features_by_counts <- log10(sce$detected + 1)
sce$pct_counts_Mito <- sce$subsets_Mito_percent
sce$subsets_Mito_percent <- NULL

h <- hist(sce$log10_total_counts, 
          breaks=100, col="grey80",
          xlab="Log-total UMI count")
cuts <- cut(h$breaks, c(0,4.0,Inf))
plot(h, col=c("red","white")[cuts])

h <- hist(sce$log10_total_features_by_counts, 
          breaks=100, col="grey80",
          xlab="Log-total number of expressed features")
cuts <- cut(h$breaks, c(-0.5,3.5,Inf))
plot(h, col=c("red","white")[cuts])


h <- hist(sce$pct_counts_Mito, 
          breaks=100, col="grey80",
          xlab="Proportion of reads in mitochondrial genes")
cuts <- cut(h$breaks, c(10,100,Inf))
plot(h, col=c("red","white")[cuts])


# Dimensionality reduction plots
sce <- runColDataPCA(sce, variables = c("sum", "detected", "subsets_Mito_sum", "subsets_Mito_detected",
                                        "log10_total_counts", "log10_total_features_by_counts", "pct_counts_Mito"))

# saveRDS(sce, "sce_E2.5_QCmetrics.rds")

# reducedDimNames(sce)
plotReducedDim(sce, dimred = "PCA_coldata")

library(RColorBrewer)
ggplot(as.data.frame(sce@int_colData$reducedDims@listData$PCA_coldata), 
       aes(x=PC1, y=PC2, color = sce$pct_counts_Mito)) +
  scale_colour_gradientn(colours = rev(brewer.pal(11, "RdBu"))) + 
  geom_point(size = 0.1) + 
  theme_bw() + 
  theme(text = element_text(size = 20),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(size = 1),
        axis.ticks = element_line(size = 1),
        legend.title = element_blank(),
        legend.key = element_blank())

##
ggplot(as.data.frame(sce@int_colData$reducedDims@listData$PCA_coldata), 
       aes(x=PC1, y=PC2, color = sce$log10_total_counts)) +
  scale_colour_gradientn(colours = rev(brewer.pal(11, "RdBu"))) + 
  geom_point(size = 0.1) + 
  theme_bw() + 
  theme(text = element_text(size = 20),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(size = 1),
        axis.ticks = element_line(size = 1),
        legend.title = element_blank(),
        legend.key = element_blank())

filtering <- numeric()
filtering <- rep(0, ncol(sce))
filtering[which(sce$log10_total_counts >= 3.5 & sce$pct_counts_Mito <= 10)] <- 1
sce$filtering <- as.integer(filtering)
print(table(sce$filtering))

ggplot(as.data.frame(sce@int_colData$reducedDims@listData$PCA_coldata), 
       aes(x=PC1, y=PC2, color = as.factor(sce$filtering))) +
  # scale_colour_gradientn(colours = rev(brewer.pal(11, "RdBu"))) + 
  geom_point(size = 1) + 
  theme_bw() + 
  theme(text = element_text(size = 20),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(size = 1),
        axis.ticks = element_line(size = 1),
        legend.title = element_blank(),
        legend.key = element_blank())

filter_by_Mt  = sce$pct_counts_Mito <= 10
filter_by_total_counts = sce$log10_total_counts >= 3.5

sce$use <- (
  filter_by_Mt &
    filter_by_total_counts
)

ggplot(as.data.frame(sce@int_colData$reducedDims@listData$PCA_coldata), 
       aes(x=PC1, y=PC2, color = as.factor(sce$use))) +
  # scale_colour_gradientn(colours = rev(brewer.pal(11, "RdBu"))) + 
  geom_point(size = 1) + 
  theme_bw() + 
  theme(text = element_text(size = 20),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(size = 1),
        axis.ticks = element_line(size = 1),
        legend.title = element_blank(),
        legend.key = element_blank())

sce_filtered <- sce[,sce$use]

# saveRDS(sce_filtered, "sce_filtered_E2.5.rds")

rownames(sce_filtered@assays@data$counts) <- rownames(sce_filtered)
colnames(sce_filtered) <- sce_filtered$Barcode
colnames(sce_filtered@assays@data$counts) <- colnames(sce_filtered)

sce_filtered@assays@data$logcounts <- log2(sce_filtered@assays@data$counts + 1)


gene_symbol <- read.table(paste0("F:/Device_D/SNU_Chicken/data/raw/B7/genes.tsv"), sep="\t")[,c(1,2)]
colnames(gene_symbol) <- c("ensembl", "symbol")
rownames(gene_symbol) <- gene_symbol$ensembl
gene_symbol$symbol <- uniquifyFeatureNames(as.vector(gene_symbol$ensembl), as.vector(gene_symbol$symbol))


rownames(sce_filtered) <- as.vector(gene_symbol[rownames(sce_filtered),2])
colnames(sce_filtered) <- paste0("B7-", substr(sce_filtered$Barcode, 1, 16))
colnames(sce_filtered@assays@data$counts) <- colnames(sce_filtered)
colnames(sce_filtered@assays@data$logcounts) <- colnames(sce_filtered)
rownames(sce_filtered@assays@data$counts) <- rownames(sce_filtered)
rownames(sce_filtered@assays@data$logcounts) <- rownames(sce_filtered)

sce_filtered@assays@data$counts["DAZL",] <- max(sce_filtered@assays@data$counts["DAZL",], sce_filtered@assays@data$counts["DAZL-EGFP",])
sce_filtered@assays@data$logcounts["DAZL",] <- log2(sce_filtered@assays@data$counts["DAZL",] + 1)
sce_filtered <- sce_filtered[rownames(sce_filtered) != "DAZL-EGFP",]

sce_filtered$sample_index <- c(rep("B7", ncol(sce_filtered)))
sce_filtered$condition <- c(rep("E2.5", ncol(sce_filtered)))

sce_filtered <- sce_filtered[rowSums(sce_filtered@assays@data$counts) != 0,]

#Normalization
clusters <- quickCluster(sce_filtered, method="igraph")
sce_filtered <- computeSumFactors(sce_filtered, clusters=clusters)
sce_filtered <- logNormCounts(sce_filtered)

# saveRDS(sce_filtered, "sce_filtered_E2.5_normalized.rds")

#HVG selection
dec <- modelGeneVar(sce_filtered)
plot(dec$mean, dec$total, xlab="Mean log-expression", ylab="Variance")
curve(metadata(dec)$trend(x), col="blue", add=TRUE)
top.hvgs <- getTopHVGs(dec, n = 1000)
length(top.hvgs)

seurat <- as.Seurat(sce_filtered)
seurat@assays$RNA@var.features <- top.hvgs

seurat <- ScaleData(seurat, features = top.hvgs)

#Dim reduction
seurat <- RunPCA(seurat, npcs = 50, weight.by.var = F)
plot(seurat@reductions$pca@stdev)
PCA = 15

set.seed(10)
seurat <- FindNeighbors(seurat, reduction = "pca", dims = 1:PCA)
seurat <- FindClusters(seurat, resolution = 1)
seurat <- RunUMAP(seurat, reduction = "pca", dims = 1:PCA)

DimPlot(seurat, group.by = "seurat_clusters", pt.size = 1)
FeaturePlot(seurat, features = "DAZL", pt.size = 1, order = T) +
  theme(axis.title = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        legend.position = "right")

#calculate Z/W chromosome scores for assigning sexes
Wchr_gene <- rownames(ensemblGenes[ensemblGenes$chromosome_name == "W", ])
Wchr_gene <- gene_symbol[Wchr_gene, "symbol"]
Wchr_gene <- rownames(seurat)[rownames(seurat) %in% Wchr_gene]

seurat <- ScaleData(seurat, features = Wchr_gene)
seurat$Wchr_score <- colMeans(seurat@assays$RNA@scale.data)
FeaturePlot(seurat, features = "Wchr_score",
            cols = c("#e0e0e0", "#b2182b"), pt.size = 1, order = T,
            min.cutoff = -min(c(abs(min(seurat$Wchr_score)), abs(max(seurat$Wchr_score)))),
            max.cutoff = min(c(abs(min(seurat$Wchr_score)), abs(max(seurat$Wchr_score))))) +
  scale_colour_gradientn(colours = rev(brewer.pal(11, "RdBu"))) +
  theme(axis.title = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        legend.position = "right")

Zchr_gene <- rownames(ensemblGenes[ensemblGenes$chromosome_name == "Z", ])
Zchr_gene <- gene_symbol[Zchr_gene, "symbol"]
Zchr_gene <- rownames(seurat)[rownames(seurat) %in% Zchr_gene]

seurat <- ScaleData(seurat, features = Zchr_gene)
seurat$Zchr_score <- colMeans(seurat@assays$RNA@scale.data)
FeaturePlot(seurat, features = "Zchr_score",
            cols = c("#e0e0e0", "#b2182b"), pt.size = 1, order = T,
            min.cutoff = -min(c(abs(min(seurat$Zchr_score)), abs(max(seurat$Zchr_score)))),
            max.cutoff = min(c(abs(min(seurat$Zchr_score)), abs(max(seurat$Zchr_score))))) +
  scale_colour_gradientn(colours = rev(brewer.pal(11, "RdBu"))) +
  theme(axis.title = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        legend.position = "right")


DimPlot(seurat, cells.highlight = names(seurat$Wchr_score[seurat$Wchr_score == min(seurat$Wchr_score)]))
DimPlot(seurat, cells.highlight = colnames(seurat[, seurat@assays$RNA@counts["HINTW",] == 0]))

#assign sexes
sex_vec <- c(rep(NA, ncol(seurat)))
sex_vec[seurat@active.ident == 1 | seurat@active.ident == 4] <- "Male"
sex_vec[seurat@active.ident == 0 | seurat@active.ident == 2 | seurat@active.ident == 3] <- "Female"
names(sex_vec) <- colnames(seurat)

seurat$Sex <- sex_vec

#Remove erythrocytes
out_rm <- list()
seurat_rm <- seurat[, seurat@active.ident != 5]

out_rm[["Male"]] <- seurat_rm[, seurat_rm$Sex == "Male"]
out_rm[["Female"]] <- seurat_rm[, seurat_rm$Sex == "Female"]

