.libPaths(Sys.getenv("R_LIBS_Chicken"))

library(SingleCellExperiment)
library(scater)
library(DropletUtils)
library(scran)

raw_path <- "D:/dgcha/projects/SNU_Chicken/data/Raw/"
data_path <- "D:/dgcha/projects/SNU_Chicken/data/"

setwd(data_path)

#load scRNA-seq data/remove empty droplets
sce <- list()
br.out <- list()
e.out <- list()
is.cell <- list()
ncell <- list()

for(i in 1:13){
  sample_index <- list.files(raw_path)[i]
  sce[[sample_index]] <- read10xCounts(paste0("D:/dgcha/projects/SNU_Chicken/data/Raw/",list.files(raw_path)[i]))
  
  br.out[[sample_index]] <- barcodeRanks(counts(sce[[sample_index]]))

  pdf(paste0("DropletUtils_RankPlot_", sample_index, ".pdf"))
    plot(br.out[[sample_index]]$rank, br.out[[sample_index]]$total, log="xy", xlab="Rank", ylab="Total")
    o <- order(br.out[[sample_index]]$rank)
    lines(br.out[[sample_index]]$rank[o], br.out[[sample_index]]$fitted[o], col="red")
    
    set.seed(100)
    e.out[[sample_index]] <- emptyDrops(counts(sce[[sample_index]]))  ## Cells that have UMI counts lower than 100 are empty cells.
    is.cell[[sample_index]] <- e.out[[sample_index]]$FDR <= 0.05
    ncell[[sample_index]] <- sum(is.cell[[sample_index]], na.rm=TRUE)
    
    abline(h=br.out[[sample_index]][br.out[[sample_index]]$rank == sum(is.cell[[sample_index]], na.rm=TRUE)-1,]$total, col="red", lty=2)
    abline(h=metadata(br.out[[sample_index]])$knee, col="dodgerblue", lty=2)
    abline(h=metadata(br.out[[sample_index]])$inflection, col="forestgreen", lty=2)
    legend("bottomleft", lty=2, col=c("dodgerblue", "forestgreen", "red"), 
           legend=c("knee", "inflection", "FDR_0.05"))
  dev.off()
}

for(i in 1:13){
  sample_index <- list.files(raw_path)[i]
  sce[[sample_index]] <- sce[[sample_index]][,which(e.out[[sample_index]]$FDR <= 0.05)]
}

save.image("QC_v1.RData")


#######################################

for(i in 1:13){
  sample_index <- list.files(raw_path)[i]
  sce[[sample_index]]@assays@data$logcounts <- log2(sce[[sample_index]]@assays@data$counts + 1)
}


#to load gene symbols/mitochondrial gene info
library(biomaRt)
date = Sys.Date()
ensembl <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset="ggallus_gene_ensembl", host="www.ensembl.org")
ensemblGenes <- getBM(attributes=c('ensembl_gene_id', 'external_gene_name',  'chromosome_name', 'gene_biotype'), mart=ensembl)
rownames(ensemblGenes) <- ensemblGenes[,1]

mtGenes = ensemblGenes[ensemblGenes[,3]=="MT",]

#calculate QC factors
for(i in 1:13){
  sample_index <- list.files(raw_path)[i]
  is.mito = rownames(sce[[sample_index]]) %in% mtGenes[,1]
  length(is.mito[is.mito == TRUE])
  sce[[sample_index]] = calculateQCMetrics(sce[[sample_index]], feature_controls=list(Mito=is.mito))
  sce[[sample_index]]@int_colData$reducedDims$PCA <- NULL
  sce[[sample_index]] = runPCA(sce[[sample_index]], use_coldata=T)
}

p <- list()
for(i in 1:13){
  sample_index <- list.files(raw_path)[i]
  p[[sample_index]] <- plotPCA(sce[[sample_index]])
}

saveRDS(sce, file = "sce_QCmetrics_calculated.rds")


#Determine QC thresholds
#F8(ZW) - 3.5, 3.0, 20
h <- hist(sce$`F8`$log10_total_counts, 
          breaks=100, col="grey80",
          xlab="Log-total UMI count")
cuts <- cut(h$breaks, c(0,3.5,Inf))
plot(h, col=c("red","white")[cuts])

h <- hist(sce$`F8`$log10_total_features_by_counts, 
          breaks=100, col="grey80",
          xlab="Log-total number of expressed features")
cuts <- cut(h$breaks, c(-0.5,3.0,Inf))
plot(h, col=c("red","white")[cuts])

h <- hist(sce$`F8`$pct_counts_Mito, 
          breaks=100, col="grey80",
          xlab="Proportion of reads in mitochondrial genes")
cuts <- cut(h$breaks, c(20,100,Inf))
plot(h, col=c("red","white")[cuts])


#F9(ZZ) - 4.0, 3.5, 15
h <- hist(sce$`F9`$log10_total_counts, 
          breaks=100, col="grey80",
          xlab="Log-total UMI count")
cuts <- cut(h$breaks, c(0,4.0,Inf))
plot(h, col=c("red","white")[cuts])

h <- hist(sce$`F9`$log10_total_features_by_counts, 
          breaks=100, col="grey80",
          xlab="Log-total number of expressed features")
cuts <- cut(h$breaks, c(-0.5,3.5,Inf))
plot(h, col=c("red","white")[cuts])

h <- hist(sce$`F9`$pct_counts_Mito, 
          breaks=100, col="grey80",
          xlab="Proportion of reads in mitochondrial genes")
cuts <- cut(h$breaks, c(15,100,Inf))
plot(h, col=c("red","white")[cuts])


#F10(ZW) - 4.0, 3.5, 10
h <- hist(sce$`F10`$log10_total_counts, 
          breaks=100, col="grey80",
          xlab="Log-total UMI count")
cuts <- cut(h$breaks, c(0,4.0,Inf))
plot(h, col=c("red","white")[cuts])

h <- hist(sce$`F10`$log10_total_features_by_counts, 
          breaks=100, col="grey80",
          xlab="Log-total number of expressed features")
cuts <- cut(h$breaks, c(-0.5,3.5,Inf))
plot(h, col=c("red","white")[cuts])

h <- hist(sce$`F10`$pct_counts_Mito, 
          breaks=100, col="grey80",
          xlab="Proportion of reads in mitochondrial genes")
cuts <- cut(h$breaks, c(10,100,Inf))
plot(h, col=c("red","white")[cuts])


#F11(ZZ) - 4.0, 3.5, 10
h <- hist(sce$`F11`$log10_total_counts, 
          breaks=100, col="grey80",
          xlab="Log-total UMI count")
cuts <- cut(h$breaks, c(0,4.0,Inf))
plot(h, col=c("red","white")[cuts])

h <- hist(sce$`F11`$log10_total_features_by_counts, 
          breaks=100, col="grey80",
          xlab="Log-total number of expressed features")
cuts <- cut(h$breaks, c(-0.5,3.5,Inf))
plot(h, col=c("red","white")[cuts])

h <- hist(sce$`F11`$pct_counts_Mito, 
          breaks=100, col="grey80",
          xlab="Proportion of reads in mitochondrial genes")
cuts <- cut(h$breaks, c(10,100,Inf))
plot(h, col=c("red","white")[cuts])


#F12(ZZ) - 4.0, 3.5, 10
h <- hist(sce$`F12`$log10_total_counts, 
          breaks=100, col="grey80",
          xlab="Log-total UMI count")
cuts <- cut(h$breaks, c(0,4.0,Inf))
plot(h, col=c("red","white")[cuts])

h <- hist(sce$`F12`$log10_total_features_by_counts, 
          breaks=100, col="grey80",
          xlab="Log-total number of expressed features")
cuts <- cut(h$breaks, c(-0.5,3.5,Inf))
plot(h, col=c("red","white")[cuts])

h <- hist(sce$`F12`$pct_counts_Mito, 
          breaks=100, col="grey80",
          xlab="Proportion of reads in mitochondrial genes")
cuts <- cut(h$breaks, c(10,100,Inf))
plot(h, col=c("red","white")[cuts])


#G1(ZW) - 3.75, 3.5, 5
h <- hist(sce$`G1`$log10_total_counts, 
          breaks=100, col="grey80",
          xlab="Log-total UMI count")
cuts <- cut(h$breaks, c(0,3.75,Inf))
plot(h, col=c("red","white")[cuts])

h <- hist(sce$`G1`$log10_total_features_by_counts, 
          breaks=100, col="grey80",
          xlab="Log-total number of expressed features")
cuts <- cut(h$breaks, c(-0.5,3.5,Inf))
plot(h, col=c("red","white")[cuts])

h <- hist(sce$`G1`$pct_counts_Mito, 
          breaks=100, col="grey80",
          xlab="Proportion of reads in mitochondrial genes")
cuts <- cut(h$breaks, c(5,100,Inf))
plot(h, col=c("red","white")[cuts])


#G2(ZZ) - 4.0, 3.5, 15
h <- hist(sce$`G2`$log10_total_counts, 
          breaks=100, col="grey80",
          xlab="Log-total UMI count")
cuts <- cut(h$breaks, c(0,4.0,Inf))
plot(h, col=c("red","white")[cuts])

h <- hist(sce$`G2`$log10_total_features_by_counts, 
          breaks=100, col="grey80",
          xlab="Log-total number of expressed features")
cuts <- cut(h$breaks, c(-0.5,3.5,Inf))
plot(h, col=c("red","white")[cuts])

h <- hist(sce$`G2`$pct_counts_Mito, 
          breaks=100, col="grey80",
          xlab="Proportion of reads in mitochondrial genes")
cuts <- cut(h$breaks, c(15,100,Inf))
plot(h, col=c("red","white")[cuts])


#G3(ZW) - 4.0, 3.5, 20
h <- hist(sce$`G3`$log10_total_counts, 
          breaks=100, col="grey80",
          xlab="Log-total UMI count")
cuts <- cut(h$breaks, c(0,4.0,Inf))
plot(h, col=c("red","white")[cuts])

h <- hist(sce$`G3`$log10_total_features_by_counts, 
          breaks=100, col="grey80",
          xlab="Log-total number of expressed features")
cuts <- cut(h$breaks, c(-0.5,3.5,Inf))
plot(h, col=c("red","white")[cuts])

h <- hist(sce$`G3`$pct_counts_Mito, 
          breaks=100, col="grey80",
          xlab="Proportion of reads in mitochondrial genes")
cuts <- cut(h$breaks, c(20,100,Inf))
plot(h, col=c("red","white")[cuts])


#G4(ZZ) - 3.5, 3.0, 20
h <- hist(sce$`G4`$log10_total_counts, 
          breaks=100, col="grey80",
          xlab="Log-total UMI count")
cuts <- cut(h$breaks, c(0,3.5,Inf))
plot(h, col=c("red","white")[cuts])

h <- hist(sce$`G4`$log10_total_features_by_counts, 
          breaks=100, col="grey80",
          xlab="Log-total number of expressed features")
cuts <- cut(h$breaks, c(-0.5,3.0,Inf))
plot(h, col=c("red","white")[cuts])

h <- hist(sce$`G4`$pct_counts_Mito, 
          breaks=100, col="grey80",
          xlab="Proportion of reads in mitochondrial genes")
cuts <- cut(h$breaks, c(20,100,Inf))
plot(h, col=c("red","white")[cuts])


#G10(ZZ) - 3.5, 3.0, 10
h <- hist(sce$`G10`$log10_total_counts, 
          breaks=100, col="grey80",
          xlab="Log-total UMI count")
cuts <- cut(h$breaks, c(0,3.5,Inf))
plot(h, col=c("red","white")[cuts])

h <- hist(sce$`G10`$log10_total_features_by_counts, 
          breaks=100, col="grey80",
          xlab="Log-total number of expressed features")
cuts <- cut(h$breaks, c(-0.5,3.0,Inf))
plot(h, col=c("red","white")[cuts])

h <- hist(sce$`G10`$pct_counts_Mito, 
          breaks=100, col="grey80",
          xlab="Proportion of reads in mitochondrial genes")
cuts <- cut(h$breaks, c(10,100,Inf))
plot(h, col=c("red","white")[cuts])


#G11(ZW) - 3.5, 3.0, 10
h <- hist(sce$`G11`$log10_total_counts, 
          breaks=100, col="grey80",
          xlab="Log-total UMI count")
cuts <- cut(h$breaks, c(0,3.5,Inf))
plot(h, col=c("red","white")[cuts])

h <- hist(sce$`G11`$log10_total_features_by_counts, 
          breaks=100, col="grey80",
          xlab="Log-total number of expressed features")
cuts <- cut(h$breaks, c(-0.5,3.0,Inf))
plot(h, col=c("red","white")[cuts])

h <- hist(sce$`G11`$pct_counts_Mito, 
          breaks=100, col="grey80",
          xlab="Proportion of reads in mitochondrial genes")
cuts <- cut(h$breaks, c(10,100,Inf))
plot(h, col=c("red","white")[cuts])s


#G12(ZW) - 3.5, 3.0, 15
h <- hist(sce$`G12`$log10_total_counts, 
          breaks=100, col="grey80",
          xlab="Log-total UMI count")
cuts <- cut(h$breaks, c(0,3.5,Inf))
plot(h, col=c("red","white")[cuts])

h <- hist(sce$`G12`$log10_total_features_by_counts, 
          breaks=100, col="grey80",
          xlab="Log-total number of expressed features")
cuts <- cut(h$breaks, c(-0.5,3.0,Inf))
plot(h, col=c("red","white")[cuts])

h <- hist(sce$`G12`$pct_counts_Mito, 
          breaks=100, col="grey80",
          xlab="Proportion of reads in mitochondrial genes")
cuts <- cut(h$breaks, c(20,100,Inf))
plot(h, col=c("red","white")[cuts])


#########################Using ColData PCA########################################
library(RColorBrewer)

#F8(ZW) - 3.5, 3.0, 20
ggplot(as.data.frame(sce[["F8"]]@int_colData$reducedDims$PCA), 
       aes(x=PC1, y=PC2, color = sce[["F8"]]$log10_total_counts)) +
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

ggplot(as.data.frame(sce[["F8"]]@int_colData$reducedDims$PCA), 
       aes(x=PC1, y=PC2, color = sce[["F8"]]$log10_total_features_by_counts)) +
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

ggplot(as.data.frame(sce[["F8"]]@int_colData$reducedDims$PCA), 
       aes(x=PC1, y=PC2, color = sce[["F8"]]$pct_counts_Mito)) +
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
filtering <- rep(0, ncol(sce[["F8"]]))
filtering[which(sce[["F8"]]$log10_total_counts >= 3.5 & sce[["F8"]]$log10_total_features_by_counts >= 3.0 &
                  sce[["F8"]]$pct_counts_Mito <= 20)] <- 1
sce[["F8"]]$filtering <- as.integer(filtering)
print(table(sce[["F8"]]$filtering))

ggplot(as.data.frame(sce[["F8"]]@int_colData$reducedDims$PCA), 
       aes(x=PC1, y=PC2, color = as.factor(sce[["F8"]]$filtering))) +
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

sce[["F8"]] <- sce[["F8"]][, sce[["F8"]]$filtering == 1]


#F9(ZZ) - 4.0, 3.5, 15
ggplot(as.data.frame(sce[["F9"]]@int_colData$reducedDims$PCA), 
       aes(x=PC1, y=PC2, color = sce[["F9"]]$log10_total_counts)) +
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

ggplot(as.data.frame(sce[["F9"]]@int_colData$reducedDims$PCA), 
       aes(x=PC1, y=PC2, color = sce[["F9"]]$log10_total_features_by_counts)) +
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

ggplot(as.data.frame(sce[["F9"]]@int_colData$reducedDims$PCA), 
       aes(x=PC1, y=PC2, color = sce[["F9"]]$pct_counts_Mito)) +
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
filtering <- rep(0, ncol(sce[["F9"]]))
filtering[which(sce[["F9"]]$log10_total_counts >= 4.0 & sce[["F9"]]$log10_total_features_by_counts >= 3.5 &
                  sce[["F9"]]$pct_counts_Mito <= 15)] <- 1
sce[["F9"]]$filtering <- as.integer(filtering)
print(table(sce[["F9"]]$filtering))

ggplot(as.data.frame(sce[["F9"]]@int_colData$reducedDims$PCA), 
       aes(x=PC1, y=PC2, color = as.factor(sce[["F9"]]$filtering))) +
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

sce[["F9"]] <- sce[["F9"]][, sce[["F9"]]$filtering == 1]


#F10
ggplot(as.data.frame(sce[["F10"]]@int_colData$reducedDims$PCA), 
       aes(x=PC1, y=PC2, color = sce[["F10"]]$log10_total_counts)) +
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

ggplot(as.data.frame(sce[["F10"]]@int_colData$reducedDims$PCA), 
       aes(x=PC1, y=PC2, color = sce[["F10"]]$log10_total_features_by_counts)) +
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

ggplot(as.data.frame(sce[["F10"]]@int_colData$reducedDims$PCA), 
       aes(x=PC1, y=PC2, color = sce[["F10"]]$pct_counts_Mito)) +
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
filtering <- rep(0, ncol(sce[["F10"]]))
filtering[which(sce[["F10"]]$log10_total_counts >= 4.0 & sce[["F10"]]$log10_total_features_by_counts >= 3.5 &
                  sce[["F10"]]$pct_counts_Mito <= 15)] <- 1
sce[["F10"]]$filtering <- as.integer(filtering)
print(table(sce[["F10"]]$filtering))

ggplot(as.data.frame(sce[["F10"]]@int_colData$reducedDims$PCA), 
       aes(x=PC1, y=PC2, color = as.factor(sce[["F10"]]$filtering))) +
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

sce[["F10"]] <- sce[["F10"]][, sce[["F10"]]$filtering == 1]


#F11(ZZ) - 4.0, 3.5, 15
ggplot(as.data.frame(sce[["F11"]]@int_colData$reducedDims$PCA), 
       aes(x=PC1, y=PC2, color = sce[["F11"]]$log10_total_counts)) +
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

ggplot(as.data.frame(sce[["F11"]]@int_colData$reducedDims$PCA), 
       aes(x=PC1, y=PC2, color = sce[["F11"]]$log10_total_features_by_counts)) +
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

ggplot(as.data.frame(sce[["F11"]]@int_colData$reducedDims$PCA), 
       aes(x=PC1, y=PC2, color = sce[["F11"]]$pct_counts_Mito)) +
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
filtering <- rep(0, ncol(sce[["F11"]]))
filtering[which(sce[["F11"]]$log10_total_counts >= 3.5 & sce[["F11"]]$log10_total_features_by_counts >= 3.5 &
                  sce[["F11"]]$pct_counts_Mito <= 15)] <- 1
sce[["F11"]]$filtering <- as.integer(filtering)
print(table(sce[["F11"]]$filtering))

ggplot(as.data.frame(sce[["F11"]]@int_colData$reducedDims$PCA), 
       aes(x=PC1, y=PC2, color = as.factor(sce[["F11"]]$filtering))) +
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

sce[["F11"]] <- sce[["F11"]][, sce[["F11"]]$filtering == 1]


#F12 - 4.0, 3.5, 15
ggplot(as.data.frame(sce[["F12"]]@int_colData$reducedDims$PCA), 
       aes(x=PC1, y=PC2, color = sce[["F12"]]$log10_total_counts)) +
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

ggplot(as.data.frame(sce[["F12"]]@int_colData$reducedDims$PCA), 
       aes(x=PC1, y=PC2, color = sce[["F12"]]$log10_total_features_by_counts)) +
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

ggplot(as.data.frame(sce[["F12"]]@int_colData$reducedDims$PCA), 
       aes(x=PC1, y=PC2, color = sce[["F12"]]$pct_counts_Mito)) +
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
filtering <- rep(0, ncol(sce[["F12"]]))
filtering[which(sce[["F12"]]$log10_total_counts >= 3.5 & sce[["F12"]]$log10_total_features_by_counts >= 3.0 &
                  sce[["F12"]]$pct_counts_Mito <= 10)] <- 1
sce[["F12"]]$filtering <- as.integer(filtering)
print(table(sce[["F12"]]$filtering))

ggplot(as.data.frame(sce[["F12"]]@int_colData$reducedDims$PCA), 
       aes(x=PC1, y=PC2, color = as.factor(sce[["F12"]]$filtering))) +
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

sce[["F12"]] <- sce[["F12"]][, sce[["F12"]]$filtering == 1]


#G1 - 4.0, 3.5, 15
ggplot(as.data.frame(sce[["G1"]]@int_colData$reducedDims$PCA), 
       aes(x=PC1, y=PC2, color = sce[["G1"]]$log10_total_counts)) +
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

ggplot(as.data.frame(sce[["G1"]]@int_colData$reducedDims$PCA), 
       aes(x=PC1, y=PC2, color = sce[["G1"]]$log10_total_features_by_counts)) +
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

ggplot(as.data.frame(sce[["G1"]]@int_colData$reducedDims$PCA), 
       aes(x=PC1, y=PC2, color = sce[["G1"]]$pct_counts_Mito)) +
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
filtering <- rep(0, ncol(sce[["G1"]]))
filtering[which(sce[["G1"]]$log10_total_counts >= 2.5 & sce[["G1"]]$log10_total_features_by_counts >= 2.5 &
                  sce[["G1"]]$pct_counts_Mito <= 20)] <- 1
sce[["G1"]]$filtering <- as.integer(filtering)
print(table(sce[["G1"]]$filtering))

ggplot(as.data.frame(sce[["G1"]]@int_colData$reducedDims$PCA), 
       aes(x=PC1, y=PC2, color = as.factor(sce[["G1"]]$filtering))) +
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

sce[["G1"]] <- sce[["G1"]][, sce[["G1"]]$filtering == 1]

#G2 - 4.0, 3.5, 15
ggplot(as.data.frame(sce[["G2"]]@int_colData$reducedDims$PCA), 
       aes(x=PC1, y=PC2, color = sce[["G2"]]$log10_total_counts)) +
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

ggplot(as.data.frame(sce[["G2"]]@int_colData$reducedDims$PCA), 
       aes(x=PC1, y=PC2, color = sce[["G2"]]$log10_total_features_by_counts)) +
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

ggplot(as.data.frame(sce[["G2"]]@int_colData$reducedDims$PCA), 
       aes(x=PC1, y=PC2, color = sce[["G2"]]$pct_counts_Mito)) +
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
filtering <- rep(0, ncol(sce[["G2"]]))
filtering[which(sce[["G2"]]$log10_total_counts >= 3.0 & sce[["G2"]]$log10_total_features_by_counts >= 2.5 &
                  sce[["G2"]]$pct_counts_Mito <= 15)] <- 1
sce[["G2"]]$filtering <- as.integer(filtering)
print(table(sce[["G2"]]$filtering))

ggplot(as.data.frame(sce[["G2"]]@int_colData$reducedDims$PCA), 
       aes(x=PC1, y=PC2, color = as.factor(sce[["G2"]]$filtering))) +
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

sce[["G2"]] <- sce[["G2"]][, sce[["G2"]]$filtering == 1]


#G3 - 4.0, 3.5, 15
ggplot(as.data.frame(sce[["G3"]]@int_colData$reducedDims$PCA), 
       aes(x=PC1, y=PC2, color = sce[["G3"]]$log10_total_counts)) +
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

ggplot(as.data.frame(sce[["G3"]]@int_colData$reducedDims$PCA), 
       aes(x=PC1, y=PC2, color = sce[["G3"]]$log10_total_features_by_counts)) +
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

ggplot(as.data.frame(sce[["G3"]]@int_colData$reducedDims$PCA), 
       aes(x=PC1, y=PC2, color = sce[["G3"]]$pct_counts_Mito)) +
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
filtering <- rep(0, ncol(sce[["G3"]]))
filtering[which(sce[["G3"]]$log10_total_counts >= 3.5 & sce[["G3"]]$log10_total_features_by_counts >= 3.0 &
                  sce[["G3"]]$pct_counts_Mito <= 25)] <- 1
sce[["G3"]]$filtering <- as.integer(filtering)
print(table(sce[["G3"]]$filtering))

ggplot(as.data.frame(sce[["G3"]]@int_colData$reducedDims$PCA), 
       aes(x=PC1, y=PC2, color = as.factor(sce[["G3"]]$filtering))) +
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

sce[["G3"]] <- sce[["G3"]][, sce[["G3"]]$filtering == 1]

#G4 - 4.0, 3.5, 15
ggplot(as.data.frame(sce[["G4"]]@int_colData$reducedDims$PCA), 
       aes(x=PC1, y=PC2, color = sce[["G4"]]$log10_total_counts)) +
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

ggplot(as.data.frame(sce[["G4"]]@int_colData$reducedDims$PCA), 
       aes(x=PC1, y=PC2, color = sce[["G4"]]$log10_total_features_by_counts)) +
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

ggplot(as.data.frame(sce[["G4"]]@int_colData$reducedDims$PCA), 
       aes(x=PC1, y=PC2, color = sce[["G4"]]$pct_counts_Mito)) +
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
filtering <- rep(0, ncol(sce[["G4"]]))
filtering[which(sce[["G4"]]$log10_total_counts >= 3.5 & sce[["G4"]]$log10_total_features_by_counts >= 3 &
                  sce[["G4"]]$pct_counts_Mito <= 20)] <- 1
sce[["G4"]]$filtering <- as.integer(filtering)
print(table(sce[["G4"]]$filtering))

ggplot(as.data.frame(sce[["G4"]]@int_colData$reducedDims$PCA), 
       aes(x=PC1, y=PC2, color = as.factor(sce[["G4"]]$filtering))) +
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

sce[["G4"]] <- sce[["G4"]][, sce[["G4"]]$filtering == 1]


#G10(ZZ) - 4.0, 3.5, 15
ggplot(as.data.frame(sce[["G10"]]@int_colData$reducedDims$PCA), 
       aes(x=PC1, y=PC2, color = sce[["G10"]]$log10_total_counts)) +
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

ggplot(as.data.frame(sce[["G10"]]@int_colData$reducedDims$PCA), 
       aes(x=PC1, y=PC2, color = sce[["G10"]]$log10_total_features_by_counts)) +
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

ggplot(as.data.frame(sce[["G10"]]@int_colData$reducedDims$PCA), 
       aes(x=PC1, y=PC2, color = sce[["G10"]]$pct_counts_Mito)) +
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
filtering <- rep(0, ncol(sce[["G10"]]))
filtering[which(sce[["G10"]]$log10_total_counts >= 3.5 & sce[["G10"]]$log10_total_features_by_counts >= 3.0 &
                  sce[["G10"]]$pct_counts_Mito <= 10)] <- 1
sce[["G10"]]$filtering <- as.integer(filtering)
print(table(sce[["G10"]]$filtering))

ggplot(as.data.frame(sce[["G10"]]@int_colData$reducedDims$PCA), 
       aes(x=PC1, y=PC2, color = as.factor(sce[["G10"]]$filtering))) +
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

sce[["G10"]] <- sce[["G10"]][, sce[["G10"]]$filtering == 1]


#G11(ZZ) - 4.0, 3.5, 15
ggplot(as.data.frame(sce[["G11"]]@int_colData$reducedDims$PCA), 
       aes(x=PC1, y=PC2, color = sce[["G11"]]$log10_total_counts)) +
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

ggplot(as.data.frame(sce[["G11"]]@int_colData$reducedDims$PCA), 
       aes(x=PC1, y=PC2, color = sce[["G11"]]$log10_total_features_by_counts)) +
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

ggplot(as.data.frame(sce[["G11"]]@int_colData$reducedDims$PCA), 
       aes(x=PC1, y=PC2, color = sce[["G11"]]$pct_counts_Mito)) +
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
filtering <- rep(0, ncol(sce[["G11"]]))
filtering[which(sce[["G11"]]$log10_total_counts >= 3.5 & sce[["G11"]]$log10_total_features_by_counts >= 3 &
                  sce[["G11"]]$pct_counts_Mito <= 10)] <- 1
sce[["G11"]]$filtering <- as.integer(filtering)
print(table(sce[["G11"]]$filtering))

ggplot(as.data.frame(sce[["G11"]]@int_colData$reducedDims$PCA), 
       aes(x=PC1, y=PC2, color = as.factor(sce[["G11"]]$filtering))) +
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

sce[["G11"]] <- sce[["G11"]][, sce[["G11"]]$filtering == 1]


#G12(ZZ) - 4.0, 3.5, 15
ggplot(as.data.frame(sce[["G12"]]@int_colData$reducedDims$PCA), 
       aes(x=PC1, y=PC2, color = sce[["G12"]]$log10_total_counts)) +
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

ggplot(as.data.frame(sce[["G12"]]@int_colData$reducedDims$PCA), 
       aes(x=PC1, y=PC2, color = sce[["G12"]]$log10_total_features_by_counts)) +
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

ggplot(as.data.frame(sce[["G12"]]@int_colData$reducedDims$PCA), 
       aes(x=PC1, y=PC2, color = sce[["G12"]]$pct_counts_Mito)) +
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
filtering <- rep(0, ncol(sce[["G12"]]))
filtering[which(sce[["G12"]]$log10_total_counts >= 3.5 & sce[["G12"]]$log10_total_features_by_counts >= 3.0 &
                  sce[["G12"]]$pct_counts_Mito <= 25)] <- 1
sce[["G12"]]$filtering <- as.integer(filtering)
print(table(sce[["G12"]]$filtering))

ggplot(as.data.frame(sce[["G12"]]@int_colData$reducedDims$PCA), 
       aes(x=PC1, y=PC2, color = as.factor(sce[["G12"]]$filtering))) +
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

sce[["G12"]] <- sce[["G12"]][, sce[["G12"]]$filtering == 1]
