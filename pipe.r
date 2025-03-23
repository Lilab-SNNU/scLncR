library(Seurat)
library(dplyr)
library(ggsci)
library(tidyverse)
library(patchwork)
# library(scPlant)

getwd()
rm(list=ls())
dir.create("~/work/scRNA/step1_data_process/SRR8257103_all")
setwd("~/work/scRNA/step1_data_process/SRR8257103_all")
getwd()

####################################***Load Data***####################################

###***CSV Format***###
# rawcount <- read.csv("./PRJNA471914/PRJNA471914_Root_rawdata_MYB46_Induced.txt", header=T, row.names=1)
# datann <- as(as.matrix(rawcount), "dgCMatrix")
# seuobj <- CreateSeuratObject(counts=datann, min.features=100)

###***10X Matrix***###
data_dir <- '~/work/scRNA/matrix/SRR8257103_all/outs/filtered_feature_bc_matrix'
data <- Read10X(data.dir = data_dir)
ath_data <- CreateSeuratObject(counts = data, project = "SRR8257103_all")

####################################***Process QC***####################################
ath_data[["percent.mt"]] <- PercentageFeatureSet(ath_data, pattern = "^MT-")
pdf("SRR8257103_all.pdf")

# Visualize QC metrics as a violin plot
plot <- VlnPlot(ath_data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ggsave('FeatureQC.pdf', plot=plot, width = 8, height = 6)
plot

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(ath_data, feature1 = "nCount_RNA", feature2 = "percent.mt")+scale_color_npg()
plot2 <- FeatureScatter(ath_data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")+scale_color_npg()
ggsave('FeatureScatter_1.pdf', plot=plot1, width = 8, height = 6)
ggsave('FeatureScatter_2.pdf', plot=plot2, width = 8, height = 6)
plot1
plot2

ath_data <- subset(ath_data, subset = nFeature_RNA > 100 & nFeature_RNA < 7000 & percent.mt < 10)  
ath_data <- NormalizeData(ath_data, normalization.method = "LogNormalize", scale.factor = 10000)

####################################***PCA Process***####################################
#Identification of highly variable features (feature selection)
ath_data <- FindVariableFeatures(ath_data, selection.method = "vst", nfeatures = 3000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(ath_data), 10)
# plot variable features with and without labels

plot1 <- VariableFeaturePlot(ath_data)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
ggsave('VariableFeaturePlot.pdf', plot=plot1, width = 8, height = 6)
ggsave('LabelPoints.pdf', plot=plot2, width = 8, height = 6)
plot1
plot2

# Scaling the data 
all.genes <- rownames(ath_data)
ath_data <- ScaleData(ath_data, features = all.genes)

#Perform linear dimensional reductionPerform linear dimensional reduction
ath_data <- RunPCA(ath_data, features = VariableFeatures(object = ath_data))

# Get the feature loadings for a given DimReduc
t(Loadings(object = ath_data[["pca"]])[1:5,1:5])
# Get the feature loadings for a specified DimReduc in a Seurat object
t(Loadings(object = ath_data, reduction = "pca")[1:5,1:5])
# Set the feature loadings for a given DimReduc
new.loadings <- Loadings(object = ath_data[["pca"]])
new.loadings <- new.loadings + 0.01
Loadings(object = ath_data[["pca"]]) <- new.loadings
VizDimLoadings(ath_data)
p <- DimPlot(ath_data, reduction = "pca")
ggsave('DimPlot.pdf', plot=p, width = 8, height = 6)
p
p <- DimHeatmap(ath_data, dims = 1:15, cells = 500, balanced = TRUE)
ggsave('DimHeatmap.pdf', plot=p, width = 8, height = 6)
p
# Determine the ‘dimensionality’ of the dataset 


####################################***Cluster Process***####################################
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
ath_data <- JackStraw(ath_data, num.replicate = 100)
ath_data <- ScoreJackStraw(ath_data, dims = 1:20)
plot1<-JackStrawPlot(ath_data, dims = 1:15)
plot2<-ElbowPlot(ath_data)
ggsave('JackStrawPlot.pdf', plot=plot1, width = 8, height = 6)
ggsave('ElbowPlot.pdf', plot=plot2, width = 8, height = 6)
plot1
plot2
# Cluster the cells 
#Identify clusters of cells by a shared nearest neighbor (SNN) modularity optimization based clustering algorithm. First calculate k-nearest neighbors and construct the SNN graph. Then optimize the modularity function to determine clusters. For a full description of the algorithms, see Waltman and van Eck (2013) The European Physical Journal B. Thanks to Nigel Delaney (evolvedmicrobe@github) for the rewrite of the Java modularity optimizer code in Rcpp!
ath_data <- FindNeighbors(ath_data, dims = 1:10)
ath_data <- FindClusters(ath_data, resolution = 0.5)

#Constructs a phylogenetic tree relating the 'average' cell from each identity class. 
# Tree is estimated based on a distance matrix constructed in either gene expression space or PCA spac

ath_data<-BuildClusterTree(ath_data)
Tool(object = ath_data, slot = 'BuildClusterTree')

plot <- PlotClusterTree(ath_data)
ggsave('ClusterTreePlot.pdf', plot=plot, width = 8, height = 6)
plot
#Calculate the Barcode Distribution Inflection
ath_data<-CalculateBarcodeInflections(ath_data)
SubsetByBarcodeInflections(ath_data)

#Run non-linear dimensional reduction (UMAP/tSNE)
# If you haven't installed UMAP, you can do so via reticulate::py_install(packages =
# 'umap-learn')
ath_data <- RunUMAP(ath_data, dims = 1:10)

#### head(ath_data@reductions$umap@cell.embeddings) # 提取UMAP坐标值。

ath_data <- RunTSNE(ath_data, dims = 1:10)


# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
plot1<-DimPlot(ath_data, reduction = "umap",label = TRUE)+scale_color_npg()
plot2<-DimPlot(ath_data, reduction = "tsne",label = TRUE)+scale_color_npg()
# CombinePlots(plots = list(plot1, plot2),legend="bottom")
ggsave("umap.pdf", plot=plot1, width = 8, height = 6)
ggsave("tsne.pdf", plot=plot2, width = 8, height = 6)
plot1
plot2
library(gridExtra) 

####################################***Marker Gene Filter***####################################
# Finding differentially expressed features (cluster biomarkers) 
# find all markers of cluster 1
cluster1.markers <- FindMarkers(ath_data, ident.1 = 1, min.pct = 0.25)
head(cluster1.markers, n = 5)

ath_data.markers <- FindAllMarkers(ath_data, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

head(ath_data.markers, n = 5)

top2 <- ath_data.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)

# ath_data[['groups']] <- sample(x = c('g1', 'g2'), size = ncol(x = ath_data), replace = TRUE)
# head(FindConservedMarkers(ath_data, ident.1 = 0, ident.2 = 1, grouping.var = "groups"), n = 10)

# plot1<-VlnPlot(ath_data, features = c("AT1G12080", "AT4G11210"))+scale_color_npg()
# plot1
# you can plot raw counts as well
# plot2<- VlnPlot(ath_data, features = c("AT1G12080", "AT4G11210"),ncol=1, same.y.lims=T,slot = "counts", log = TRUE)+scale_color_npg()
# plot2
# CombinePlots(plots = list(plot1, plot2))

plot1<- FeaturePlot(ath_data, features = top2$gene, min.cutoff = 0, max.cutoff = 4)
top10 <- ath_data.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
plot2<- DoHeatmap(ath_data, features = top10$gene) + NoLegend()+scale_color_npg()
ggsave("FeaturePlot.pdf", plot=plot1, width = 20, height = 30)
ggsave("DoHeatmap.pdf", plot=plot2, width = 8, height = 6)
plot1
plot2


####################################***Cell Annatation***####################################
# singleR 
#Assigning cell type identity to clusters
# new.cluster.ids <- c("aa", "ss", "dd", "ff", "gg", "hh", "jj", "kk", "ll", "QQ", "WW", "EE")
# names(new.cluster.ids) <- levels(ath_data)
# ath_data <- RenameIdents(ath_data, new.cluster.ids)
# plot1<-DimPlot(ath_data, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
# plot2<-DimPlot(ath_data, reduction = "tsne", label = TRUE, pt.size = 0.5) + NoLegend()
# ggsave("umap_new_name.pdf", plot=plot1, width = 8, height = 6)
# ggsave("tsne_new_name.pdf", plot=plot2, width = 8, height = 6)
# plot1
# plot2


####################################***Conection Analysis***####################################
AverageExp<-AverageExpression(ath_data,features=unique(top10$gene))
typeof(AverageExp)
head(AverageExp$RNA)

library(psych)
library(pheatmap)
coorda<-corr.test(AverageExp$RNA,AverageExp$RNA,method="spearman")
plot <- pheatmap(coorda$r)
ggsave("pheatmap.pdf", plot=plot, width = 8, height = 6)
dev.off()
ath_data
saveRDS(ath_data, file = "./R_result.rds")

getwd()
