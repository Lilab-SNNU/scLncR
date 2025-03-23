library(Seurat)
library(dplyr)
library(ggsci)
library(tidyverse)
library(patchwork)
# library(scPlant)

getwd()
rm(list=ls())
dir.create("/home/data/scLncR/test_0707/step3/tmp")
setwd("/home/data/scLncR/test_0707/step3/tmp")
getwd()



load_sc <- function(data_dir, project_name){
    data <- Read10X(data.dir = data_dir)
    ath_lnc <- CreateSeuratObject(counts = data, project = project_name)

####################################***Process QC***####################################
    ath_lnc[["percent.mt"]] <- PercentageFeatureSet(ath_lnc, pattern = "^MT-")
    ath_lnc[["percent.mt"]] <- PercentageFeatureSet(ath_lnc, pattern = "^MT-")
    ath_lnc <- subset(ath_lnc, subset = nFeature_RNA > 100 & nFeature_RNA < 7000 & percent.mt < 10) 
}

data_00 <- '/home/data/scLncR/test_0707/SRR8257100/outs/filtered_feature_bc_matrix'
data_01 <- '/home/data/scLncR/test_0707/SRR8257101/outs/filtered_feature_bc_matrix'
data_02 <- '/home/data/scLncR/test_0707/SRR8257102/outs/filtered_feature_bc_matrix'
data_03 <- '/home/data/scLncR/test_0707/SRR8257103/outs/filtered_feature_bc_matrix'

ath_list <- list()
ath_list[["sample_00"]] <- load_sc(data_00, "SRR8257100")
ath_list[["sample_01"]] <- load_sc(data_01, "SRR8257101")
ath_list[["sample_02"]] <- load_sc(data_02, "SRR8257102")
ath_list[["sample_03"]] <- load_sc(data_03, "SRR8257103")

for (i in 1:length(ath_list)) {
  ath_list[[i]] <- NormalizeData(ath_list[[i]], verbose = F)
  ath_list[[i]] <- FindVariableFeatures(ath_list[[i]], selection.method = "vst", nfeatures = 4000, verbose = F)
}


features <- SelectIntegrationFeatures(object.list=ath_list)
ath_anchors <- FindIntegrationAnchors(object.list=ath_list, anchor.features=features)
ath_combined <- IntegrateData(anchorset=ath_anchors)
DefaultAssay(ath_combined) <- "integrated"
ath_combined


all.genes <- rownames(ath_combined)
ath_lnc <- ScaleData(ath_combined, features = all.genes, verbose=FALSE)
ath_lnc <- RunPCA(ath_lnc, npcs=30, features = VariableFeatures(object = ath_lnc), verbose=FALSE)

# Get the feature loadings for a given DimReduc
t(Loadings(object = ath_lnc[["pca"]])[1:5,1:5])
# Get the feature loadings for a specified DimReduc in a Seurat object
t(Loadings(object = ath_lnc, reduction = "pca")[1:5,1:5])
# Set the feature loadings for a given DimReduc
new.loadings <- Loadings(object = ath_lnc[["pca"]])
new.loadings <- new.loadings + 0.01
Loadings(object = ath_lnc[["pca"]]) <- new.loadings
VizDimLoadings(ath_lnc)
p <- DimPlot(ath_lnc, reduction = "pca")
ggsave('DimPlot.pdf', plot=p, width = 8, height = 6)
p
p <- DimHeatmap(ath_lnc, dims = 1:15, cells = 500, balanced = TRUE)
ggsave('DimHeatmap.pdf', plot=p, width = 8, height = 6)
p


ath_lnc <- JackStraw(ath_lnc, num.replicate = 100)
ath_lnc <- ScoreJackStraw(ath_lnc, dims = 1:20)
plot1<-JackStrawPlot(ath_lnc, dims = 1:20)
plot2<-ElbowPlot(ath_lnc)
ggsave('JackStrawPlot.pdf', plot=plot1, width = 8, height = 6)
ggsave('ElbowPlot.pdf', plot=plot2, width = 8, height = 6)
plot1
plot2


ath_lnc <- FindNeighbors(ath_lnc, dims = 1:30)
ath_lnc <- FindClusters(ath_lnc, resolution = 0.4)

ath_lnc<-BuildClusterTree(ath_lnc)
Tool(object = ath_lnc, slot = 'BuildClusterTree')

plot <- PlotClusterTree(ath_lnc)
ggsave('ClusterTreePlot.pdf', plot=plot, width = 8, height = 6)
plot

ath_lnc<-CalculateBarcodeInflections(ath_lnc)
# SubsetByBarcodeInflections(ath_lnc)

ath_lnc <- RunUMAP(ath_lnc, dims = 1:30)

# head(ath_lnc@reductions$umap@cell.embeddings) # 提取UMAP坐标值。

ath_lnc <- RunTSNE(ath_lnc, dims = 1:30)


# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
plot1<-DimPlot(ath_lnc, reduction = "umap", group.by='orig.ident') + scale_color_npg()
plot2<-DimPlot(ath_lnc, reduction = "umap", label = TRUE, repel=TRUE) + scale_color_npg()
plot3<-DimPlot(ath_lnc, reduction = "umap", split.by='orig.ident', label = TRUE, repel=TRUE) + scale_color_npg()
ggsave("umap_sample.pdf", plot=plot1, width = 8, height = 6)
ggsave("umap_CD8.pdf", plot=plot2, width = 8, height = 6)
ggsave("umap_split.pdf", plot=plot3, width = 8, height = 4)
plot1
plot2
plot3


plot1<-DimPlot(ath_lnc, reduction = "tsne", group.by='orig.ident') + scale_color_npg()
plot2<-DimPlot(ath_lnc, reduction = "tsne", label = TRUE, repel=TRUE) + scale_color_npg()
plot3<-DimPlot(ath_lnc, reduction = "tsne", split.by='orig.ident', label = TRUE, repel=TRUE) + scale_color_npg()
ggsave("tsne_sample.pdf", plot=plot1, width = 8, height = 6)
ggsave("tsne_CD8.pdf", plot=plot2, width = 8, height = 6)
ggsave("tsne_split.pdf", plot=plot3, width = 16, height = 4)
plot1
plot2
plot3
library(gridExtra) 

cluster1.markers <- FindMarkers(ath_lnc, ident.1 = 1, min.pct = 0.25)
head(cluster1.markers, n = 5)

ath_lnc.markers <- FindAllMarkers(ath_lnc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
head(ath_lnc.markers, n = 5)
write.csv(ath_lnc.markers, file = "Ath_Marker.csv")
top10 <- ath_lnc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
top2 <- ath_lnc.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
plot1<- FeaturePlot(ath_lnc, features = top2$gene, min.cutoff = 0, max.cutoff = 4)
plot2<- DoHeatmap(ath_lnc, features = top10$gene) + NoLegend()+scale_color_npg()
ggsave("FeaturePlot.pdf", plot=plot1, width = 20, height = 30)
ggsave("DoHeatmap.pdf", plot=plot2, width = 8, height = 6)
plot1
plot2

AverageExp<-AverageExpression(ath_lnc,features=unique(top10$gene))
typeof(AverageExp)
head(AverageExp$RNA)

library(psych)
library(pheatmap)
coorda<-corr.test(AverageExp$RNA,AverageExp$RNA,method="spearman")
plot <- pheatmap(coorda$r)
ggsave("pheatmap.pdf", plot=plot, width = 8, height = 6)
dev.off()
saveRDS(ath_lnc, file = "./R_result.rds")
