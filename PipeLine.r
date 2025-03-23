library(Seurat)
library(dplyr)
library(ggsci)
library(tidyverse)
library(patchwork)

getwd()
rm(list=ls())


########################################***Cut Data From All Matrix***########################################
cut_cell_matrix <- function(data_matrix, pattern, project_name){
    row_indexs <- grep(pattern, rownames(data_matrix))
    data_cut <- data_matrix[row_indexs, ]
    seu_obj <- CreateSeuratObject(counts = data_cut, project = project_name)
    return(seu_obj)
}


#############################################***Data Normalize ***#############################################
data_normalize <- function(seu_obj, filter=FALSE){
    seu_obj[["percent.mt"]] <- PercentageFeatureSet(seu_obj, pattern = "^MT-")
    
    if (filter){
        pdf("filter.pdf")

        # Visualize QC metrics as a violin plot
        plot <- VlnPlot(seu_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
        ggsave('FeatureQC.pdf', plot=plot, width = 8, height = 6)
        plot

        # FeatureScatter is typically used to visualize feature-feature relationships, but can be used
        # for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
        plot1 <- FeatureScatter(seu_obj, feature1 = "nCount_RNA", feature2 = "percent.mt")+scale_color_npg()
        plot2 <- FeatureScatter(seu_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")+scale_color_npg()
        ggsave('FeatureScatter_1.pdf', plot=plot1, width = 8, height = 6)
        ggsave('FeatureScatter_2.pdf', plot=plot2, width = 8, height = 6)
        plot1
        plot2
        dev.off()
        seu_obj <- subset(seu_obj, subset = nFeature_RNA > 100 & nFeature_RNA < 7000 & percent.mt < 10)
    }
    seu_obj <- NormalizeData(seu_obj, normalization.method = "LogNormalize", scale.factor = 10000)
    return(seu_obj)
}


#############################################***Data Process***#############################################
data_process <- function(seu_obj, res_dir="", merge=FALSE){
    dir.create(res_dir)
    setwd(res_dir)
    getwd()
    pdf("res.pdf")
    if(merge){
        ####################################***PCA Process***####################################
        #Identification of highly variable features (feature selection)
        seu_obj <- FindVariableFeatures(seu_obj, selection.method = "vst", nfeatures = 3000)

        # Identify the 10 most highly variable genes
        top10 <- head(VariableFeatures(seu_obj), 10)
        # plot variable features with and without labels

        plot1 <- VariableFeaturePlot(seu_obj)
        plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
        ggsave('VariableFeaturePlot.pdf', plot=plot1, width = 8, height = 6)
        ggsave('LabelPoints.pdf', plot=plot2, width = 8, height = 6)
        plot1
        plot2
    }

    # Scaling the data 
    all.genes <- rownames(seu_obj)
    seu_obj <- ScaleData(seu_obj, features = all.genes)

    #Perform linear dimensional reductionPerform linear dimensional reduction
    seu_obj <- RunPCA(seu_obj, features = VariableFeatures(object = seu_obj))

    # Get the feature loadings for a given DimReduc
    t(Loadings(object = seu_obj[["pca"]])[1:5,1:5])
    # Get the feature loadings for a specified DimReduc in a Seurat object
    t(Loadings(object = seu_obj, reduction = "pca")[1:5,1:5])
    # Set the feature loadings for a given DimReduc
    new.loadings <- Loadings(object = seu_obj[["pca"]])
    new.loadings <- new.loadings + 0.01
    Loadings(object = seu_obj[["pca"]]) <- new.loadings
    VizDimLoadings(seu_obj)
    p <- DimPlot(seu_obj, reduction = "pca")
    ggsave('DimPlot.pdf', plot=p, width = 8, height = 6)
    p
    p <- DimHeatmap(seu_obj, dims = 1:15, cells = 500, balanced = TRUE)
    ggsave('DimHeatmap.pdf', plot=p, width = 8, height = 6)
    p
    # Determine the ‘dimensionality’ of the dataset 


    ####################################***Cluster Process***####################################
    # NOTE: This process can take a long time for big datasets, comment out for expediency. More
    # approximate techniques such as those implemented in ElbowPlot() can be used to reduce
    # computation time
    seu_obj <- JackStraw(seu_obj, num.replicate = 100)
    seu_obj <- ScoreJackStraw(seu_obj, dims = 1:20)
    plot1<-JackStrawPlot(seu_obj, dims = 1:15)
    plot2<-ElbowPlot(seu_obj)
    ggsave('JackStrawPlot.pdf', plot=plot1, width = 8, height = 6)
    ggsave('ElbowPlot.pdf', plot=plot2, width = 8, height = 6)
    plot1
    plot2
    # Cluster the cells 
    #Identify clusters of cells by a shared nearest neighbor (SNN) modularity optimization based clustering algorithm. First calculate k-nearest neighbors and construct the SNN graph. Then optimize the modularity function to determine clusters. For a full description of the algorithms, see Waltman and van Eck (2013) The European Physical Journal B. Thanks to Nigel Delaney (evolvedmicrobe@github) for the rewrite of the Java modularity optimizer code in Rcpp!
    seu_obj <- FindNeighbors(seu_obj, dims = 1:10)
    seu_obj <- FindClusters(seu_obj, resolution = 0.5)

    #Constructs a phylogenetic tree relating the 'average' cell from each identity class. 
    # Tree is estimated based on a distance matrix constructed in either gene expression space or PCA spac

    seu_obj<-BuildClusterTree(seu_obj)
    Tool(object = seu_obj, slot = 'BuildClusterTree')

    plot <- PlotClusterTree(seu_obj)
    ggsave('ClusterTreePlot.pdf', plot=plot, width = 8, height = 6)
    plot
    #Calculate the Barcode Distribution Inflection
    seu_obj<-CalculateBarcodeInflections(seu_obj)
    SubsetByBarcodeInflections(seu_obj)

    #Run non-linear dimensional reduction (UMAP/tSNE)
    # If you haven't installed UMAP, you can do so via reticulate::py_install(packages =
    # 'umap-learn')
    seu_obj <- RunUMAP(seu_obj, dims = 1:10)

    #### head(seu_obj@reductions$umap@cell.embeddings) # 提取UMAP坐标值。

    seu_obj <- RunTSNE(seu_obj, dims = 1:10)


    # note that you can set `label = TRUE` or use the LabelClusters function to help label
    # individual clusters
    plot1<-DimPlot(seu_obj, reduction = "umap",label = TRUE)+scale_color_npg()
    plot2<-DimPlot(seu_obj, reduction = "tsne",label = TRUE)+scale_color_npg()
    # CombinePlots(plots = list(plot1, plot2),legend="bottom")
    ggsave("umap.pdf", plot=plot1, width = 8, height = 6)
    ggsave("tsne.pdf", plot=plot2, width = 8, height = 6)
    plot1
    plot2
    library(gridExtra) 

    ####################################***Marker Gene Filter***####################################
    # Finding differentially expressed features (cluster biomarkers) 
    # find all markers of cluster 1
    cluster1.markers <- FindMarkers(seu_obj, ident.1 = 1, min.pct = 0.25)
    head(cluster1.markers, n = 5)

    seu_obj.markers <- FindAllMarkers(seu_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

    head(seu_obj.markers, n = 5)

    top2 <- seu_obj.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)

    # seu_obj[['groups']] <- sample(x = c('g1', 'g2'), size = ncol(x = seu_obj), replace = TRUE)
    # head(FindConservedMarkers(seu_obj, ident.1 = 0, ident.2 = 1, grouping.var = "groups"), n = 10)

    # plot1<-VlnPlot(seu_obj, features = c("AT1G12080", "AT4G11210"))+scale_color_npg()
    # plot1
    # you can plot raw counts as well
    # plot2<- VlnPlot(seu_obj, features = c("AT1G12080", "AT4G11210"),ncol=1, same.y.lims=T,slot = "counts", log = TRUE)+scale_color_npg()
    # plot2
    # CombinePlots(plots = list(plot1, plot2))

    plot1<- FeaturePlot(seu_obj, features = top2$gene, min.cutoff = 0, max.cutoff = 4)
    top10 <- seu_obj.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
    plot2<- DoHeatmap(seu_obj, features = top10$gene) + NoLegend()+scale_color_npg()
    ggsave("FeaturePlot.pdf", plot=plot1, width = 20, height = 30)
    ggsave("DoHeatmap.pdf", plot=plot2, width = 8, height = 6)
    plot1
    plot2
    AverageExp<-AverageExpression(seu_obj,features=unique(top10$gene))
    typeof(AverageExp)
    head(AverageExp$RNA)

    library(psych)
    library(pheatmap)
    coorda<-corr.test(AverageExp$RNA,AverageExp$RNA,method="spearman")
    plot <- pheatmap(coorda$r)
    ggsave("pheatmap.pdf", plot=plot, width = 8, height = 6)
    dev.off()
    return(seu_obj)
}


#############################################***Get Renew Data***#############################################
pipe <- function(seu_obj, sample_name){
    df_count <- seu_obj@assays$RNA@counts
    seu_lnc <- cut_cell_matrix(data_matrix=df_count, pattern="^AthLNC.*", project_name=paste(sample_name, "lnc", sep="_"))
    seu_lnc <- data_normalize(seu_lnc, filter=FALSE)
    seu_gene <- cut_cell_matrix(data_matrix=df_count, pattern="^AT.*", project_name=paste(sample_name, "gene", sep="_"))
    seu_gene <- data_normalize(seu_gene, filter=FALSE)
    df_bind <- rbind(seu_gene@assays$RNA@data, seu_lnc@assays$RNA@data)
    seu_bind <- seu_obj
    seu_bind@assays$RNA@data <- df_bind
    return(seu_bind)
}


#############################################***Data Merge***#############################################
pipe_merge <- function(dir1, dir2, dir3, dir4){
    
    load_sc <- function(data_dir, project_name){
        data <- Read10X(data.dir = data_dir)
        seu_obj <- CreateSeuratObject(counts = data, project = project_name)
        seu_obj[["percent.mt"]] <- PercentageFeatureSet(seu_obj, pattern = "^MT-")
        seu_obj <- subset(seu_obj, subset = nFeature_RNA > 100 & nFeature_RNA < 7000 & percent.mt < 10)
        seu_obj <- NormalizeData(seu_obj, normalization.method = "LogNormalize", scale.factor = 10000)
        return(seu_obj)
    }
    
    data_merge <- function(data_list){
        features <- SelectIntegrationFeatures(object.list=data_list)
        data_anchors <- FindIntegrationAnchors(object.list=data_list, anchor.features=features)
        data_combined <- IntegrateData(anchorset=data_anchors)
        DefaultAssay(data_combined) <- "integrated"
        data_combined
        return(data_combined)
    }
    dir.create("merge_13_16")
    setwd("merge_13_16")
    path <- getwd()
    project_name <- c("sample_13", "sample_14", "sample_15", "sample_16")
    ath_list <- list()
    ath_list[["sample_13"]] <- load_sc(dir1, "sample_13")
    ath_list[["sample_14"]] <- load_sc(dir2, "sample_14")
    ath_list[["sample_15"]] <- load_sc(dir3, "sample_15")
    ath_list[["sample_16"]] <- load_sc(dir4, "sample_16")
    
    renew_list <- list()
    for (i in 1:length(ath_list)) {
        ath_renew <- pipe(ath_list[[i]], sample_name=project_name[i])
        ath_list[[i]] <- FindVariableFeatures(ath_list[[i]], selection.method = "vst", nfeatures = 4000, verbose = F)
        ath_renew <- FindVariableFeatures(ath_renew, selection.method = "vst", nfeatures = 4000, verbose = F)
        renew_list[[i]] <- ath_renew
    }
    
    ath_combined <- data_merge(ath_list)
    renew_combined <- data_merge(renew_list)
    
    ath_combined <- data_process(ath_combined, res_dir="old")
    renew_combined <- data_process(renew_combined, res_dir="new")
    saveRDS(ath_combined, file = "./R_result.rds")
    saveRDS(renew_combined, file = "./R_result_renew.rds")
}

data_00 <- '/home/li/anwser/ysw/work/scRNA/matrix/GSE141730_all/SRR10620013'
data_01 <- '/home/li/anwser/ysw/work/scRNA/matrix/GSE141730_all/SRR10620014'
data_02 <- '/home/li/anwser/ysw/work/scRNA/matrix/GSE141730_all/SRR10620015'
data_03 <- '/home/li/anwser/ysw/work/scRNA/matrix/GSE141730_all/SRR10620016'

pipe_merge(data_00, data_01, data_02, data_03)
# dir.create("~/work/scRNA/step1_data_process/GSE141730_all/SRR8257103")
# setwd("~/work/scRNA/step1_data_process/GSE141730_all/SRR8257103")

# data_dir <- '~/work/scRNA/matrix/SRR8257103_all/outs/filtered_feature_bc_matrix'
# data <- Read10X(data.dir = data_dir)
# ath_data <- CreateSeuratObject(counts = data, project = "SRR8257103_all")
# ath_data <- data_normalize(ath_data, filter=TRUE)
# ath_bind <- pipe(ath_data, sample_name="SRR8257103")
# ath_data <- data_process(ath_data, res_dir="old")
# ath_bind <- data_process(ath_bind, res_dir="new")


# saveRDS(ath_data, file = "./R_result.rds")
# saveRDS(ath_bind, file = "./R_result_renew.rds")
