library(getopt)

spec <- matrix(
    c("sample_dir",  "d", 2, "character", "The dirs contant count matrix of samples, you should use the absolute path. \
                        If you have multiple sample data, you can use commas to separate them.",
      "lnc_name", "n", 2, "character", "LncRNA name pre, same as with step1",
      "output_path",  "o", 2, "character",  "Output dir of res",
      "help",   "h", 0, "logical",  "This is Help!"),
    byrow=TRUE, ncol=5)


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

data_process <- function(seu_obj, res_dir="", lnc=TRUE, lnc_name=lnc_name){
    dir.create(res_dir)
    setwd(res_dir)
    getwd()
    
    ####################################***PCA Process***####################################
    #Identification of highly variable features (feature selection)
    seu_obj <- FindVariableFeatures(seu_obj, selection.method = "vst", nfeatures = 6000)
    print(seu_obj)

    # Identify the 10 most highly variable genes
    top10 <- head(VariableFeatures(seu_obj), 10)
    # plot variable features with and without labels

    plot1 <- VariableFeaturePlot(seu_obj)
    plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
    ggsave('VariableFeaturePlot.pdf', plot=plot1, width = 8, height = 6)
    ggsave('LabelPoints.pdf', plot=plot2, width = 8, height = 6)

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
    
    plot <- DimPlot(seu_obj, reduction = "pca")
    ggsave('DimPlot.pdf', plot=plot, width = 8, height = 6)
    
    plot <- DimHeatmap(seu_obj, dims = 1:15, cells = 500, balanced = TRUE)
    ggsave('DimHeatmap.pdf', plot=plot, width = 8, height = 6)

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


    # Cluster the cells 
    # Identify clusters of cells by a shared nearest neighbor (SNN) modularity optimization based clustering algorithm. 
    # First calculate k-nearest neighbors and construct the SNN graph. Then optimize the modularity function to determine clusters. 
    # For a full description of the algorithms, see Waltman and van Eck (2013) The European Physical Journal B. 
    seu_obj <- FindNeighbors(seu_obj, dims = 1:10)
    seu_obj <- FindClusters(seu_obj, resolution = 0.5)

    #Constructs a phylogenetic tree relating the 'average' cell from each identity class. 
    # Tree is estimated based on a distance matrix constructed in either gene expression space or PCA spac

    seu_obj<-BuildClusterTree(seu_obj)
    Tool(object = seu_obj, slot = 'BuildClusterTree')

    plot <- PlotClusterTree(seu_obj)

    # Calculate the Barcode Distribution Inflection
    seu_obj <- CalculateBarcodeInflections(seu_obj)
    SubsetByBarcodeInflections(seu_obj)
    
    # Run non-linear dimensional reduction (UMAP/tSNE)
      
    seu_obj <- RunUMAP(seu_obj, dims = 1:10)

    seu_obj <- RunTSNE(seu_obj, dims = 1:10)



    # note that you can set `label = TRUE` or use the LabelClusters function to help label
    # individual clusters
    plot1<-DimPlot(seu_obj, reduction = "umap",label = TRUE)+scale_color_npg()
    plot2<-DimPlot(seu_obj, reduction = "tsne",label = TRUE)+scale_color_npg()
    # CombinePlots(plots = list(plot1, plot2),legend="bottom")
    ggsave("umap.pdf", plot=plot1, width = 8, height = 6)
    ggsave("tsne.pdf", plot=plot2, width = 8, height = 6)
    
    library(gridExtra) 

    ####################################***Marker Gene Filter***####################################
    # Finding differentially expressed features (cluster biomarkers) 
    # find all markers of cluster 1

    seu_obj.markers <- FindAllMarkers(seu_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

    if (lnc){
        lnc_trans <- grep("^AthLnc.*", seu_obj.markers$gene)
        lncRNAs <- seu_obj.markers[lnc_trans, ]
        seu_obj@assays$sig_lncRNA <- lncRNAs
        write.csv(lncRNAs, file="markers_lncRNAs_matrix.csv")
        
        lncs <- rownames(lncRNAs)
        plot1<- FeaturePlot(seu_obj, features = lncs, min.cutoff = 0, max.cutoff = 4)
        plot2<- DoHeatmap(seu_obj, features = lncs) + NoLegend()+scale_color_npg()
        ggsave("FeaturePlot_lncRNAs.pdf", plot=plot1, width = 20, height = 30)
        ggsave("DoHeatmap_lncRNAs.pdf", plot=plot2, width = 8, height = 6)
    }
    print("\n###########################***Debug***##########################")
    print("Find Marker Done")
    top2 <- seu_obj.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
    top10 <- seu_obj.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
    
    print("\n###########################***Debug***##########################")
    print("Top Gene Choose Done")
    saveRDS(seu_obj, file = "./tmp_result.rds")
    plot1<- FeaturePlot(seu_obj, features = top2$gene, min.cutoff = 0, max.cutoff = 4)
    ggsave("FeaturePlot.pdf", plot=plot1, width = 20, height = 30)


    AverageExp<-AverageExpression(seu_obj,features=unique(top10$gene))
    typeof(AverageExp)
    head(AverageExp$RNA)

    library(psych)
    library(pheatmap)
    coorda<-corr.test(AverageExp$RNA,AverageExp$RNA,method="spearman")
    plot <- pheatmap(coorda$r)
    ggsave("pheatmap.pdf", plot=plot, width = 8, height = 6)
    
    return(seu_obj)
}

#############################################***Get Renew Data***#############################################

pipe <- function(seu_obj, sample_name){

    print("\n#############################################***Get Renew Data***#############################################\n")
    df_count <- seu_obj@assays$RNA@counts
    seu_lnc <- cut_cell_matrix(data_matrix=df_count, pattern="^AthLnc.*", project_name=paste(sample_name, "lnc", sep="_"))
    seu_lnc <- data_normalize(seu_lnc, filter=FALSE)
    print(seu_lnc)
    seu_gene <- cut_cell_matrix(data_matrix=df_count, pattern="^AT.*", project_name=paste(sample_name, "gene", sep="_"))
    seu_gene <- data_normalize(seu_gene, filter=FALSE)
    print(seu_gene)
    df_bind <- rbind(seu_gene@assays$RNA@data, seu_lnc@assays$RNA@data)
    seu_bind <- seu_obj
    seu_bind@assays$RNA@data <- df_bind
    print("\n#############################################***Get Renew Data Done***#############################################\n")
    return(seu_bind)
}

#############################################***Data Merge***#############################################

pipe_merge <- function(sample_dirs, project_names, res_dir){

    library(Seurat)
    library(dplyr)
    library(ggsci)
    library(tidyverse)
    library(patchwork)
    
    load_sc <- function(data_dir, project_name){
        print(data_dir)
        data <- Read10X(data.dir = data_dir)
        seu_obj <- CreateSeuratObject(counts = data, project = project_name)
        seu_obj[["percent.mt"]] <- PercentageFeatureSet(seu_obj, pattern = "^MT-")
        seu_obj <- subset(seu_obj, subset = nFeature_RNA > 100 & nFeature_RNA < 7000 & percent.mt < 10)
        seu_obj <- NormalizeData(seu_obj, normalization.method = "LogNormalize", scale.factor = 10000)
        return(seu_obj)
    }
    
    data_merge <- function(data_list){
        features <- SelectIntegrationFeatures(object.list=data_list, nfeatures=6000)
        data_anchors <- FindIntegrationAnchors(object.list=data_list, anchor.features=features)
        data_combined <- IntegrateData(anchorset=data_anchors)
        DefaultAssay(data_combined) <- "integrated"
        data_combined
        return(data_combined)
    }
    dir.create(res_dir)
    setwd(res_dir)
    ath_list <- list()

    for(i in 1:length(sample_dirs)){
        ath_list[[project_names[i]]] <- load_sc(data_di=sample_dirs[i], project_name=project_names[i])
    }

#    renew_list <- list()
#    for (i in 1:length(ath_list)) {
#        ath_renew <- pipe(ath_list[[i]], sample_name=project_names[i])
#        ath_list[[i]] <- FindVariableFeatures(ath_list[[i]], selection.method = "vst", nfeatures = 6000, verbose = F)
#        ath_renew <- FindVariableFeatures(ath_renew, selection.method = "vst", nfeatures = 6000, verbose = F)
#        renew_list[[i]] <- ath_renew
#        print(ath_renew)
#    }

    ath_combined <- data_merge(ath_list)
#    renew_combined <- data_merge(renew_list)
#    print(renew_combined)
    
    ath_combined <- data_process(ath_combined, res_dir="old", lnc=FALSE)
    saveRDS(ath_combined, file = "./R_result.rds")
#    setwd("../")

    print("no lnc done")

#    renew_combined <- data_process(renew_combined, res_dir="new")
#    saveRDS(renew_combined, file = "./R_result_renew.rds")
#    print("all pipe done")
}

#############################################***Function Test***#############################################
# data_00 <- '/home/li/anwser/ysw/work/scRNA/step2_cellranger/SRR10620013/outs/filtered_feature_bc_matrix'
# data_01 <- '/home/li/anwser/ysw/work/scRNA/step2_cellranger/SRR10620014/outs/filtered_feature_bc_matrix'
library(stringr)
opt <- getopt(spec=spec)

if( !is.null(opt$help) || is.null(opt$sample_dir)){
    cat(paste(getopt(spec=spec, usage = T), "\n"))
    quit()
}

print(opt)

all_dirs <- list.files(opt$sample_dir, full.names=T)
samples_dir <- all_dirs[grep("SRR", all_dirs)]
sample_names <- str_split_fixed(samples_dir, "//", 2)[,2]
sample_matrix_dirs <- paste0(samples_dir, "/outs/filtered_feature_bc_matrix")
res_dir <- opt$output_path
print(sample_matrix_dirs)
pipe_merge(sample_dirs=sample_matrix_dirs, project_names=sample_names, res_dir=res_dir)