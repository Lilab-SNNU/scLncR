#################################################***Part3 Single Cell Data PreProcess***##################################################
suppressPackageStartupMessages({
    suppressWarnings({
        suppressMessages({
            library(optparse)
            library(stringr)
            library(gridExtra) 
            library(psych)
            library(pheatmap)
            library(Seurat)
            library(dplyr)
            library(ggsci)
            library(tidyverse)
            library(patchwork)
            library(ggplot2)
            library(reshape2)
            library(scMayoMap)
            library(shiny)
            library(shinyjs)
            library(SingleR)
            library(ggrepel)
        })
    })
})

########################################*** Data load and Create Seu_obj***########################################
load_sc <- function(data_dir="", project_name="", mt_name="", min.RNAs=100, max.RNAs=7000, pct.mt=10, lnc_name="", nfeatures=6000){
    cat("=============== Data loading ...... \n")
    print(data_dir)
    data <- Read10X(data.dir = data_dir)
    seu_obj <- CreateSeuratObject(counts = data, project = project_name)
    
    if(mt_name!=""){
        seu_obj[["percent.mt"]] <- PercentageFeatureSet(seu_obj, pattern = paste0("^", mt_name))
    }
    
    saveRDS(seu_obj, file=paste0(project_name, "_seu_obj_raw.rds"))
    cat("=============== Data load and create seu_obj success. \n")
    
    seu_obj <- data_QC(seu_obj, min.RNAs=min.RNAs, max.RNAs=max.RNAs, pct.mt=pct.mt)
    seu_obj <- Inde_Normalize(seu_obj, lnc_name=lnc_name)
    seu_obj <- FindVariableFeatures(seu_obj, selection.method = "vst", nfeatures = nfeatures, verbose = F)
    
    return(seu_obj)
}

########################################*** Data QC***########################################

data_QC <- function(seu_obj, min.RNAs=100, max.RNAs=7000, pct.mt=10){
    cat("=============== Data QC ...... \n")
    # Visualize QC metrics as a violin plot
    plot <- VlnPlot(seu_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
    ggsave('FeatureQC.pdf', plot=plot, width = 8, height = 6)

    # FeatureScatter is typically used to visualize feature-feature relationships, but can be used
    # for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
    plot1 <- FeatureScatter(seu_obj, feature1 = "nCount_RNA", feature2 = "percent.mt")
    plot2 <- FeatureScatter(seu_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
    ggsave('FeatureScatter_1.pdf', plot=plot1, width = 8, height = 6)
    ggsave('FeatureScatter_2.pdf', plot=plot2, width = 8, height = 6)
    
    seu_obj <- subset(seu_obj, subset = nFeature_RNA > min.RNAs & nFeature_RNA < max.RNAs & percent.mt < pct.mt)
    cat("=============== Data QC Done \n")
    return(seu_obj)
}

########################################***Cut Data From All Matrix***########################################

cut_cell_matrix <- function(data_matrix, pattern){
    cat("=============== Split mRNA and lncRNA matrix...... \n")
    row_indexs <- grep(pattern, rownames(data_matrix))
    data_lnc <- data_matrix[row_indexs, ]
    data_gene <- data_matrix[-row_indexs, ]
    seu_lnc <- CreateSeuratObject(counts = data_lnc, project = "lncRNA")
    seu_gene <- CreateSeuratObject(counts = data_gene, project = "mRNA")
    seu_objs <- list(seu_lnc=seu_lnc, seu_gene=seu_gene)
    cat("=============== Split mRNA and lncRNA matrix success...... \n")
    return(seu_objs)
}


#############################################***Data Normalize ***#############################################

Inde_Normalize <- function(seu_obj, lnc_name=""){
    cat("=============== Independent normalization of mRNA and lncRNA matrix ...... \n")
    df_count <- seu_obj@assays$RNA@counts
    lnc_pattern <- sprintf("^%s.*", lnc_name)
    seu_objs <- cut_cell_matrix(data_matrix=df_count, pattern=lnc_pattern)
    seu_lnc <- seu_objs$seu_lnc
    seu_gene <- seu_objs$seu_gene
    
    seu_lnc <- NormalizeData(seu_lnc, normalization.method = "LogNormalize", scale.factor = 10000)
    seu_gene <- NormalizeData(seu_gene, normalization.method = "LogNormalize", scale.factor = 10000)
    
    df_bind <- rbind(seu_gene@assays$RNA@data, seu_lnc@assays$RNA@data)
    seu_bind <- seu_obj
    seu_bind@assays$RNA@data <- df_bind
    cat("=============== Independent normalization of mRNA and lncRNA matrix success \n")
    return(seu_bind)
}


#############################################***Data Merge ***#############################################

data_merge <- function(data_list, nfeatures=nfeatures){
    cat("=============== Multi samples Integrating ... \n")
    features <- SelectIntegrationFeatures(object.list=data_list, nfeatures=nfeatures)
    data_anchors <- FindIntegrationAnchors(object.list=data_list, anchor.features=features)
    data_combined <- IntegrateData(anchorset=data_anchors)
    DefaultAssay(data_combined) <- "integrated"
    cat("=============== Multi samples Integrate Done \n")
    return(data_combined)
}


#############################################*** Add Meta info to Seu_obj ***#############################################

data_group <- function(seu_obj, samples_info=""){
    cat("=============== Add samples meta information ... \n")
    meta_df <- read.table(samples_info, header = TRUE, sep = "\t")

    if ("orig.ident" %in% colnames(seu_obj@meta.data)) {
      # extract sample names
      samples <- as.character(seu_obj@meta.data$orig.ident)
      
      # map meta info by samples
      for (col_name in c("Group", "Deg_analysis", "site", "location_analysis")) {
        if (col_name %in% colnames(meta_df)) {
          sample_to_value <- setNames(meta_df[[col_name]], meta_df$sample)
          # add meta info to seu_obj
          seu_obj@meta.data[[col_name]] <- sample_to_value[samples]
        }
      }
    }
    cat("=============== Add samples meta information success \n")
    return(seu_obj)
}


#############################################***Data Process***#############################################

data_preprocess <- function(seu_obj, lnc=TRUE, lnc_name=lnc_name, nfeatures=nfeatures, dims=15, resolution=0.8, output_dir="", colour=colour){
    
    ####################################***PCA Process***####################################
    #Identification of highly variable features (feature selection)
    seu_obj <- FindVariableFeatures(seu_obj, selection.method = "vst", nfeatures = nfeatures)

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
    
    pdf('DimHeatmap.pdf', width = 10, height = 8)
    DimHeatmap(seu_obj, dims = 1:15, cells = 500, balanced = TRUE)
    dev.off()

    ####################################***Cluster Process***####################################
    # NOTE: This process can take a long time for big datasets, comment out for expediency. More
    # approximate techniques such as those implemented in ElbowPlot() can be used to reduce
    # computation time
    seu_obj <- JackStraw(seu_obj, num.replicate = 100)
    seu_obj <- ScoreJackStraw(seu_obj, dims = 1:dims)
    plot1<-JackStrawPlot(seu_obj, dims = 1:dims)
    plot2<-ElbowPlot(seu_obj)
    ggsave('JackStrawPlot.pdf', plot=plot1, width = 8, height = 6)
    ggsave('ElbowPlot.pdf', plot=plot2, width = 8, height = 6)


    # Cluster the cells 
    # Identify clusters of cells by a shared nearest neighbor (SNN) modularity optimization based clustering algorithm. 
    # First calculate k-nearest neighbors and construct the SNN graph. Then optimize the modularity function to determine clusters. 
    # For a full description of the algorithms, see Waltman and van Eck (2013) The European Physical Journal B. 
    seu_obj <- FindNeighbors(seu_obj, dims = 1:dims)
    seu_obj <- FindClusters(seu_obj, resolution = resolution)

    #Constructs a phylogenetic tree relating the 'average' cell from each identity class. 
    # Tree is estimated based on a distance matrix constructed in either gene expression space or PCA spac

    seu_obj<-BuildClusterTree(seu_obj)
    Tool(object = seu_obj, slot = 'BuildClusterTree')

    pdf('ClusterTreePlot.pdf', width = 8, height = 6)
    PlotClusterTree(seu_obj)
    dev.off()
    saveRDS(seu_obj, file = "./tmp.rds")
    
    # Run non-linear dimensional reduction (UMAP/tSNE)
      
    seu_obj <- RunUMAP(seu_obj, dims = 1:dims)

    seu_obj <- RunTSNE(seu_obj, dims = 1:dims)



    # note that you can set `label = TRUE` or use the LabelClusters function to help label
    # individual clusters
    plot1<-DimPlot(seu_obj, reduction = "umap",label = TRUE, raster=FALSE)+scale_color_manual(values = colour)
    plot2<-DimPlot(seu_obj, reduction = "tsne",label = TRUE, raster=FALSE)+scale_color_manual(values = colour)
    # CombinePlots(plots = list(plot1, plot2),legend="bottom")
    ggsave(file.path(output_dir, "umap.pdf"), plot=plot1, width = 8, height = 6, dpi=600)
    ggsave(file.path(output_dir, "tsne.pdf"), plot=plot2, width = 8, height = 6, dpi=600)
    ggsave(file.path(output_dir, "umap.png"), plot=plot1, width = 8, height = 6, dpi=600)
    ggsave(file.path(output_dir, "tsne.png"), plot=plot2, width = 8, height = 6, dpi=600)
    
    ####################################***Marker Gene Filter***####################################
    # Finding differentially expressed features (cluster biomarkers) 
    # find all markers of cluster 1

    # seu_obj.markers <- FindAllMarkers(seu_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

    if (lnc){
        lnc_pattern <- sprintf("^%s.*", lnc_name)
        lnc_trans <- grep(lnc_pattern, seu_obj.markers$gene)
        lncRNAs <- seu_obj.markers[lnc_trans, ]
        seu_obj@assays$sig_lncRNA <- lncRNAs
        write.csv(lncRNAs, file="markers_lncRNAs_matrix.csv")
        
        lncs <- rownames(lncRNAs)
        plot1<- FeaturePlot(seu_obj, features = lncs, min.cutoff = 0, max.cutoff = 4)
        plot2<- DoHeatmap(seu_obj, features = lncs) + NoLegend()+scale_color_manual(values = colour)
        ggsave("FeaturePlot_lncRNAs.pdf", plot=plot1, width = 20, height = 30, dpi=600)
        ggsave("DoHeatmap_lncRNAs.pdf", plot=plot2, width = 8, height = 6, dpi=600)
        ggsave("FeaturePlot_lncRNAs.png", plot=plot1, width = 20, height = 30, dpi=600)
        ggsave("DoHeatmap_lncRNAs.png", plot=plot2, width = 8, height = 6, dpi=600)
    }
    # top2 <- seu_obj.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
    
    return(seu_obj)
}

#############################################*** Data Auto Annotate***#############################################

sc_anno_SingleR <- function(seu_obj, ref_file="", output_dir="") {
    ref <- readRDS(ref_file)
    sce_ref <- Seurat::as.SingleCellExperiment(ref)
    norm_count <- GetAssayData(seu_obj, slot = "data")
    
    pred <- SingleR(test = norm_count, 
                    ref = sce_ref, 
                    labels = ref$Celltype, 
                    de.method = "wilcox")

    seu_obj$cell_type <- pred$labels
    saveRDS(seu_obj, file = file.path(output_dir, "singleR_anno.RDS"))
    return (seu_obj)
}

sc_anno_scMM <- function(seu_obj, db_path="", output_dir="", tissue = 'root', n_top_markers = 1) {

    need_info <- read.csv(db_path)
    need_info <- distinct(need_info)
    scMayoMapDatabase <- spread(need_info, key = c('celltype'), value = 'value')
    scMayoMapDatabase[is.na(scMayoMapDatabase)] <- 0

    # 找marker
    # seurat.markers <- FindAllMarkers(seu_obj, method = 'MAST')
    seurat.markers <- FindAllMarkers(seu_obj, only.pos = TRUE, logfc.threshold = 0.25)
    write.csv(seurat.markers, file="seu_obj_markers.csv")
    top_n_markers <- seurat.markers %>% group_by(cluster) %>% top_n(n = n_top_markers, wt = avg_log2FC)
  
    vln_p <- VlnPlot(seu_obj, features = unique(top_n_markers$gene), pt.size = 0, ncol = 1) +
    
                     scale_x_discrete("") +
                     theme(axis.text.x.bottom = element_blank())
  
    dot_p <- DotPlot(seu_obj, features = unique(top_n_markers$gene)) + RotatedAxis() +
                     scale_x_discrete("") + scale_y_discrete("")
  
    ggsave(file.path(output_dir, "VlnPlot.pdf"), plot = vln_p, width = 20, height = 70, limitsize=F, dpi=600)
    ggsave(file.path(output_dir, "Dotmap.pdf"), plot = dot_p, width = 60, height = 15, limitsize=F, dpi=600)
    ggsave(file.path(output_dir, "VlnPlot.png"), plot = vln_p, width = 20, height = 70, limitsize=F, dpi=600)
    ggsave(file.path(output_dir, "Dotmap.png"), plot = dot_p, width = 60, height = 15, limitsize=F, dpi=600)
  
    # cluster 相关性查看
    top_10_markers <- seurat.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
    AverageExp <- AverageExpression(seu_obj, features = unique(top_10_markers$gene))

    coorda <- corr.test(as.matrix(AverageExp$RNA), as.matrix(AverageExp$RNA), method = "spearman")
    plot_corr <- pheatmap(coorda$r)

    ggsave(file.path(output_dir, "pheatmap.pdf"), plot = plot_corr, width = 8, height = 6, , dpi=600)
    ggsave(file.path(output_dir, "pheatmap.png"), plot = plot_corr, width = 8, height = 6, , dpi=600)

    # cluster 注释
    scMayoMap.obj <- scMayoMap(data = seurat.markers, database = scMayoMapDatabase, tissue = tissue)
    plt <- scMayoMap.plot(scMayoMap.object = scMayoMap.obj)

    ggsave(plt, file = file.path(output_dir, "plt_scMM.pdf"), width = 16, height = 12, dpi = 600)
    ggsave(plt, file = file.path(output_dir, "plt_scMM.png"), width = 16, height = 12, dpi = 600)
    saveRDS(scMayoMap.obj, file = file.path(output_dir, "scMayoMap.RDS"))
  
    celltype <- scMayoMap.obj$markers %>%
                group_by(cluster) %>%
                slice_max(n = 1, order_by = score)
  
    result <- celltype %>%
               group_by(cluster) %>%
               filter(score == max(score)) %>%
               slice_tail(n = 1) %>%
               ungroup()
  
    result <- result[, c(1, 2)]
    result$cluster <- as.character(result$cluster)

    complete_annotation <- result
    seu_obj$seurat_clusters <- as.character(seu_obj$seurat_clusters)

    seu_obj$scMM <- complete_annotation$celltype[match(seu_obj$seurat_clusters, complete_annotation$cluster)]
    saveRDS(seu_obj, file = file.path(output_dir, "scMM_anno.RDS"))
    return(seu_obj)
}

Anno_plot <- function(seu_obj, output_dir="", colour=colour){
    
    plot1 <- DimPlot(seu_obj, reduction = "umap", label = TRUE) + scale_color_manual(values = colour)
    plot2 <- DimPlot(seu_obj, reduction = "tsne", label = TRUE) + scale_color_manual(values = colour)
    plot3 <- DimPlot(seu_obj, reduction = "umap", label = TRUE, group.by = "cell_type") + scale_color_manual(values = colour)
    plot4 <- DimPlot(seu_obj, reduction = "tsne", label = TRUE, group.by = "cell_type") + scale_color_manual(values = colour)

    ggsave(file.path(output_dir, "/umap.pdf"), plot = plot1, width = 8, height = 6, dpi = 600)
    ggsave(file.path(output_dir, "/tsne.pdf"), plot = plot2, width = 8, height = 6, dpi = 600)
    ggsave(file.path(output_dir, "/umap_anno.pdf"), plot = plot3, width = 8, height = 6, dpi = 600)
    ggsave(file.path(output_dir, "/tsne_anno.pdf"), plot = plot4, width = 8, height = 6, dpi = 600)
    ggsave(file.path(output_dir, "/umap.png"), plot = plot1, width = 8, height = 6, dpi = 600)
    ggsave(file.path(output_dir, "/tsne.png"), plot = plot2, width = 8, height = 6, dpi = 600)
    ggsave(file.path(output_dir, "/umap_anno.png"), plot = plot3, width = 8, height = 6, dpi = 600)
    ggsave(file.path(output_dir, "/tsne_anno.png"), plot = plot4, width = 8, height = 6, dpi = 600)
    
    if("Group" %in% colnames(seu_obj@meta.data)){
        split_p1 <- DimPlot(seu_obj, reduction = "umap", label = TRUE, split.by = "Group", raster=FALSE) + scale_color_manual(values = colour)
        split_p2 <- DimPlot(seu_obj, reduction = "umap", label = TRUE, group.by = "cell_type", split.by = "Group") + scale_color_manual(values = colour)

        ggsave(file.path(output_dir, "/umap_group.pdf"), plot = split_p1, width = 18, height = 6, dpi = 600)
        ggsave(file.path(output_dir, "/umap_anno_group.pdf"), plot = split_p2, width = 18, height = 6, dpi = 600)
        ggsave(file.path(output_dir, "/umap_group.png"), plot = split_p1, width = 18, height = 6, dpi = 600)
        ggsave(file.path(output_dir, "/umap_anno_group.png"), plot = split_p2, width = 18, height = 6, dpi = 600)
    }
    
    if("site" %in% colnames(seu_obj@meta.data)){
        split_p1 <- DimPlot(seu_obj, reduction = "umap", label = TRUE, split.by = "site", raster=FALSE) + scale_color_manual(values = colour)
        split_p2 <- DimPlot(seu_obj, reduction = "umap", label = TRUE, group.by = "cell_type", split.by = "site") + scale_color_manual(values = colour)

        ggsave(file.path(output_dir, "/umap_site.pdf"), plot = split_p1, width = 18, height = 6, dpi = 600)
        ggsave(file.path(output_dir, "/umap_anno_site.pdf"), plot = split_p2, width = 18, height = 6, dpi = 600)
        ggsave(file.path(output_dir, "/umap_site.png"), plot = split_p1, width = 18, height = 6, dpi = 600)
        ggsave(file.path(output_dir, "/umap_anno_site.png"), plot = split_p2, width = 18, height = 6, dpi = 600)
    }

}

#############################################*** Cell Type Statistics***#############################################
cell_type_stats <- function(seu_obj, colour=colour, output_dir=""){
    
    if("Group" %in% colnames(seu_obj@meta.data)){
        # 1. Extract and analyze data
        celltype_df <- seu_obj@meta.data %>%
          dplyr::select(Group, cell_type) %>%
          dplyr::count(Group, cell_type)

        # 2. Calculate percentages and label positions
        celltype_percent_df <- celltype_df %>%
          group_by(Group) %>%
          mutate(
            Percentage = n / sum(n) * 100,
            cumulative_percentage = cumsum(Percentage),
            label_position = cumulative_percentage - Percentage/2,
            label_text = paste0(round(Percentage, 1), "%")
          )
          
        # 3. Create a stacked bar chart with percentage labels
        stacked_bar_chart <- ggplot(celltype_percent_df, aes(x = Group, y = Percentage, fill = cell_type)) +
            geom_bar(stat = "identity", position = "stack") +
            geom_text(
                aes(label = ifelse(Percentage > 5, paste0(round(Percentage, 1), "%"), "")),
                position = position_stack(vjust = 0.5),
                size = 3,
                color = "white",
                fontface = "bold"
            ) +
            labs(title = "Cell Type Composition by Group (Stacked Bar Chart)",
                x = "Group",
                y = "Percentage of Cells (%)",
                fill = "Cell Type") +
            theme_classic() +
            theme(axis.text.x = element_text(angle = 45, hjust = 1))
            
        ggsave(file.path(output_dir, "/celltype_stacked_bar_with_Group.pdf"), 
          plot = stacked_bar_chart,
          width = 10, height = 8, units = "in"
        )
    }
    
    if("site" %in% colnames(seu_obj@meta.data)){
        # 1. Extract and analyze data
        celltype_df <- seu_obj@meta.data %>%
            dplyr::select(site, cell_type) %>%
            dplyr::count(site, cell_type)

        # 2. Calculate percentages and label positions
        celltype_percent_df <- celltype_df %>%
            group_by(site) %>%
            mutate(
                Percentage = n / sum(n) * 100,
                cumulative_percentage = cumsum(Percentage),
                label_position = cumulative_percentage - Percentage/2,
                label_text = paste0(round(Percentage, 1), "%")
            )
          
        # 3. Create a stacked bar chart with percentage labels
        stacked_bar_chart <- ggplot(celltype_percent_df, aes(x = site, y = Percentage, fill = cell_type)) +
            geom_bar(stat = "identity", position = "stack") +
            geom_text(
                aes(label = ifelse(Percentage > 5, paste0(round(Percentage, 1), "%"), "")),
                position = position_stack(vjust = 0.5),
                size = 3,
                color = "white",
                fontface = "bold"
            ) +
            labs(title = "Cell Type Composition by site (Stacked Bar Chart)",
                x = "site",
                y = "Percentage of Cells (%)",
                fill = "Cell Type") +
            theme_classic() +
            theme(axis.text.x = element_text(angle = 45, hjust = 1))
            
        ggsave(file.path(output_dir, "/celltype_stacked_bar_with_site.pdf"), 
            plot = stacked_bar_chart,
            width = 10, height = 8, units = "in"
        )
    }
}


#############################################*** Data Process Pipeline ***#############################################

run_dataProcess <- function(user_args, script_dir) {
    # Load required package
    if (!requireNamespace("optparse", quietly = TRUE)) {
    stop("Package 'optparse' is required. Please install it: install.packages('optparse')")
    }

    # Define option list
    option_list <- list(
        optparse::make_option(
            c("-c", "--config"),
            type = "character",
            help = "Path to the YAML configuration file (e.g., config.yaml)",
            metavar = "FILE"
        )
    )

    # Create parser
    parser <- optparse::OptionParser(
        usage = "scLncR dataProcess [options]",
        option_list = option_list,
        description = "dataProcess: Integrate, preprocess, and annotate single-cell RNA-seq data"
    )

    # Parse arguments
    opt <- optparse::parse_args(parser, args = user_args, print_help_and_exit = FALSE)

    # If help requested or config missing
    if (is.null(opt$config)) {
        cat("\n")
        optparse::print_help(parser)
        q()
    }

    config_file <- opt$config
    if (!file.exists(config_file)) {
        stop("Config file not found: ", config_file)
    }

    # Load and validate config
    cfg <- load_dataProcess_config(config_file)

    # Extract parameters (optional: you can just pass cfg directly to functions)
    counts_dir <- cfg$counts_dir
    mt_name  <- cfg$mt_name
    min.RNAs    <- cfg$min.RNAs
    max.RNAs <- cfg$max.RNAs
    pct.mt  <- cfg$percent.mt
    lnc_name      <- cfg$lnc_name
    nfeatures  <- cfg$nfeatures
    dims <- cfg$dims
    resolution <- cfg$resolution
    samples_info <- cfg$samples_info
    anno_method <- cfg$anno_method
    ref_file <- cfg$ref_file
    marker_gene_file <- cfg$marker_gene_file
    tissue <- cfg$tissue
    n_top_markers <- cfg$n_top_markers
    output_dir <- cfg$output_dir

    # Build local parameters 
    rawdata_dir <- paste0(output_dir, "/sample_seurat_raw")
    preprocess_dir <- paste0(output_dir, "/data_preprocess")
    anno_dir <- paste0(output_dir, "/data_annotation")
    dir.create(output_dir, recursive = TRUE)
    dir.create(rawdata_dir, recursive = TRUE)
    dir.create(preprocess_dir, recursive = TRUE)
    dir.create(anno_dir, recursive = TRUE)
    setwd(output_dir)
    
    count_dirs <- system(sprintf("find %s -name filtered_feature_bc_matrix", counts_dir), intern=T)
    project_names <- basename(dirname(dirname(count_dirs)))
    seu_list <- list()
    colour=c("#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB","#98FB98","#F08080","#1E90FF","#7CFC00","#FFFF00",  
         "#808000","#FF00FF","#FA8072","#7B68EE","#9400D3","#800080","#A0522D","#D2B48C","#D2691E","#87CEEB","#40E0D0","#5F9EA0",
         "#FF1493","#0000CD","#008B8B","#FFE4B5","#8A2BE2","#228B22","#E9967A","#4682B4","#32CD32","#F0E68C","#FFFFE0","#EE82EE",
         "#FF6347","#6A5ACD","#9932CC","#8B008B","#8B4513","#DEB887")
    
    for(i in 1:length(count_dirs)){
        print(count_dirs[i])
        print(project_names[i])
        
        seu_list[[project_names[i]]] <- load_sc(data_dir=count_dirs[i], 
                                                project_name=project_names[i], 
                                                mt_name=mt_name, 
                                                min.RNAs=min.RNAs, 
                                                max.RNAs=max.RNAs, 
                                                pct.mt=pct.mt, 
                                                lnc_name=lnc_name, 
                                                nfeatures = nfeatures)
    }
    
    system(paste0("mv *.rds ", rawdata_dir))
    
    seu_obj <- data_merge(seu_list, nfeatures=nfeatures)
    saveRDS(seu_obj, file=paste0(rawdata_dir, "/seu_obj_inte_raw.rds"))
    
    if(samples_info != ""){
        seu_obj <- data_group(seu_obj, samples_info=samples_info)
    }
    
    seu_obj <- data_preprocess(seu_obj, lnc=F, lnc_name=lnc_name, nfeatures=nfeatures, dims=dims, resolution=resolution, output_dir=output_dir, colour=colour)
    saveRDS(seu_obj, file=file.path(preprocess_dir, "/preprocessed_result.rds"))
    
    if(anno_method == "SingleR"){
        seu_obj <- sc_anno_SingleR(seu_obj, ref_file=ref_file, output_dir=anno_dir)
    }
    
    if(anno_method == "scMM"){
        seu_obj <- sc_anno_scMM(seu_obj, db_path=marker_gene_file, output_dir=anno_dir, tissue=tissue, n_top_markers=n_top_markers)
    }
    
    Anno_plot(seu_obj, output_dir=anno_dir, colour=colour)
    cell_type_stats(seu_obj, colour=colour, output_dir=anno_dir)
    
    saveRDS(seu_obj, file=file.path(anno_dir, "/anno_result.rds"))
}
