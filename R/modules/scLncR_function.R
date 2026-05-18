########################################***snRNA/scRNA enrichment analysis***########################################

SnScEnrichmentAnalysis <- function(seu_obj, LOG2FC_THRESH=0.25, PADJ_THRESH=0.05, output_path="./", lncRNA_name="scLncR"){
    dir.create(output_path, recursive = TRUE, showWarnings = FALSE)
    setwd(output_path)
    if (!("cell_type" %in% colnames(seu_obj@meta.data))) {
      seu_obj$cell_type <- as.character(seu_obj$seurat_clusters)
    }

    cell_types <- sort(unique(seu_obj$cell_type))

    sc_enriched_genes_list <- list()
    sn_enriched_genes_list <- list()
    balanced_genes_list <- list()
    deg_full_list <- list()

    for (ct in cell_types) {
      message("Processing cell type: ", ct)

      subset_obj <- subset(seu_obj, subset = cell_type == ct)
      Idents(subset_obj) <- subset_obj$site

      group_table <- table(Idents(subset_obj))
      if (length(group_table) < 2 || any(group_table == 0)) {
        message("  -> Skipping: only one group present.")
        sc_enriched_genes_list[[ct]] <- character(0)
        sn_enriched_genes_list[[ct]] <- character(0)
        balanced_genes_list[[ct]] <- character(0)
        next
      }

      deg <- FindMarkers(
        object = subset_obj,
        ident.1 = "scRNA",
        ident.2 = "snRNA",
        test.use = "wilcox",
        logfc.threshold = 0,
        min.pct = 0.1,
        return.thresh = 1
      )

      if (nrow(deg) == 0) {
        sc_enriched_genes_list[[ct]] <- character(0)
        sn_enriched_genes_list[[ct]] <- character(0)
        balanced_genes_list[[ct]] <- character(0)
        deg_full_list[[ct]] <- deg
        next
      }

      deg_df <- as.data.frame(deg)
      if ("gene" %in% colnames(deg_df)) deg_df$gene <- NULL
      deg_df <- deg_df %>%
        rownames_to_column("gene") %>%
        mutate(
          category = case_when(
            avg_log2FC > LOG2FC_THRESH & p_val_adj < PADJ_THRESH ~ "scRNA-enriched",
            avg_log2FC < -LOG2FC_THRESH & p_val_adj < PADJ_THRESH ~ "snRNA-enriched",
            TRUE ~ "Balanced/Non-differential"
          )
        )

      sc_enriched_genes_list[[ct]] <- deg_df$gene[deg_df$category == "scRNA-enriched"]
      sn_enriched_genes_list[[ct]] <- deg_df$gene[deg_df$category == "snRNA-enriched"]
      balanced_genes_list[[ct]] <- deg_df$gene[deg_df$category == "Balanced/Non-differential"]

      deg_full_list[[ct]] <- deg_df

      message("  -> scRNA-enriched:", length(sc_enriched_genes_list[[ct]]),
              "; snRNA-enriched:", length(sn_enriched_genes_list[[ct]]),
              "; Balanced:", length(balanced_genes_list[[ct]]))
    }

    summary_list <- tibble(
      cell_type = character(),
      category = character(),
      count = numeric(),
      proportion = numeric()
    )

    for (ct in names(sc_enriched_genes_list)) {
      total <- length(sc_enriched_genes_list[[ct]]) + length(sn_enriched_genes_list[[ct]]) + length(balanced_genes_list[[ct]])
      if (total == 0) next

      df_temp <- tibble(
        cell_type = ct,
        category = c("scRNA-enriched", "snRNA-enriched", "Balanced/Non-differential"),
        count = c(length(sc_enriched_genes_list[[ct]]), length(sn_enriched_genes_list[[ct]]), length(balanced_genes_list[[ct]]))
      ) %>%
        mutate(proportion = count / total)

      summary_list <- bind_rows(summary_list, df_temp)
    }

    dir.create("sn_vs_sc_enrichment_results", showWarnings = FALSE, recursive = TRUE)
    saveRDS(list(
      sc = sc_enriched_genes_list,
      sn = sn_enriched_genes_list,
      balanced = balanced_genes_list,
      deg_full = deg_full_list
    ), file = "sn_vs_sc_enrichment_results/gene_classification.rds")

    write.csv(summary_list, "sn_vs_sc_enrichment_results/gene_category_summary.csv", row.names = FALSE)
    write.csv(deg_full_list, "sn_vs_sc_enrichment_results.csv")

    sn_sc_summary <- tibble(
      cell_type = character(),
      category = character(),
      count = numeric()
    )

    for (ct in names(sc_enriched_genes_list)) {
      n_sc <- length(sc_enriched_genes_list[[ct]])
      n_sn <- length(sn_enriched_genes_list[[ct]])

      df_temp <- tibble(
        cell_type = ct,
        category = c("scRNA-enriched", "snRNA-enriched"),
        count = c(n_sc, n_sn)
      )
      sn_sc_summary <- bind_rows(sn_sc_summary, df_temp)
    }

    sn_sc_summary <- sn_sc_summary %>%
      group_by(cell_type) %>%
      mutate(percentage = count / sum(count) * 100) %>%
      ungroup()

    sn_sc_summary$category <- factor(
      sn_sc_summary$category,
      levels = c("snRNA-enriched", "scRNA-enriched")
    )

    p_stack_pct <- ggplot(sn_sc_summary, aes(x = cell_type, y = percentage, fill = category)) +
      geom_col(position = "fill") +
      scale_y_continuous(labels = scales::percent_format()) +
      scale_fill_manual(values = c(
        "snRNA-enriched" = "#E41A1C",
        "scRNA-enriched" = "#377EB8"
      )) +
      labs(
        title = "snRNA/scRNA Enriched Genes by Cell Type",
        x = "Cell Type",
        y = "Proportion",
        fill = "Expression Enrichment Category"
      ) +
      theme_minimal(base_size = 10) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))

    ggsave("sn_vs_sc_enrichment_results/stacked_pct_sn_vs_sc_enrichment.png", p_stack_pct, width = 8, height = 5, dpi=600)
    ggsave("sn_vs_sc_enrichment_results/stacked_pct_sn_vs_sc_enrichment.pdf", p_stack_pct, width = 8, height = 5, dpi=600)

    all_sc_lnc <- unlist(sc_enriched_genes_list)[startsWith(unlist(sc_enriched_genes_list), lncRNA_name)]
    all_sn_lnc <- unlist(sn_enriched_genes_list)[startsWith(unlist(sn_enriched_genes_list), lncRNA_name)]
    all_lnc_genes <- unique(c(all_sc_lnc, all_sn_lnc))
    lnc_freq_df <- tibble(gene = character(), frequency = integer())

    if (length(all_lnc_genes) == 0) {
      message("Warning: No lncRNA found in enrichment gene lists.")
    } else {
      lnc_freq <- integer(length(all_lnc_genes))
      names(lnc_freq) <- all_lnc_genes

      for (gene in all_lnc_genes) {
        count_cts <- 0
        for (ct in names(sc_enriched_genes_list)) {
          if (gene %in% sc_enriched_genes_list[[ct]] || gene %in% sn_enriched_genes_list[[ct]]) {
            count_cts <- count_cts + 1
          }
        }
        lnc_freq[gene] <- count_cts
      }

      lnc_freq_df <- data.frame(
        gene = names(lnc_freq),
        frequency = as.integer(lnc_freq)
      ) %>%
        arrange(desc(frequency))

      write.csv(lnc_freq_df, "sn_vs_sc_enrichment_results/lncRNA_frequency.csv", row.names = FALSE)

      top_lnc <- head(lnc_freq_df, 20)
      p_lnc <- ggplot(top_lnc, aes(x = reorder(gene, frequency), y = frequency)) +
        geom_col(fill = "purple") +
        coord_flip() +
        labs(
          title = "Top 20 lncRNAs by Number of Cell Types Detected",
          x = "lncRNA",
          y = "Number of Cell Types"
        ) +
        theme_minimal()

      ggsave("sn_vs_sc_enrichment_results/top_lncRNAs_frequency.png", p_lnc, width = 6, height = 5, dpi=600)
      ggsave("sn_vs_sc_enrichment_results/top_lncRNAs_frequency.pdf", p_lnc, width = 6, height = 5, dpi=600)

      message("✅ lncRNA enrichment analysis complete. Top lncRNAs plot saved.")
    }

    lnc_sn_sc_summary <- tibble(
      cell_type = character(),
      category = character(),
      count = numeric()
    )

    for (ct in names(sc_enriched_genes_list)) {
      n_sc_lnc <- length(sc_enriched_genes_list[[ct]][startsWith(sc_enriched_genes_list[[ct]], lncRNA_name)])
      n_sn_lnc <- length(sn_enriched_genes_list[[ct]][startsWith(sn_enriched_genes_list[[ct]], lncRNA_name)])

      if (n_sc_lnc + n_sn_lnc == 0) next

      df_temp <- tibble(
        cell_type = ct,
        category = c("scRNA-enriched", "snRNA-enriched"),
        count = c(n_sc_lnc, n_sn_lnc)
      )
      lnc_sn_sc_summary <- bind_rows(lnc_sn_sc_summary, df_temp)
    }

    if (nrow(lnc_sn_sc_summary) == 0) {
      message("Warning: No lncRNA found in any cell type for lncRNA-specific enrichment plot.")
    } else {
      lnc_sn_sc_summary <- lnc_sn_sc_summary %>%
        group_by(cell_type) %>%
        mutate(percentage = count / sum(count) * 100) %>%
        ungroup()

      lnc_sn_sc_summary$category <- factor(
        lnc_sn_sc_summary$category,
        levels = c("snRNA-enriched", "scRNA-enriched")
      )

      p_lnc_stack_pct <- ggplot(lnc_sn_sc_summary, aes(x = cell_type, y = percentage, fill = category)) +
        geom_col(position = "fill") +
        scale_y_continuous(labels = scales::percent_format()) +
        scale_fill_manual(values = c(
          "snRNA-enriched" = "#E41A1C",
          "scRNA-enriched" = "#377EB8"
        )) +
        labs(
          title = paste("snRNA/scRNA Enriched", lncRNA_name, "by Cell Type"),
          x = "Cell Type",
          y = paste("Proportion of", lncRNA_name, "Genes"),
          fill = "Expression Enrichment Category"
        ) +
        theme_minimal(base_size = 10) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))

      ggsave("sn_vs_sc_enrichment_results/stacked_pct_lncRNA_sn_vs_sc_enrichment.png", p_lnc_stack_pct, width = 8, height = 5, dpi=600)
      ggsave("sn_vs_sc_enrichment_results/stacked_pct_lncRNA_sn_vs_sc_enrichment.pdf", p_lnc_stack_pct, width = 8, height = 5, dpi=600)

      message("✅ lncRNA snRNA/scRNA enrichment plot saved.")
    }

    if (nrow(lnc_freq_df) == 0) {
      message("Warning: No lncRNA frequency table available; skipping lncRNA snRNA/scRNA enrichment heatmap.")
      return(invisible(NULL))
    }

    top_lnc_genes <- head(lnc_freq_df, 20)$gene
    lnc_matrix <- matrix(NA, nrow = length(top_lnc_genes), ncol = length(cell_types))
    rownames(lnc_matrix) <- top_lnc_genes
    colnames(lnc_matrix) <- cell_types

    for (gene in top_lnc_genes) {
      for (ct in cell_types) {
        if (gene %in% sc_enriched_genes_list[[ct]]) {
          lnc_matrix[gene, ct] <- "scRNA-enriched"
        } else if (gene %in% sn_enriched_genes_list[[ct]]) {
          lnc_matrix[gene, ct] <- "snRNA-enriched"
        }
      }
    }

    lnc_matrix <- lnc_matrix[rowSums(!is.na(lnc_matrix)) > 0, , drop = FALSE]
    if (nrow(lnc_matrix) == 0) {
      message("Warning: No lncRNA has snRNA/scRNA enrichment labels across cell types; skipping heatmap.")
      return(invisible(NULL))
    }
    lnc_heatmap_df <- reshape2::melt(lnc_matrix)
    colnames(lnc_heatmap_df) <- c("lncRNA", "cell_type", "enrichment_category")

    p_lnc_heatmap <- ggplot(lnc_heatmap_df, aes(x = cell_type, y = lncRNA, fill = enrichment_category)) +
      geom_tile() +
      scale_fill_manual(
        values = c(
          "snRNA-enriched" = "#D2691E",
          "scRNA-enriched" = "#32CD32"
        ),
        breaks = c("snRNA-enriched", "scRNA-enriched"),
        na.value = "#87CEEB"
      ) +
      labs(
        title = paste("snRNA/scRNA Expression Enrichment of Top 20 lncRNAs (", lncRNA_name, ") across Cell Types", sep = ""),
        x = "Cell Type",
        y = "lncRNA",
        fill = "Expression Enrichment Category"
      ) +
      theme_minimal(base_size = 10) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        axis.text.y = element_text(size = 7),
        panel.grid.major = element_line(color = "white"),
        legend.position = "right"
      )

    ggsave("sn_vs_sc_enrichment_results/lncRNA_sn_sc_enrichment_heatmap.pdf", p_lnc_heatmap, width = 10, height = 6, dpi=600)
    ggsave("sn_vs_sc_enrichment_results/lncRNA_sn_sc_enrichment_heatmap.png", p_lnc_heatmap, width = 10, height = 6, dpi=600)

    message("✅ lncRNA snRNA/scRNA enrichment heatmap saved (Top 20 lncRNAs) - NA cells are light blue")
}

LocationAnalysis <- function(...) {
    warning(
      "LocationAnalysis() is deprecated. Use SnScEnrichmentAnalysis() for snRNA/scRNA enrichment analysis. ",
      "The result reflects relative expression enrichment between snRNA-seq and scRNA-seq, not direct subcellular localization."
    )
    SnScEnrichmentAnalysis(...)
}

########################################***Trajectory analysis***########################################

# 定义主函数
monocle2_analysis <- function(seu_obj, qval=0.05, reduceDimensionMethod="DDRTree", output_path="", hub_genes="all"){
    
    # seu_obj <- subset(seu_obj, subset = cell_type == cell_types)
    #Extract data, phenotype data, and feature data from the SeuratObject
    data <- as.matrix(seu_obj@assays$RNA@data)
    pd <- new('AnnotatedDataFrame', data = seu_obj@meta.data)
    fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
    fd <- new('AnnotatedDataFrame', data = fData)

    #构建S4对象，CellDataSet
    seu_cds <- newCellDataSet(data,
                            phenoData = pd,
                            featureData = fd,
                            expressionFamily = negbinomial.size())
    cat("=========### Estimate the size factor for each cell \n")
    seu_cds <- estimateSizeFactors(seu_cds) 
    
    cat("=========### Estimating the dispersion of gene expression \n")
    seu_cds <- estimateDispersions(seu_cds)
    
    cat("=========### Differential gene analysis \n")
    # seu_cds <- detectGenes(seu_cds, min_expr = 0.1 )

    # expressed_genes <- row.names(subset(fData(seu_cds),
    #                                     num_cells_expressed >= 5))
    # diff_genes <- VariableFeatures(seu_obj)
    
    diff_test_res <- differentialGeneTest(seu_cds[hub_genes, ],
                                          fullModelFormulaStr = "~ cell_type + Group", 
                                          cores=24)
    
    cat("=========### Choose Differential gene as ordering_genes")
    ordering_genes <- row.names (subset(diff_test_res, qval < qval)) ## 不要也写0.1 ，而是要写0.01。
    seu_cds <- setOrderingFilter(seu_cds, ordering_genes)
    
    # disp_table <- dispersionTable(seu_cds)
    # ordering_genes_temp <- subset(disp_table, mean_expression >= 0.1) 
    # ordering_genes<-ordering_genes_temp$gene_id
    # seu_cds <- setOrderingFilter(seu_cds, ordering_genes)


    p1 <- plot_ordering_genes(seu_cds)
    ggsave(p1, file=paste0(output_path, "/plot_ordering_genes.pdf"), dpi=600)
    ggsave(p1, file=paste0(output_path, "/plot_ordering_genes.png"), dpi=600)

    seu_cds <- reduceDimension(seu_cds, max_components = 2, method = reduceDimensionMethod) # DDRTree方式
    seu_cds <- orderCells(seu_cds)


    colour=c("#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB","#98FB98","#F08080","#1E90FF","#7CFC00","#FFFF00",  
            "#808000","#FF00FF","#FA8072","#7B68EE","#9400D3","#800080","#A0522D","#D2B48C","#D2691E","#87CEEB","#40E0D0","#5F9EA0",
            "#FF1493","#0000CD","#008B8B","#FFE4B5","#8A2BE2","#228B22","#E9967A","#4682B4","#32CD32","#F0E68C","#FFFFE0","#EE82EE",
            "#FF6347","#6A5ACD","#9932CC","#8B008B","#8B4513","#DEB887")

    a1 <- plot_cell_trajectory(seu_cds, color_by = "seurat_clusters") + scale_color_manual(values = colour)

    a2 <- plot_cell_trajectory(seu_cds, color_by = "State") + scale_color_manual(values = colour)

    a3 <- plot_cell_trajectory(seu_cds, color_by = "Pseudotime") 

    a4 <- plot_cell_trajectory(seu_cds, color_by = "cell_type")  + scale_color_manual(values = colour)

    p2 <- (a1 + a2) / (a3 + a4)
    ggsave(p2, file=paste0(output_path, "/Pseudotime_test.pdf"), width=12, height=10, dpi=600)


    p3 <- plot_cell_trajectory(seu_cds, color_by = "cell_type") +
    facet_wrap(~cell_type, nrow = 2) #设置几行几列展示
    ggsave(p3, file=paste0(output_path, "/cell_trajectory.pdf"), dpi=600)
    
    b1 <- plot_cell_trajectory(seu_cds, x = 1, y = 2, color_by = "cell_type") + 
    theme(legend.position='none',panel.border = element_blank()) + #去掉第一个的legend
    scale_color_manual(values = colour) 
    b2 <- plot_complex_cell_trajectory(seu_cds, x = 1, y = 2, color_by = "cell_type")+
                                       scale_color_manual(values = colour) +
                                       theme(legend.title = element_blank()) 

    p4 <- b1 / b2
    ggsave(p4, file=paste0(output_path, "/complex_cell_trajectory.pdf"), dpi=600)
    saveRDS(seu_cds, file=paste0(output_path, "/seu_cds.rds"))

    plt_lnc <- 
    ppp <- plot_genes_jitter(seu_cds[plt_lnc,],
                  grouping = "Group",
                  color_by = "cell_type",
                  nrow= 2,
                  ncol = NULL,
                  plot_trend = TRUE)
    return(seu_cds)
}
########################################***WGCNA Analysis***########################################

#######################**********Set Seurat object for hdWGCNA
DO_WGCNA <- function(seu_obj, output_path="", cell_types=c(), pro_name="scLncR", lnc_name="AthLnc", gene_select_method="fraction"){
    dir.create(output_path, recursive = TRUE)
    DefaultAssay(object = seu_obj) <- "RNA"
    
    if(gene_select_method=="fraction"){
        seu_obj <- SetupForWGCNA(
            seu_obj,
            gene_select = "fraction", # the gene selection approach
            fraction = 0.01, 
            wgcna_name = pro_name # the name of the hdWGCNA experiment
        )
    }
    
    if(gene_select_method=="variable"){
        seu_obj <- FindVariableFeatures(seu_obj, selection.method = "vst", nfeatures = 6000)
        seu_obj <- SetupForWGCNA(
            seu_obj,
            gene_select = "variable", # the gene selection approach
            wgcna_name = pro_name # the name of the hdWGCNA experiment
        )
    }


    #######################********** Building Metacells

    # Building Metacells for each group info
    seu_obj <- MetacellsByGroups(
        seurat_obj = seu_obj,
        group.by = c("cell_type", "Group"), 
        k = 25, 
        max_shared = 10, 
        ident.group = 'cell_type'
    )
    
    # normalize metacell expression matrix:
    wgcna_name <- seu_obj@misc$active_wgcna
    metacell_obj <- GetMetacellObject(seu_obj, wgcna_name)
    metacell_obj <- Inde_Normalize(metacell_obj, lnc_name=lnc_name)
    seu_obj@misc[[wgcna_name]]$wgcna_metacell_obj <- metacell_obj

    #######################********** Co-expression network analysis

    ##########**********Set up the expression matrix
    seu_obj <- SetDatExpr(
        seu_obj,
        group_name = cell_types, # the name of the group of interest in the group.by column
        group.by='cell_type', # the metadata column containing the cell type info. This same column should have also been used in MetacellsByGroups
        assay = 'RNA', # using RNA assay
        slot = 'data' # using normalized data
    )
    ##########**********Select soft-power threshold
    # Test different soft powers:
    seu_obj <- TestSoftPowers(
        seu_obj,
        networkType = 'unsigned' # you can also use "signed" or "signed hybrid"
    )
    
    # plot the results:
    plot_list <- PlotSoftPowers(seu_obj)
    
    # assemble with patchwork
    p1 <- wrap_plots(plot_list, ncol=2)
    ggsave(p1, file=paste0(output_path,"/wrap_plots.pdf"), height=8, width=10, dpi=450)
    # check data
    power_table <- GetPowerTable(seu_obj)

    ##########**********Construct co-expression network
    # If no soft threshold is specified, ConstrucNetwork will automatically specify a soft threshold
    # construct co-expression network:
    tom_outdir <- file.path(output_path, "TOM")
    seu_obj <- ConstructNetwork(
        seu_obj,
        tom_name = pro_name, # name of the topoligical overlap matrix written to disk
        tom_outdir = tom_outdir, 
        overwrite_tom = TRUE 
    )
    seu_obj@misc[[wgcna_name]]$wgcna_net$TOMFiles <- tom_outdir

    pdf(paste0(output_path,"/PlotDendrogram.pdf"), height=8, width=10)
    PlotDendrogram(seu_obj, main=paste0(pro_name, ' hdWGCNA Dendrogram'))
    dev.off()

    ##########**********Compute harmonized module eigengenes
    # need to run ScaleData first or else harmony throws an error:
    #seu_obj <- ScaleData(seu_obj, features=VariableFeatures(seu_obj))
    
    seu_obj <- ScaleData(seu_obj, features=VariableFeatures(seu_obj))
    # compute all MEs in the full single-cell dataset
    seu_obj <- ModuleEigengenes(
        seu_obj,
        group.by.vars="Group"
    )
    
    # harmonized module eigengenes:
    # allow the user to apply Harmony batch correction to the MEs, yielding harmonized module eigengenes (hMEs)
    hMEs <- GetMEs(seu_obj)
    
    # module eigengenes:
    #MEs <- GetMEs(seu_obj, harmonized=FALSE)

    ##########**********Compute module connectivity
    # compute eigengene-based connectivity (kME):
    # focus on the “hub genes”
    seu_obj <- ModuleConnectivity(
        seu_obj,
        group.by = 'cell_type', 
        group_name = cell_types
    )
    
    # rename the modules
    seu_obj <- ResetModuleNames(
        seu_obj,
        new_name = paste0(pro_name, "_Module")
    )
    
    # plot genes ranked by kME for each module
    p <- PlotKMEs(seu_obj, ncol=5)
    p
    ggsave(p, file=paste0(output_path, "/PCgenes.pdf") ,width = 24,height =8)

    ##########**********Getting the module assignment table
    # get the module assignment table:
    modules <- GetModules(seu_obj) %>% 
    subset(module != 'grey')
    write.csv(modules, file=paste0(output_path, "/wgcna_modules.csv"))
    
    mods <- levels(modules$module)
    mods <- mods[mods != "grey"]
    
    hub_df <- do.call(rbind, lapply(mods, function(cur_mod) {
        cur <- subset(modules, module == cur_mod)
        cur <- cur[, c("gene_name", "module", paste0("kME_", 
            cur_mod))]
        names(cur)[3] <- "kME"
        cur <- dplyr::arrange(cur, desc(kME))
    }))
    write.csv(hub_df, file=paste0(output_path, "/wgcna_hubGenes.csv"))


    ##########**********Compute hub gene signature scores
    # compute gene scoring for the top 25 hub genes by kME for each module
    # with UCell method
    library(UCell)
    seu_obj <- ModuleExprScore(
        seu_obj,
        n_genes = 25,
        method='UCell' # Seurat方法(AddModuleScore)
    )

    #######################**********Basic Visualization
    # Create hMEs feature maps for each module
    plot_list <- ModuleFeaturePlot(
        seu_obj,
        features='hMEs', # plot the hMEs
        order=TRUE # order so the points with highest hMEs are on top
    )
    
    # stitch together with patchwork
    p3 <- wrap_plots(plot_list, ncol=5)
    ggsave(p3, file=paste0(output_path, "/combinePlot.pdf"),width = 20,height = 4)

    # The situation of each module in different cell subpopulations
    seu_obj$cluster <- do.call(rbind, strsplit(as.character(seu_obj$Group), ' '))[,1]
    
    p4 <- ModuleRadarPlot(
        seu_obj,
        group.by = 'cluster',
        barcodes = seu_obj@meta.data %>% 
            subset(cell_type == cell_types) %>% 
            rownames(),
        axis.label.size=4,
        grid.label.size=4
    )
    ggsave(p4, file=paste0(output_path, "/radarPlot_group.pdf") ,width = 24,height = 18)
    
    # View module related diagrams
    pdf(paste0(output_path, "/ModuleCorrelogram.pdf"))
    ModuleCorrelogram(seu_obj)
    dev.off()

    # Visualize the top 50 (customizable values) hub genes for each module using Module Network Plot
    # get hMEs from seurat object
    MEs <- GetMEs(seu_obj, harmonized=TRUE)
    modules <- GetModules(seu_obj)
    mods <- levels(modules$module); mods <- mods[mods != 'grey']

    # add hMEs to Seurat meta-data:
    seu_obj@meta.data <- cbind(seu_obj@meta.data, MEs)
    p <- DotPlot(seu_obj, features=mods, group.by = 'cell_type')

    # flip the x/y axes, rotate the axis labels, and change color scheme:
    p <- p +
      RotatedAxis() +
      scale_color_gradient2(high='red', mid='grey95', low='blue')
    ggsave(p, file=file.path(output_path,"dot_modules.pdf"), dpi=450, width=10)
    saveRDS(seu_obj, file=paste0(output_path, "/seu_wgcna.RDS"))
    return(seu_obj)
}

# =============================================================================
# scLncR - Main Runner for LncExplore Module
# =============================================================================

run_function <- function(args, script_dir) {
  # Load required packages
  if (!requireNamespace("optparse", quietly = TRUE)) {
    stop("Package 'optparse' is required. Install with: install.packages('optparse')")
  }
  
  # -----------------------------
  # 1. Define CLI options
  # -----------------------------
  option_list <- list(
    optparse::make_option(
      c("-c", "--config"),
      type = "character",
      help = "Path to LncExplore YAML configuration file",
      metavar = "FILE"
    )
  )
  
  parser <- optparse::OptionParser(
    usage = "[options]",
    option_list = option_list,
    description = "function explore: snRNA/scRNA enrichment, trajectory inference, and co-expression network analysis for lncRNAs"
  )
  
  # Print help if requested or no args
  if (length(args) == 0 || "--help" %in% args || "-h" %in% args) {
    cat("\nUsage: scLncR LncExplore [options]\n")
    cat("LncExplore: Multi-modal lncRNA functional exploration pipeline\n\n")
    optparse::print_help(parser)
    q(status = 0)
  }
  
  # Parse arguments
  opt <- optparse::parse_args(parser, args = args, print_help_and_exit = FALSE)
  
  if (is.null(opt$config)) {
    stop("Configuration file must be specified via -c/--config")
  }
  
  config_path <- normalizePath(opt$config, mustWork = TRUE)
  
  # -----------------------------
  # 2. Load config (with dynamic validation)
  # -----------------------------
  config_loader_path <- file.path(script_dir, "R", "utils", "config_loader_function.R")
  legacy_loader_path <- file.path(script_dir, "R", "utils", "config_loader_LncExplore.R")
  if (file.exists(config_loader_path)) {
    source(config_loader_path, local = TRUE)
  } else if (file.exists(legacy_loader_path)) {
    warning("Using legacy function config loader: ", legacy_loader_path)
    source(legacy_loader_path, local = TRUE)
  } else if (!exists("load_LncExplore_config", mode = "function")) {
    stop("Config loader not found. Expected: ", config_loader_path)
  }
  
  cfg <- load_LncExplore_config(config_path)
  
  # -----------------------------
  # 3. Load input Seurat object
  # -----------------------------
  message("\n📦 Loading input Seurat object...")
  if (!file.exists(cfg$input_seurat)) {
    stop("Input Seurat object not found: ", cfg$input_seurat)
  }
  seu_obj <- readRDS(cfg$input_seurat)
  message("✅ Loaded Seurat object with ", ncol(seu_obj), " cells and ", nrow(seu_obj), " features")
  
  # -----------------------------
  # 4. Module execution dispatcher (with error isolation)
  # -----------------------------
  results <- list(
    location = NULL,
    monocle2 = NULL,
    wgcna = NULL,
    errors = character()
  )
  
  # Helper: Safe module executor
  run_module_safely <- function(module_name, func_name, params, cfg_section) {
    message("\n", strrep("=", 60))
    message("🚀 STARTING MODULE: ", toupper(module_name))
    message(strrep("=", 60))
    
    # Check function availability
    if (!exists(func_name, mode = "function")) {
      err_msg <- sprintf(
        "❌ MODULE FAILED: '%s' function not found. Ensure analysis functions are sourced.",
        func_name
      )
      message(err_msg)
      results$errors <<- c(results$errors, paste(module_name, ":", err_msg))
      return(FALSE)
    }
    
    # Execute with error capture
    tryCatch({
      do.call(func_name, params)
      message("\n✅ MODULE COMPLETED: ", toupper(module_name))
      TRUE
    }, error = function(e) {
      err_detail <- sprintf(
        "❌ MODULE FAILED: %s | Error: %s",
        module_name,
        conditionMessage(e)
      )
      message("\n", err_detail)
      results$errors <<- c(results$errors, paste(module_name, ":", conditionMessage(e)))
      FALSE
    })
  }
  
  # -----------------------------
  # 5. Execute selected modules
  # -----------------------------
  modules_to_run <- cfg$.run_modules  # Injected by config loader
  
  if ("location" %in% modules_to_run) {
    params_loc <- list(
      seu_obj = seu_obj,
      LOG2FC_THRESH = cfg$location$LOG2FC_THRESH,
      PADJ_THRESH = cfg$location$PADJ_THRESH,
      output_path = cfg$location$output_path,
      lncRNA_name = cfg$location$lncRNA_name
    )
    results$location <- run_module_safely(
      "snRNA/scRNA Enrichment Analysis",
      "SnScEnrichmentAnalysis",
      params_loc,
      cfg$location
    )
  }
  
  if ("monocle2" %in% modules_to_run) {
    # Parse hub_genes: "all" → "all", else split to vector
    hub_genes_val <- if (tolower(cfg$monocle2$hub_genes) == "all") {
      "all"
    } else {
      trimws(unlist(strsplit(cfg$monocle2$hub_genes, ",")))
    }
    
    params_mono <- list(
      seu_obj = seu_obj,
      qval = cfg$monocle2$qval,
      reduceDimensionMethod = cfg$monocle2$reduceDimensionMethod,
      output_path = cfg$monocle2$output_path,
      hub_genes = hub_genes_val
    )
    results$monocle2 <- run_module_safely(
      "Monocle2 Trajectory", 
      "monocle2_analysis", 
      params_mono,
      cfg$monocle2
    )
  }
  
  if ("wgcna" %in% modules_to_run) {
    # Normalize cell_types: handle empty list/vector
    cell_types_val <- if (is.null(cfg$wgcna$cell_types) || length(cfg$wgcna$cell_types) == 0) {
      character(0)
    } else {
      as.character(unlist(cfg$wgcna$cell_types))
    }
    
    params_wgcna <- list(
      seu_obj = seu_obj,
      output_path = cfg$wgcna$output_path,
      cell_types = cell_types_val,
      pro_name = cfg$wgcna$pro_name,
      lnc_name = cfg$wgcna$lnc_name,
      gene_select_method = cfg$wgcna$gene_select_method
    )
    results$wgcna <- run_module_safely(
      "WGCNA Network", 
      "DO_WGCNA", 
      params_wgcna,
      cfg$wgcna
    )
  }
  
  # -----------------------------
  # 6. Final summary report
  # -----------------------------
  message("\n", strrep("=", 60))
  message("📊 LNC EXPLORE PIPELINE SUMMARY")
  message(strrep("=", 60))
  
  total <- length(modules_to_run)
  succeeded <- sum(unlist(results[modules_to_run]), na.rm = TRUE)
  failed <- total - succeeded
  
  message(sprintf("• Modules scheduled : %d", total))
  message(sprintf("• Modules succeeded : %d ✅", succeeded))
  message(sprintf("• Modules failed    : %d ❌", failed))
  
  if (failed > 0) {
    message("\n⚠️  FAILED MODULES DETAILS:")
    for (err in results$errors) {
      message("  - ", err)
    }
    stop("Pipeline completed with errors. Check logs above for details.")
  } else {
    message("\n🎉 ALL SELECTED MODULES COMPLETED SUCCESSFULLY!")
    message("   Results saved to respective output directories specified in config.")
  }
}
