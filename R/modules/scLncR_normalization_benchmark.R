`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0 || (is.character(x) && nchar(x[1]) == 0)) {
    return(y)
  }
  x
}

safe_get_assay_data <- function(seu_obj, assay = "RNA", layer = "data") {
  out <- tryCatch({
    Seurat::GetAssayData(seu_obj, assay = assay, layer = layer)
  }, error = function(e) NULL)
  if (!is.null(out)) return(out)

  out <- tryCatch({
    Seurat::GetAssayData(seu_obj, assay = assay, slot = layer)
  }, error = function(e) NULL)
  if (!is.null(out)) return(out)

  out <- tryCatch({
    seu_obj[[assay]][layer]
  }, error = function(e) NULL)
  if (!is.null(out)) return(out)

  stop("Cannot read assay data: assay='", assay, "', layer='", layer, "'.")
}

safe_set_assay_data <- function(seu_obj, assay = "RNA", layer = "data", new_data) {
  out <- tryCatch({
    Seurat::SetAssayData(seu_obj, assay = assay, layer = layer, new.data = new_data)
  }, error = function(e) NULL)
  if (!is.null(out)) return(out)

  out <- tryCatch({
    Seurat::SetAssayData(seu_obj, assay = assay, slot = layer, new.data = new_data)
  }, error = function(e) NULL)
  if (!is.null(out)) return(out)

  out <- tryCatch({
    seu_obj[[assay]][layer] <- new_data
    seu_obj
  }, error = function(e) NULL)
  if (!is.null(out)) return(out)

  stop("Cannot set assay data: assay='", assay, "', layer='", layer, "'.")
}

safe_average_expression <- function(seu_obj, assay, group_by) {
  avg <- tryCatch({
    Seurat::AverageExpression(
      seu_obj,
      assays = assay,
      group.by = group_by,
      layer = "data",
      verbose = FALSE
    )
  }, error = function(e) NULL)

  if (is.null(avg)) {
    avg <- tryCatch({
      Seurat::AverageExpression(
        seu_obj,
        assays = assay,
        group.by = group_by,
        slot = "data",
        verbose = FALSE
      )
    }, error = function(e) NULL)
  }

  if (is.null(avg)) return(NULL)
  if (is.list(avg) && assay %in% names(avg)) return(as.matrix(avg[[assay]]))
  as.matrix(avg)
}

configure_benchmark_future <- function(cfg) {
  if (!requireNamespace("future", quietly = TRUE)) return(invisible(FALSE))
  size_gb <- as.numeric(cfg$future_globals_max_size_gb %||% 50)
  if (is.finite(size_gb) && size_gb > 0) {
    options(future.globals.maxSize = size_gb * 1024^3)
  }
  plan_name <- tolower(as.character(cfg$future_plan %||% "sequential"))
  if (plan_name == "sequential") {
    future::plan(future::sequential)
  }
  invisible(TRUE)
}

summarize_values_sparse_aware <- function(mat) {
  total <- length(mat)
  if (total == 0) {
    return(list(mean = NA_real_, median = NA_real_, sd = NA_real_, zero_fraction = NA_real_))
  }
  if (inherits(mat, "sparseMatrix")) {
    x <- mat@x
    n_nonzero <- length(x)
    n_zero <- total - n_nonzero
    mean_val <- if (n_nonzero == 0) 0 else sum(x, na.rm = TRUE) / total
    zero_fraction <- n_zero / total
    median_val <- if (zero_fraction >= 0.5) 0 else stats::median(as.numeric(mat), na.rm = TRUE)
    sum_sq <- if (n_nonzero == 0) 0 else sum(x^2, na.rm = TRUE)
    sd_val <- if (total > 1) sqrt(max((sum_sq - total * mean_val^2) / (total - 1), 0)) else NA_real_
    return(list(mean = mean_val, median = median_val, sd = sd_val, zero_fraction = zero_fraction))
  }
  vals <- as.numeric(mat)
  list(
    mean = mean(vals, na.rm = TRUE),
    median = stats::median(vals, na.rm = TRUE),
    sd = stats::sd(vals, na.rm = TRUE),
    zero_fraction = mean(vals == 0, na.rm = TRUE)
  )
}

read_normalization_benchmark_config <- function(config_path) {
  if (!requireNamespace("yaml", quietly = TRUE)) {
    stop("Package 'yaml' is required. Install with: install.packages('yaml')")
  }
  if (!file.exists(config_path)) stop("Config file not found: ", config_path)

  cfg <- yaml::read_yaml(config_path)
  required_fields <- c("input_seurat", "output_dir", "lnc_name")
  missing_fields <- setdiff(required_fields, names(cfg))
  if (length(missing_fields) > 0) {
    stop("Missing required configuration fields: ", paste(missing_fields, collapse = ", "))
  }

  if (!file.exists(cfg$input_seurat)) {
    stop("Input Seurat RDS not found: ", cfg$input_seurat)
  }
  if (!is.character(cfg$lnc_name) || nchar(cfg$lnc_name) == 0) {
    stop("'lnc_name' must be a non-empty character string.")
  }

  cfg$group_by <- cfg$group_by %||% "cell_type"
  cfg$condition_col <- cfg$condition_col %||% "Group"
  cfg$cluster_col <- cfg$cluster_col %||% "seurat_clusters"
  cfg$strategies <- cfg$strategies %||% c("joint_lognormalize", "separate_lognormalize")
  cfg$marker_test <- cfg$marker_test %||% "wilcox"
  cfg$logfc_threshold <- as.numeric(cfg$logfc_threshold %||% 0.25)
  cfg$min_pct <- as.numeric(cfg$min_pct %||% 0.1)
  cfg$padj_threshold <- as.numeric(cfg$padj_threshold %||% 0.05)
  cfg$top_n_markers <- as.integer(cfg$top_n_markers %||% 50)
  cfg$random_seed <- as.integer(cfg$random_seed %||% 1234)
  cfg$future_plan <- cfg$future_plan %||% "sequential"
  cfg$future_globals_max_size_gb <- as.numeric(cfg$future_globals_max_size_gb %||% 50)
  cfg$max_cells_per_ident <- as.integer(cfg$max_cells_per_ident %||% 0)
  cfg$save_normalized_objects <- isTRUE(cfg$save_normalized_objects %||% TRUE)

  cfg$strategies <- unique(as.character(unlist(cfg$strategies)))
  cfg
}

normalize_with_strategy <- function(seu_obj, strategy, lnc_name, random_seed = 1234) {
  lnc_pattern <- paste0("^", lnc_name)

  if (strategy == "joint_lognormalize") {
    obj <- seu_obj
    Seurat::DefaultAssay(obj) <- "RNA"
    obj <- Seurat::NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
    return(list(object = obj, assay = "RNA", message = "OK"))
  }

  if (strategy == "separate_lognormalize") {
    obj <- seu_obj
    counts_mat <- safe_get_assay_data(obj, assay = "RNA", layer = "counts")
    lnc_idx <- grepl(lnc_pattern, rownames(counts_mat))
    if (sum(lnc_idx) == 0) {
      stop("No lncRNA features found for pattern: ", lnc_pattern)
    }
    if (sum(!lnc_idx) == 0) {
      stop("No non-lncRNA features left after lncRNA split.")
    }

    seu_lnc <- Seurat::CreateSeuratObject(counts = counts_mat[lnc_idx, , drop = FALSE], meta.data = obj@meta.data)
    seu_gene <- Seurat::CreateSeuratObject(counts = counts_mat[!lnc_idx, , drop = FALSE], meta.data = obj@meta.data)
    seu_lnc <- Seurat::NormalizeData(seu_lnc, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
    seu_gene <- Seurat::NormalizeData(seu_gene, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)

    data_bind <- rbind(
      safe_get_assay_data(seu_gene, assay = "RNA", layer = "data"),
      safe_get_assay_data(seu_lnc, assay = "RNA", layer = "data")
    )
    data_bind <- data_bind[rownames(counts_mat), colnames(counts_mat), drop = FALSE]
    obj <- safe_set_assay_data(obj, assay = "RNA", layer = "data", new_data = data_bind)
    return(list(object = obj, assay = "RNA", message = "OK"))
  }

  if (strategy == "sctransform") {
    if (!"SCTransform" %in% getNamespaceExports("Seurat")) {
      stop("Seurat::SCTransform is not available in current Seurat version.")
    }
    obj <- seu_obj
    set.seed(random_seed)
    obj <- Seurat::SCTransform(
      obj,
      assay = "RNA",
      new.assay.name = "SCT",
      return.only.var.genes = FALSE,
      verbose = FALSE
    )
    Seurat::DefaultAssay(obj) <- "SCT"
    return(list(object = obj, assay = "SCT", message = "OK"))
  }

  if (strategy == "scran") {
    if (!requireNamespace("SingleCellExperiment", quietly = TRUE) ||
        !requireNamespace("scran", quietly = TRUE) ||
        !requireNamespace("scuttle", quietly = TRUE) ||
        !requireNamespace("SummarizedExperiment", quietly = TRUE)) {
      stop("scran strategy requires SingleCellExperiment, scran, scuttle, and SummarizedExperiment.")
    }

    obj <- seu_obj
    counts_mat <- safe_get_assay_data(obj, assay = "RNA", layer = "counts")
    sce <- SingleCellExperiment::SingleCellExperiment(list(counts = counts_mat))
    set.seed(random_seed)
    quick_clusters <- scran::quickCluster(sce)
    sce <- scran::computeSumFactors(sce, clusters = quick_clusters)
    sce <- scuttle::logNormCounts(sce)
    logcounts_mat <- SummarizedExperiment::assay(sce, "logcounts")
    logcounts_mat <- logcounts_mat[rownames(counts_mat), colnames(counts_mat), drop = FALSE]
    obj <- safe_set_assay_data(obj, assay = "RNA", layer = "data", new_data = logcounts_mat)
    return(list(object = obj, assay = "RNA", message = "OK"))
  }

  stop("Unknown strategy: ", strategy)
}

safe_find_all_markers <- function(
  seu_obj,
  assay_used,
  group_by,
  marker_test = "wilcox",
  logfc_threshold = 0.25,
  min_pct = 0.1,
  random_seed = 1234,
  max_cells_per_ident = 0
) {
  if (!(group_by %in% colnames(seu_obj@meta.data))) {
    return(data.frame())
  }

  idents <- as.character(seu_obj@meta.data[[group_by]])
  idents[is.na(idents)] <- "NA_group"
  if (length(unique(idents)) < 2) {
    return(data.frame())
  }
  if (is.finite(max_cells_per_ident) && max_cells_per_ident > 0) {
    set.seed(random_seed)
    keep_idx <- unlist(
      lapply(split(seq_along(idents), idents), function(idx) {
        sample(idx, size = min(length(idx), max_cells_per_ident))
      }),
      use.names = FALSE
    )
    keep_idx <- sort(unique(keep_idx))
    seu_obj <- subset(seu_obj, cells = colnames(seu_obj)[keep_idx])
    idents <- idents[keep_idx]
  }
  Seurat::Idents(seu_obj) <- idents

  markers <- tryCatch({
    Seurat::FindAllMarkers(
      object = seu_obj,
      assay = assay_used,
      slot = "data",
      only.pos = FALSE,
      test.use = marker_test,
      logfc.threshold = logfc_threshold,
      min.pct = min_pct,
      return.thresh = 1,
      random.seed = random_seed,
      verbose = FALSE
    )
  }, error = function(e) {
    data.frame()
  })

  if (nrow(markers) == 0) return(markers)
  markers <- as.data.frame(markers)
  if (!("gene" %in% colnames(markers))) {
    markers$gene <- rownames(markers)
  }
  markers
}

safe_find_lnc_hvg <- function(seu_obj, assay_used, lnc_genes) {
  obj <- seu_obj
  out <- tryCatch({
    obj <- Seurat::FindVariableFeatures(
      object = obj,
      assay = assay_used,
      selection.method = "vst",
      nfeatures = min(2000, nrow(obj)),
      verbose = FALSE
    )
    intersect(Seurat::VariableFeatures(obj, assay = assay_used), lnc_genes)
  }, error = function(e) character(0))
  unique(out)
}

compute_jaccard <- function(set_a, set_b) {
  ua <- unique(set_a)
  ub <- unique(set_b)
  uni <- union(ua, ub)
  if (length(uni) == 0) return(NA_real_)
  length(intersect(ua, ub)) / length(uni)
}

pairwise_jaccard_table <- function(sets_named) {
  nms <- names(sets_named)
  out <- data.frame(
    strategy_1 = character(),
    strategy_2 = character(),
    size_1 = integer(),
    size_2 = integer(),
    intersection = integer(),
    union = integer(),
    jaccard = numeric(),
    stringsAsFactors = FALSE
  )

  for (i in seq_along(nms)) {
    for (j in seq_along(nms)) {
      s1 <- unique(sets_named[[nms[i]]])
      s2 <- unique(sets_named[[nms[j]]])
      out <- rbind(
        out,
        data.frame(
          strategy_1 = nms[i],
          strategy_2 = nms[j],
          size_1 = length(s1),
          size_2 = length(s2),
          intersection = length(intersect(s1, s2)),
          union = length(union(s1, s2)),
          jaccard = compute_jaccard(s1, s2),
          stringsAsFactors = FALSE
        )
      )
    }
  }
  out
}

build_top_lnc_rank <- function(markers_df, lnc_pattern, padj_threshold, top_n_markers) {
  if (nrow(markers_df) == 0) return(data.frame(gene = character(), rank = integer(), score = numeric()))
  fc_col <- if ("avg_log2FC" %in% colnames(markers_df)) {
    "avg_log2FC"
  } else if ("avg_logFC" %in% colnames(markers_df)) {
    "avg_logFC"
  } else {
    NULL
  }
  if (is.null(fc_col) || !("p_val_adj" %in% colnames(markers_df))) {
    return(data.frame(gene = character(), rank = integer(), score = numeric()))
  }

  d <- markers_df[grepl(lnc_pattern, markers_df$gene), , drop = FALSE]
  d <- d[is.finite(d$p_val_adj) & d$p_val_adj <= padj_threshold, , drop = FALSE]
  if (nrow(d) == 0) return(data.frame(gene = character(), rank = integer(), score = numeric()))

  d$abs_fc <- abs(as.numeric(d[[fc_col]]))
  d <- d[order(d$p_val_adj, -d$abs_fc), , drop = FALSE]
  d <- d[!duplicated(d$gene), , drop = FALSE]
  d <- head(d, top_n_markers)
  d$rank <- seq_len(nrow(d))
  d$score <- d$abs_fc
  d[, c("gene", "rank", "score"), drop = FALSE]
}

pairwise_rank_correlation <- function(rank_lists_named) {
  nms <- names(rank_lists_named)
  out <- data.frame(
    strategy_1 = character(),
    strategy_2 = character(),
    n_shared_top_lnc = integer(),
    spearman_rank_cor = numeric(),
    stringsAsFactors = FALSE
  )

  for (i in seq_along(nms)) {
    for (j in seq_along(nms)) {
      d1 <- rank_lists_named[[nms[i]]]
      d2 <- rank_lists_named[[nms[j]]]
      common <- intersect(d1$gene, d2$gene)
      rho <- NA_real_
      if (length(common) >= 3) {
        r1 <- d1$rank[match(common, d1$gene)]
        r2 <- d2$rank[match(common, d2$gene)]
        rho <- suppressWarnings(stats::cor(r1, r2, method = "spearman"))
      }
      out <- rbind(
        out,
        data.frame(
          strategy_1 = nms[i],
          strategy_2 = nms[j],
          n_shared_top_lnc = length(common),
          spearman_rank_cor = rho,
          stringsAsFactors = FALSE
        )
      )
    }
  }
  out
}

matrix_pairwise_correlation <- function(mat_a, mat_b) {
  if (is.null(mat_a) || is.null(mat_b)) return(list(correlation = NA_real_, n_genes = 0L, n_groups = 0L))
  common_genes <- intersect(rownames(mat_a), rownames(mat_b))
  common_groups <- intersect(colnames(mat_a), colnames(mat_b))
  if (length(common_genes) < 10 || length(common_groups) < 2) {
    return(list(correlation = NA_real_, n_genes = length(common_genes), n_groups = length(common_groups)))
  }
  v1 <- as.vector(mat_a[common_genes, common_groups, drop = FALSE])
  v2 <- as.vector(mat_b[common_genes, common_groups, drop = FALSE])
  list(
    correlation = suppressWarnings(stats::cor(v1, v2, method = "spearman", use = "pairwise.complete.obs")),
    n_genes = length(common_genes),
    n_groups = length(common_groups)
  )
}

write_normalization_benchmark_report <- function(
  out_path,
  cfg,
  input_summary,
  qc_summary,
  marker_summary,
  skipped_table
) {
  lines <- c(
    "# Normalization Benchmark Report",
    "",
    "## Input",
    paste0("- Input Seurat object: `", cfg$input_seurat, "`"),
    paste0("- Cells: ", input_summary$n_cells),
    paste0("- Features: ", input_summary$n_features),
    paste0("- lncRNA prefix: `", cfg$lnc_name, "`"),
    paste0("- group_by: `", cfg$group_by, "`"),
    paste0("- condition_col: `", cfg$condition_col, "`"),
    paste0("- cluster_col: `", cfg$cluster_col, "`"),
    paste0("- random_seed: ", cfg$random_seed),
    paste0("- future_plan: `", cfg$future_plan, "`"),
    paste0("- future_globals_max_size_gb: ", cfg$future_globals_max_size_gb),
    paste0("- max_cells_per_ident: ", cfg$max_cells_per_ident),
    paste0("- save_normalized_objects: ", cfg$save_normalized_objects),
    "",
    "## Strategy Status",
    "- See `stability/strategy_qc_summary.csv` for per-strategy status and error messages.",
    "",
    "## Marker Summary",
    "- See `markers/lncRNA_marker_summary.csv` for total markers, lncRNA markers, and significant lncRNA markers.",
    "- Use these values to evaluate whether separate normalization changes lncRNA marker detection.",
    "",
    "## Stability Summary",
    "- `stability/marker_overlap_jaccard.csv`: overall marker set overlap across strategies.",
    "- `stability/lnc_marker_overlap_jaccard.csv`: lncRNA marker overlap across strategies.",
    "- `stability/top_lnc_marker_rank_correlation.csv`: top lncRNA marker rank consistency.",
    "- `stability/strategy_pairwise_correlation_matrix.csv`: pseudo-bulk/cluster-level expression correlation proxies.",
    "- `stability/lnc_hvg_overlap_jaccard.csv`: lncRNA HVG overlap across strategies.",
    "",
    "## Interpretation Guidance",
    "- This benchmark compares normalization strategies and quantifies stability of lncRNA-associated signals.",
    "- Results should be interpreted as sensitivity/stability evidence, not proof that one strategy is universally superior.",
    "- snRNA/scRNA enrichment, lncRNA markers, pseudotime and WGCNA interpretations should remain conservative and data-dependent.",
    "",
    "## Skipped Strategies",
    if (nrow(skipped_table) == 0) "- None" else paste0("- ", skipped_table$strategy, ": ", skipped_table$error_message),
    "",
    "## Reproducibility",
    "- Configuration snapshot: `config_snapshot.yaml`",
    "- Session info: `sessionInfo.txt`"
  )
  writeLines(lines, con = out_path)
}

run_normalization_benchmark <- function(cfg) {
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Package 'Seurat' is required.")
  }
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required.")
  }
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Package 'dplyr' is required.")
  }
  if (!requireNamespace("reshape2", quietly = TRUE)) {
    stop("Package 'reshape2' is required.")
  }
  if (!requireNamespace("yaml", quietly = TRUE)) {
    stop("Package 'yaml' is required.")
  }

  configure_benchmark_future(cfg)
  set.seed(cfg$random_seed)
  output_dir <- cfg$output_dir
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(output_dir, "normalized_objects"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(output_dir, "markers"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(output_dir, "stability"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(output_dir, "figures"), recursive = TRUE, showWarnings = FALSE)

  yaml::write_yaml(cfg, file.path(output_dir, "config_snapshot.yaml"))
  writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo.txt"))

  seu_input <- readRDS(cfg$input_seurat)
  if (!inherits(seu_input, "Seurat")) {
    stop("Input RDS is not a Seurat object.")
  }
  if (!("RNA" %in% names(seu_input@assays))) {
    stop("Input Seurat object does not contain RNA assay.")
  }

  counts_input <- safe_get_assay_data(seu_input, assay = "RNA", layer = "counts")
  lnc_pattern <- paste0("^", cfg$lnc_name)
  lnc_idx_global <- grepl(lnc_pattern, rownames(counts_input))
  if (sum(lnc_idx_global) == 0) {
    stop("No lncRNA features matched prefix pattern: ", lnc_pattern)
  }
  lnc_genes_global <- rownames(counts_input)[lnc_idx_global]

  group_by <- cfg$group_by
  if (!(group_by %in% colnames(seu_input@meta.data))) {
    if (cfg$cluster_col %in% colnames(seu_input@meta.data)) {
      warning("group_by='", group_by, "' not found; fallback to cluster_col='", cfg$cluster_col, "'.")
      group_by <- cfg$cluster_col
    } else {
      stop("Neither group_by nor cluster_col exists in metadata.")
    }
  }

  strategy_qc <- data.frame(
    strategy = character(),
    status = character(),
    assay_used = character(),
    n_cells = integer(),
    n_features = integer(),
    n_lnc_features = integer(),
    error_message = character(),
    stringsAsFactors = FALSE
  )

  marker_summary <- data.frame(
    strategy = character(),
    total_markers = integer(),
    lnc_markers = integer(),
    lnc_marker_fraction = numeric(),
    significant_lnc_markers = integer(),
    stringsAsFactors = FALSE
  )

  lnc_signal_summary <- data.frame(
    strategy = character(),
    feature_type = character(),
    mean = numeric(),
    median = numeric(),
    sd = numeric(),
    zero_fraction = numeric(),
    stringsAsFactors = FALSE
  )

  strategy_markers <- list()
  strategy_sig_marker_sets <- list()
  strategy_sig_lnc_marker_sets <- list()
  strategy_top_lnc_ranks <- list()
  strategy_lnc_hvg_sets <- list()
  strategy_group_avg <- list()
  strategy_cluster_avg <- list()

  for (strategy in cfg$strategies) {
    res <- tryCatch({
      norm_result <- normalize_with_strategy(
        seu_obj = seu_input,
        strategy = strategy,
        lnc_name = cfg$lnc_name,
        random_seed = cfg$random_seed
      )

      obj_norm <- norm_result$object
      assay_used <- norm_result$assay
      Seurat::DefaultAssay(obj_norm) <- assay_used
      if (isTRUE(cfg$save_normalized_objects)) {
        saveRDS(obj_norm, file.path(output_dir, "normalized_objects", paste0(strategy, ".rds")))
      }

      data_mat <- safe_get_assay_data(obj_norm, assay = assay_used, layer = "data")
      rownames_data <- rownames(data_mat)
      lnc_idx <- grepl(lnc_pattern, rownames_data)
      lnc_stats <- summarize_values_sparse_aware(data_mat[lnc_idx, , drop = FALSE])
      mrna_stats <- summarize_values_sparse_aware(data_mat[!lnc_idx, , drop = FALSE])

      lnc_signal_summary <- rbind(
        lnc_signal_summary,
        data.frame(
          strategy = strategy,
          feature_type = "lncRNA",
          mean = lnc_stats$mean,
          median = lnc_stats$median,
          sd = lnc_stats$sd,
          zero_fraction = lnc_stats$zero_fraction,
          stringsAsFactors = FALSE
        ),
        data.frame(
          strategy = strategy,
          feature_type = "mRNA",
          mean = mrna_stats$mean,
          median = mrna_stats$median,
          sd = mrna_stats$sd,
          zero_fraction = mrna_stats$zero_fraction,
          stringsAsFactors = FALSE
        )
      )

      markers_df <- safe_find_all_markers(
        seu_obj = obj_norm,
        assay_used = assay_used,
        group_by = group_by,
        marker_test = cfg$marker_test,
        logfc_threshold = cfg$logfc_threshold,
        min_pct = cfg$min_pct,
        random_seed = cfg$random_seed,
        max_cells_per_ident = cfg$max_cells_per_ident
      )

      write.csv(
        markers_df,
        file.path(output_dir, "markers", paste0(strategy, "_all_markers.csv")),
        row.names = FALSE
      )

      if (nrow(markers_df) > 0) {
        sig_markers <- unique(markers_df$gene[is.finite(markers_df$p_val_adj) & markers_df$p_val_adj <= cfg$padj_threshold])
        lnc_markers <- unique(markers_df$gene[grepl(lnc_pattern, markers_df$gene)])
        sig_lnc_markers <- unique(markers_df$gene[
          is.finite(markers_df$p_val_adj) &
            markers_df$p_val_adj <= cfg$padj_threshold &
            grepl(lnc_pattern, markers_df$gene)
        ])
      } else {
        sig_markers <- character(0)
        lnc_markers <- character(0)
        sig_lnc_markers <- character(0)
      }

      marker_summary <- rbind(
        marker_summary,
        data.frame(
          strategy = strategy,
          total_markers = nrow(markers_df),
          lnc_markers = length(lnc_markers),
          lnc_marker_fraction = ifelse(nrow(markers_df) == 0, NA_real_, length(lnc_markers) / nrow(markers_df)),
          significant_lnc_markers = length(sig_lnc_markers),
          stringsAsFactors = FALSE
        )
      )

      strategy_markers[[strategy]] <- markers_df
      strategy_sig_marker_sets[[strategy]] <- sig_markers
      strategy_sig_lnc_marker_sets[[strategy]] <- sig_lnc_markers
      strategy_top_lnc_ranks[[strategy]] <- build_top_lnc_rank(
        markers_df = markers_df,
        lnc_pattern = lnc_pattern,
        padj_threshold = cfg$padj_threshold,
        top_n_markers = cfg$top_n_markers
      )
      strategy_lnc_hvg_sets[[strategy]] <- safe_find_lnc_hvg(obj_norm, assay_used, lnc_genes_global)
      strategy_group_avg[[strategy]] <- safe_average_expression(obj_norm, assay_used, group_by)
      if (cfg$cluster_col %in% colnames(obj_norm@meta.data)) {
        strategy_cluster_avg[[strategy]] <- safe_average_expression(obj_norm, assay_used, cfg$cluster_col)
      } else {
        strategy_cluster_avg[[strategy]] <- NULL
      }

      strategy_qc <- rbind(
        strategy_qc,
        data.frame(
          strategy = strategy,
          status = "success",
          assay_used = assay_used,
          n_cells = ncol(obj_norm),
          n_features = nrow(obj_norm),
          n_lnc_features = sum(lnc_idx),
          error_message = "",
          stringsAsFactors = FALSE
        )
      )

      TRUE
    }, error = function(e) {
      strategy_qc <<- rbind(
        strategy_qc,
        data.frame(
          strategy = strategy,
          status = "skipped",
          assay_used = NA_character_,
          n_cells = ncol(seu_input),
          n_features = nrow(seu_input),
          n_lnc_features = sum(lnc_idx_global),
          error_message = conditionMessage(e),
          stringsAsFactors = FALSE
        )
      )
      FALSE
    })
  }

  write.csv(marker_summary, file.path(output_dir, "markers", "lncRNA_marker_summary.csv"), row.names = FALSE)
  write.csv(lnc_signal_summary, file.path(output_dir, "stability", "per_strategy_lnc_signal_summary.csv"), row.names = FALSE)
  write.csv(strategy_qc, file.path(output_dir, "stability", "strategy_qc_summary.csv"), row.names = FALSE)

  successful <- strategy_qc$strategy[strategy_qc$status == "success"]
  if (length(successful) == 0) {
    write_normalization_benchmark_report(
      out_path = file.path(output_dir, "normalization_benchmark_report.md"),
      cfg = cfg,
      input_summary = list(n_cells = ncol(seu_input), n_features = nrow(seu_input)),
      qc_summary = strategy_qc,
      marker_summary = marker_summary,
      skipped_table = strategy_qc[strategy_qc$status == "skipped", , drop = FALSE]
    )
    stop("No normalization strategy completed successfully.")
  }

  sig_marker_sets <- strategy_sig_marker_sets[successful]
  sig_lnc_sets <- strategy_sig_lnc_marker_sets[successful]
  lnc_rank_lists <- strategy_top_lnc_ranks[successful]
  lnc_hvg_sets <- strategy_lnc_hvg_sets[successful]

  marker_jaccard <- pairwise_jaccard_table(sig_marker_sets)
  lnc_marker_jaccard <- pairwise_jaccard_table(sig_lnc_sets)
  lnc_hvg_jaccard <- pairwise_jaccard_table(lnc_hvg_sets)
  rank_cor <- pairwise_rank_correlation(lnc_rank_lists)

  write.csv(marker_jaccard, file.path(output_dir, "stability", "marker_overlap_jaccard.csv"), row.names = FALSE)
  write.csv(lnc_marker_jaccard, file.path(output_dir, "stability", "lnc_marker_overlap_jaccard.csv"), row.names = FALSE)
  write.csv(rank_cor, file.path(output_dir, "stability", "top_lnc_marker_rank_correlation.csv"), row.names = FALSE)
  write.csv(lnc_hvg_jaccard, file.path(output_dir, "stability", "lnc_hvg_overlap_jaccard.csv"), row.names = FALSE)

  pair_corr <- data.frame(
    strategy_1 = character(),
    strategy_2 = character(),
    metric = character(),
    correlation = numeric(),
    n_genes = integer(),
    n_groups = integer(),
    stringsAsFactors = FALSE
  )

  for (s1 in successful) {
    for (s2 in successful) {
      cg <- matrix_pairwise_correlation(strategy_group_avg[[s1]], strategy_group_avg[[s2]])
      pair_corr <- rbind(
        pair_corr,
        data.frame(
          strategy_1 = s1,
          strategy_2 = s2,
          metric = "group_by_average_expression",
          correlation = cg$correlation,
          n_genes = as.integer(cg$n_genes),
          n_groups = as.integer(cg$n_groups),
          stringsAsFactors = FALSE
        )
      )

      cc <- matrix_pairwise_correlation(strategy_cluster_avg[[s1]], strategy_cluster_avg[[s2]])
      pair_corr <- rbind(
        pair_corr,
        data.frame(
          strategy_1 = s1,
          strategy_2 = s2,
          metric = "cluster_average_expression",
          correlation = cc$correlation,
          n_genes = as.integer(cc$n_genes),
          n_groups = as.integer(cc$n_groups),
          stringsAsFactors = FALSE
        )
      )
    }
  }
  write.csv(pair_corr, file.path(output_dir, "stability", "strategy_pairwise_correlation_matrix.csv"), row.names = FALSE)

  # figure 1: lnc marker counts by strategy
  marker_plot_df <- marker_summary
  if (nrow(marker_plot_df) > 0) {
    marker_long <- rbind(
      data.frame(strategy = marker_plot_df$strategy, metric = "lnc_markers", value = marker_plot_df$lnc_markers),
      data.frame(strategy = marker_plot_df$strategy, metric = "significant_lnc_markers", value = marker_plot_df$significant_lnc_markers)
    )

    p_counts <- ggplot2::ggplot(marker_long, ggplot2::aes(x = strategy, y = value, fill = metric)) +
      ggplot2::geom_col(position = "dodge") +
      ggplot2::labs(
        title = "lncRNA Marker Counts by Normalization Strategy",
        x = "Strategy",
        y = "Count",
        fill = "Metric"
      ) +
      ggplot2::theme_minimal(base_size = 11) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30, hjust = 1))

    ggplot2::ggsave(file.path(output_dir, "figures", "lnc_marker_counts_by_strategy.pdf"), p_counts, width = 8, height = 5, dpi = 300)
    ggplot2::ggsave(file.path(output_dir, "figures", "lnc_marker_counts_by_strategy.png"), p_counts, width = 8, height = 5, dpi = 300)
  }

  # figure 2: lnc marker overlap heatmap
  if (nrow(lnc_marker_jaccard) > 0) {
    p_overlap <- ggplot2::ggplot(
      lnc_marker_jaccard,
      ggplot2::aes(x = strategy_1, y = strategy_2, fill = jaccard)
    ) +
      ggplot2::geom_tile() +
      ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f", jaccard)), size = 3) +
      ggplot2::scale_fill_gradient(low = "#F7FBFF", high = "#08306B", na.value = "grey90") +
      ggplot2::labs(
        title = "lncRNA Marker Overlap (Jaccard)",
        x = "Strategy",
        y = "Strategy",
        fill = "Jaccard"
      ) +
      ggplot2::theme_minimal(base_size = 11) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30, hjust = 1))

    ggplot2::ggsave(file.path(output_dir, "figures", "lnc_marker_overlap_heatmap.pdf"), p_overlap, width = 6.5, height = 5.5, dpi = 300)
    ggplot2::ggsave(file.path(output_dir, "figures", "lnc_marker_overlap_heatmap.png"), p_overlap, width = 6.5, height = 5.5, dpi = 300)
  }

  # figure 3: lnc signal distribution by strategy
  if (nrow(lnc_signal_summary) > 0) {
    p_signal <- ggplot2::ggplot(
      lnc_signal_summary,
      ggplot2::aes(x = strategy, y = mean, fill = feature_type)
    ) +
      ggplot2::geom_col(position = "dodge") +
      ggplot2::labs(
        title = "Mean Normalized Expression by Strategy",
        x = "Strategy",
        y = "Mean normalized expression",
        fill = "Feature type"
      ) +
      ggplot2::theme_minimal(base_size = 11) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30, hjust = 1))

    ggplot2::ggsave(file.path(output_dir, "figures", "lnc_signal_distribution_by_strategy.pdf"), p_signal, width = 8, height = 5, dpi = 300)
    ggplot2::ggsave(file.path(output_dir, "figures", "lnc_signal_distribution_by_strategy.png"), p_signal, width = 8, height = 5, dpi = 300)
  }

  # figure 4: normalization qc summary
  qc_plot_df <- strategy_qc
  qc_plot_df$success_flag <- ifelse(qc_plot_df$status == "success", 1, 0)
  p_qc <- ggplot2::ggplot(qc_plot_df, ggplot2::aes(x = strategy, y = n_lnc_features, fill = as.factor(success_flag))) +
    ggplot2::geom_col() +
    ggplot2::scale_fill_manual(values = c("0" = "#D7301F", "1" = "#1A9850"), labels = c("0" = "skipped", "1" = "success")) +
    ggplot2::labs(
      title = "Normalization Strategy QC Summary",
      x = "Strategy",
      y = "Detected lncRNA features",
      fill = "Status"
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30, hjust = 1))

  ggplot2::ggsave(file.path(output_dir, "figures", "normalization_qc_summary.pdf"), p_qc, width = 8, height = 5, dpi = 300)
  ggplot2::ggsave(file.path(output_dir, "figures", "normalization_qc_summary.png"), p_qc, width = 8, height = 5, dpi = 300)

  write_normalization_benchmark_report(
    out_path = file.path(output_dir, "normalization_benchmark_report.md"),
    cfg = cfg,
    input_summary = list(n_cells = ncol(seu_input), n_features = nrow(seu_input)),
    qc_summary = strategy_qc,
    marker_summary = marker_summary,
    skipped_table = strategy_qc[strategy_qc$status == "skipped", , drop = FALSE]
  )

  invisible(
    list(
      strategy_qc_summary = strategy_qc,
      marker_summary = marker_summary,
      marker_overlap_jaccard = marker_jaccard,
      lnc_marker_overlap_jaccard = lnc_marker_jaccard,
      rank_correlation = rank_cor
    )
  )
}
