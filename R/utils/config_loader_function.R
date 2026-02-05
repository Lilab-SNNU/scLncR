# =============================================================================
# scLncR - Configuration Loader & Validator for LncExplore Module
# Author: [Yin SW]
# Description: Safely loads and validates nested YAML config for Location/Monocle2/WGCNA analyses.
# Dependencies: yaml
# =============================================================================


load_LncExplore_config <- function(config_path) {
  # Load required packages
  if (!requireNamespace("yaml", quietly = TRUE)) {
    stop("Package 'yaml' is required but not installed. Please run: install.packages('yaml')")
  }

  # Check if config file exists
  if (!file.exists(config_path)) {
    stop("Configuration file not found: ", config_path)
  }

  # Read YAML
  cfg <- yaml::read_yaml(config_path)

  # -----------------------------
  # Required top-level fields
  # -----------------------------
  required_top <- c("input_seurat", "location", "monocle2", "wgcna")
  missing_top <- setdiff(required_top, names(cfg))
  if (length(missing_top) > 0) {
    stop("Missing required top-level configuration sections: ", paste(missing_top, collapse = ", "))
  }

  # -----------------------------
  # Global input validation
  # -----------------------------
  if (!file.exists(cfg$input_seurat)) {
    stop("Input Seurat object not found: ", cfg$input_seurat)
  }
  if (!grepl("\\.rds$", cfg$input_seurat, ignore.case = TRUE)) {
    warning("Input file does not have .rds extension. Ensure it is a valid Seurat object.")
  }

  # -----------------------------
  # Helper: Validate numeric in range
  # -----------------------------
  validate_numeric_range <- function(x, name, min_val = -Inf, max_val = Inf, allow_zero = FALSE) {
    if (!is.numeric(x) || length(x) != 1 || is.na(x)) {
      stop("'", name, "' must be a single numeric value.")
    }
    lower_bound <- if (allow_zero) min_val else max(min_val, .Machine$double.eps)
    if (x < lower_bound || x > max_val) {
      stop("'", name, "' must be in range (", 
           ifelse(allow_zero, "[", "("), lower_bound, ", ", max_val, 
           ifelse(max_val == Inf, ")", "]"), ". Got: ", x)
    }
  }

  # -----------------------------
  # Location Analysis Validation
  # -----------------------------
  loc <- cfg$location
  required_loc <- c("LOG2FC_THRESH", "PADJ_THRESH", "output_path", "lncRNA_name")
  missing_loc <- setdiff(required_loc, names(loc))
  if (length(missing_loc) > 0) stop("Missing 'location' parameters: ", paste(missing_loc, collapse = ", "))
  
  validate_numeric_range(loc$LOG2FC_THRESH, "location$LOG2FC_THRESH", min_val = 0)
  validate_numeric_range(loc$PADJ_THRESH, "location$PADJ_THRESH", min_val = 0, max_val = 1)
  
  if (!is.character(loc$output_path) || nchar(loc$output_path) == 0) {
    stop("'location$output_path' must be a non-empty string.")
  }
  if (!dir.exists(dirname(loc$output_path))) {
    message("Creating location output directory parent: ", dirname(loc$output_path))
    dir.create(dirname(loc$output_path), recursive = TRUE)
  }
  
  if (!is.character(loc$lncRNA_name) || nchar(loc$lncRNA_name) == 0) {
    stop("'location$lncRNA_name' must be a non-empty string.")
  }
  if (grepl("[^a-zA-Z0-9-]", loc$lncRNA_name)) {
    stop(
      "Invalid characters in 'location$lncRNA_name': '", loc$lncRNA_name, "'. ",
      "Only letters, digits, and hyphens ('-') are allowed."
    )
  }

  # -----------------------------
  # Monocle2 Analysis Validation
  # -----------------------------
  mono <- cfg$monocle2
  required_mono <- c("qval", "reduceDimensionMethod", "output_path", "hub_genes")
  missing_mono <- setdiff(required_mono, names(mono))
  if (length(missing_mono) > 0) stop("Missing 'monocle2' parameters: ", paste(missing_mono, collapse = ", "))
  
  validate_numeric_range(mono$qval, "monocle2$qval", min_val = 0, max_val = 1)
  
  allowed_methods <- c("DDRTree", "ICA", "tSNE")
  if (!mono$reduceDimensionMethod %in% allowed_methods) {
    stop(
      "Invalid 'monocle2$reduceDimensionMethod': '", mono$reduceDimensionMethod, "'. ",
      "Must be one of: ", paste(allowed_methods, collapse = ", ")
    )
  }
  
  if (!is.character(mono$output_path) || nchar(mono$output_path) == 0) {
    stop("'monocle2$output_path' must be a non-empty string.")
  }
  if (!dir.exists(dirname(mono$output_path))) {
    message("Creating Monocle2 output directory parent: ", dirname(mono$output_path))
    dir.create(dirname(mono$output_path), recursive = TRUE)
  }
  
  if (!is.character(mono$hub_genes) || nchar(mono$hub_genes) == 0) {
    stop("'monocle2$hub_genes' must be a non-empty string ('all' or comma-separated genes).")
  }

  # -----------------------------
  # WGCNA Analysis Validation
  # -----------------------------
  wgcna <- cfg$wgcna
  required_wgcna <- c("output_path", "cell_types", "pro_name", "lnc_name", "gene_select_method")
  missing_wgcna <- setdiff(required_wgcna, names(wgcna))
  if (length(missing_wgcna) > 0) stop("Missing 'wgcna' parameters: ", paste(missing_wgcna, collapse = ", "))
  
  if (!is.character(wgcna$output_path) || nchar(wgcna$output_path) == 0) {
    stop("'wgcna$output_path' must be a non-empty string.")
  }
  if (!dir.exists(dirname(wgcna$output_path))) {
    message("Creating WGCNA output directory parent: ", dirname(wgcna$output_path))
    dir.create(dirname(wgcna$output_path), recursive = TRUE)
  }
  
  # cell_types: accept empty list or character vector
  if (!is.null(wgcna$cell_types)) {
    if (is.character(wgcna$cell_types)) {
      wgcna$cell_types <- trimws(unlist(strsplit(wgcna$cell_types, ",")))
    } else if (is.list(wgcna$cell_types)) {
      wgcna$cell_types <- unlist(wgcna$cell_types)
    }
    if (length(wgcna$cell_types) > 0 && !is.character(wgcna$cell_types)) {
      stop("'wgcna$cell_types' must be a character vector or empty list.")
    }
    cfg$wgcna$cell_types <- wgcna$cell_types  # Update normalized value
  }
  
  if (!is.character(wgcna$pro_name) || nchar(wgcna$pro_name) == 0) {
    stop("'wgcna$pro_name' must be a non-empty string.")
  }
  if (!is.character(wgcna$lnc_name)) {
    stop("'wgcna$lnc_name' must be a string (can be empty).")
  }
  
  allowed_select <- c("fraction", "variance")
  if (!wgcna$gene_select_method %in% allowed_select) {
    stop(
      "Invalid 'wgcna$gene_select_method': '", wgcna$gene_select_method, "'. ",
      "Must be one of: ", paste(allowed_select, collapse = ", ")
    )
  }

  # -----------------------------
  # Finalize
  # -----------------------------
  message("✅ LncExplore configuration loaded successfully from: ", config_path)
  return(cfg)
}

# =============================================================================
# Example usage:
# =============================================================================
# config <- load_LncExplore_config("config_LncExplore.yaml")
# print(config$location$LOG2FC_THRESH)
# =============================================================================