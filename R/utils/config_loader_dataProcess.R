# =============================================================================
# scLncR - Configuration Loader & Validator for dataProcess Module
# Author: [Yin SW]
# Description: Safely loads and validates the YAML config file for the dataProcess step.
# Dependencies: yaml
# =============================================================================


load_dataProcess_config <- function(config_path) {
  # Load required packages
  if (!requireNamespace("yaml", quietly = TRUE)) {
    stop("Package 'yaml' is required. Install it in the active scLncR R/conda environment.")
  }

  # Check if config file exists
  if (!file.exists(config_path)) {
    stop("Configuration file not found: ", config_path)
  }

  # Read YAML
  cfg <- yaml::read_yaml(config_path)
  if (is.null(cfg) || !is.list(cfg)) {
    stop("Configuration file is empty or malformed: ", config_path)
  }

  `%||%` <- function(x, y) {
    if (is.null(x) || length(x) == 0) return(y)
    x
  }

  cfg$input_format <- tolower(cfg$input_format %||% "10x")
  cfg$sequencing_platform <- tolower(cfg$sequencing_platform %||% "10x")
  cfg$counts_dir <- cfg$counts_dir %||% ""
  cfg$counts_matrix <- cfg$counts_matrix %||% ""
  cfg$project_name <- cfg$project_name %||% "scLncR"
  cfg$samples_info <- cfg$samples_info %||% ""

  # -----------------------------
  # Required fields check
  # -----------------------------
  required_fields <- c(
    "mt_name", "min.RNAs", "max.RNAs", "percent.mt",
    "lnc_name", "nfeatures", "dims", "resolution",
    "samples_info", "anno_method", "ref_file", "marker_gene_file", 
    "tissue", "n_top_markers", "output_dir"
  )

  missing_fields <- setdiff(required_fields, names(cfg))
  if (length(missing_fields) > 0) {
    stop("Missing required configuration fields: ", paste(missing_fields, collapse = ", "))
  }

  # -----------------------------
  # Path existence validation
  # -----------------------------

  allowed_input_formats <- c("10x", "featurecounts_matrix")
  if (!(cfg$input_format %in% allowed_input_formats)) {
    stop(
      "Invalid 'input_format': '", cfg$input_format, "'. ",
      "Must be one of: ", paste(allowed_input_formats, collapse = ", ")
    )
  }

  allowed_platforms <- c("10x", "smartseq2", "generic")
  if (!(cfg$sequencing_platform %in% allowed_platforms)) {
    stop(
      "Invalid 'sequencing_platform': '", cfg$sequencing_platform, "'. ",
      "Must be one of: ", paste(allowed_platforms, collapse = ", ")
    )
  }

  if (cfg$input_format == "10x") {
    if (!is.character(cfg$counts_dir) || !nzchar(cfg$counts_dir)) {
      stop("'counts_dir' is required when input_format='10x'.")
    }
    if (!dir.exists(cfg$counts_dir)) {
      stop("Directory in 'counts_dir' does not exist: ", cfg$counts_dir)
    }
  }

  if (cfg$input_format == "featurecounts_matrix") {
    if (!is.character(cfg$counts_matrix) || !nzchar(cfg$counts_matrix)) {
      stop("'counts_matrix' is required when input_format='featurecounts_matrix'.")
    }
    if (!file.exists(cfg$counts_matrix)) {
      stop("featureCounts matrix file not found: ", cfg$counts_matrix)
    }
    if (cfg$sequencing_platform != "smartseq2") {
      warning("input_format='featurecounts_matrix' is primarily intended for sequencing_platform='smartseq2'.")
    }
  }
  
  # samples_info is optional, but must be readable when provided.
  if (nzchar(cfg$samples_info) && !file.exists(cfg$samples_info)) {
    stop("Sample metadata file not found: ", cfg$samples_info)
  }

  # ref_file and marker_gene_file: conditional validation (see below)

  # Ensure output_dir is a directory (create if needed)
  if (!dir.exists(cfg$output_dir)) {
    message("Output directory does not exist. Creating: ", cfg$output_dir)
    dir.create(cfg$output_dir, recursive = TRUE)
  }

  # -----------------------------
  # String field validations
  # -----------------------------

  # mt_name
  if (!is.character(cfg$mt_name) || nchar(cfg$mt_name) == 0) {
    stop("'mt_name' must be a non-empty character string (e.g., 'MT').")
  }

  # lnc_name (can be empty, but must be character)
  if (!is.character(cfg$lnc_name)) {
    stop("'lnc_name' must be a character string (can be empty).")
  }

  # anno_method
  allowed_methods <- c("SingleR", "scMM", "none")
  if (tolower(cfg$anno_method) == "none") cfg$anno_method <- "none"
  if (!cfg$anno_method %in% allowed_methods) {
    stop(
      "Invalid 'anno_method': '", cfg$anno_method, "'. ",
      "Must be one of: ", paste(allowed_methods, collapse = ", ")
    )
  }

  # Conditional: ref_file required for SingleR
  if (cfg$anno_method == "SingleR") {
    if (!is.character(cfg$ref_file) || nchar(cfg$ref_file) == 0) {
      stop("'ref_file' is required when 'anno_method' is 'SingleR'.")
    }
    if (!file.exists(cfg$ref_file)) {
      stop("Reference file for SingleR not found: ", cfg$ref_file)
    }
  }

  # Conditional: marker_gene_file and tissue required for scMM
  if (cfg$anno_method == "scMM") {
    if (!is.character(cfg$marker_gene_file) || nchar(cfg$marker_gene_file) == 0) {
      stop("'marker_gene_file' is required when 'anno_method' is 'scMM'.")
    }
    if (!file.exists(cfg$marker_gene_file)) {
      stop("Marker gene file for scMM not found: ", cfg$marker_gene_file)
    }

    if (!is.character(cfg$tissue) || nchar(cfg$tissue) == 0) {
      stop("'tissue' is required when 'anno_method' is 'scMM'.")
    }
  }

  # colour (optional, but if provided, must be character)
  if (!is.null(cfg$colour) && !is.character(cfg$colour)) {
    stop("'colour' must be a character string or empty.")
  }

  # -----------------------------
  # Numeric parameter validations
  # -----------------------------

  # Helper: validate positive integer or numeric
  validate_positive_num <- function(x, name) {
    if (!is.numeric(x) || length(x) != 1 || is.na(x) || x <= 0) {
      stop("'", name, "' must be a positive number.")
    }
  }

  validate_nonneg_int <- function(x, name) {
    if (!is.numeric(x) || length(x) != 1 || is.na(x) || x < 0 || x != as.integer(x)) {
      stop("'", name, "' must be a non-negative integer.")
    }
  }

  validate_positive_num(cfg$min.RNAs, "min.RNAs")
  validate_positive_num(cfg$max.RNAs, "max.RNAs")
  if (cfg$min.RNAs >= cfg$max.RNAs) {
    stop("'min.RNAs' must be less than 'max.RNAs'.")
  }

  validate_positive_num(cfg$percent.mt, "percent.mt")
  if (cfg$percent.mt > 100) {
    warning("'percent.mt' is greater than 100; this will filter out all cells.")
  }

  validate_positive_num(cfg$nfeatures, "nfeatures")
  validate_positive_num(cfg$dims, "dims")
  validate_positive_num(cfg$resolution, "resolution")

  validate_nonneg_int(cfg$n_top_markers, "n_top_markers")
  if (cfg$n_top_markers == 0) {
    warning("'n_top_markers' is 0; annotation may lack marker support.")
  }

  # -----------------------------
  # Finalize
  # -----------------------------
  message("✅ Configuration loaded successfully from: ", config_path)
  return(cfg)
}

# =============================================================================
# Example usage:
# =============================================================================
# config <- load_dataProcess_config("path/to/dataProcess_config.yaml")
# print(config$output_dir)
# =============================================================================
