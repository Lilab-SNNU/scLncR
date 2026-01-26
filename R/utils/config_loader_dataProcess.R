# =============================================================================
# scLncR - Configuration Loader & Validator for dataProcess Module
# Author: [Yin SW]
# Description: Safely loads and validates the YAML config file for the dataProcess step.
# Dependencies: yaml
# =============================================================================


load_dataProcess_config <- function(config_path) {
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
  # Required fields check
  # -----------------------------
  required_fields <- c(
    "counts_dir", "mt_name", "min.RNAs", "max.RNAs", "percent.mt",
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

  # counts_dir: comma-separated paths; each must be a directory
  # Ensure output_dir is a directory (create if needed)
  if (!dir.exists(cfg$counts_dir)) {
    stop("Directory in 'counts_dir' does not exist: ", cfg$counts_dir)
  }
  
  # samples_info must be a readable file
  if (!file.exists(cfg$samples_info)) {
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
  allowed_methods <- c("SingleR", "scMM")
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
