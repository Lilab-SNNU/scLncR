# =============================================================================
# scLncR - Configuration Loader & Validator for prelnc Module
# Author: [Yin SW]
# Description: Safely loads and validates the YAML config file for the prelnc step.
# Dependencies: yaml, assertthat (optional but recommended)
# =============================================================================


load_prelnc_config <- function(config_path) {
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
    "samples_dirs", "project_name", "output_path",
    "genome_file", "gtf_file", "lncrna_name", "threads"
  )

  missing_fields <- setdiff(required_fields, names(cfg))
  if (length(missing_fields) > 0) {
    stop("Missing required configuration fields: ", paste(missing_fields, collapse = ", "))
  }

  # -----------------------------
  # Path existence validation
  # -----------------------------

  # samples_dirs must be a directory
  if (!dir.exists(cfg$samples_dirs)) {
    warning("samples_dirs does not exist or is not a directory: ", cfg$samples_dirs)
  }

  # genome_file and gtf_file must be readable files
  if (!file.exists(cfg$genome_file)) {
    stop("Reference genome file not found: ", cfg$genome_file)
  }
  if (!file.exists(cfg$gtf_file)) {
    stop("GTF annotation file not found: ", cfg$gtf_file)
  }

  # Ensure output_path is a directory (create if needed?)
  if (!dir.exists(cfg$output_path)) {
    message("Output directory does not exist. Creating: ", cfg$output_path)
    dir.create(cfg$output_path, recursive = TRUE)
  }

  # -----------------------------
  # lncrna_name naming convention check
  # -----------------------------
  lnc_name <- cfg$lncrna_name
  if (!is.character(lnc_name) || nchar(lnc_name) == 0) {
    stop("'lncrna_name' must be a non-empty character string.")
  }

  # Disallow underscores, spaces, or special characters (only letters, digits, hyphens allowed)
  if (grepl("[^a-zA-Z0-9-]", lnc_name)) {
    stop(
      "Invalid characters in 'lncrna_name': '", lnc_name, "'. ",
      "Only letters, digits, and hyphens ('-') are allowed. ",
      "Underscores ('_') and other symbols will cause Seurat object mismatches."
    )
  }

  # -----------------------------
  # threads validation
  # -----------------------------
  threads <- cfg$threads
  if (!is.numeric(threads) || threads <= 0 || threads != as.integer(threads)) {
    stop("'threads' must be a positive integer.")
  }
  cfg$threads <- as.integer(threads)

  # -----------------------------
  # Finalize
  # -----------------------------
  message("✅ Configuration loaded successfully from: ", config_path)
  return(cfg)
}

# =============================================================================
# Example usage:
# =============================================================================
# config <- load_prelnc_config("path/to/prelnc_config.yaml")
# print(config$project_name)
# =============================================================================

