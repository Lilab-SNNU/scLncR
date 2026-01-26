# =============================================================================
# scLncR - Configuration Loader & Validator for count Module
# Author: [Yin SW]
# Description: Safely loads and validates the YAML config file for the count step.
# Dependencies: yaml
# =============================================================================


load_count_config <- function(config_path) {
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
    "samples_dir", "project_name", "output_path",
    "genome", "gtf", "lnc_gtf"
  )

  missing_fields <- setdiff(required_fields, names(cfg))
  if (length(missing_fields) > 0) {
    stop("Missing required configuration fields: ", paste(missing_fields, collapse = ", "))
  }

  # -----------------------------
  # Path existence validation
  # -----------------------------

  # samples_dir must be a directory
  if (!dir.exists(cfg$samples_dir)) {
    stop("samples_dir does not exist or is not a directory: ", cfg$samples_dir)
  }

  # genome, gtf, lnc_gtf must be readable files
  if (!file.exists(cfg$genome)) {
    stop("Reference genome file not found: ", cfg$genome)
  }
  if (!file.exists(cfg$gtf)) {
    stop("GTF annotation file not found: ", cfg$gtf)
  }
  if (!file.exists(cfg$lnc_gtf)) {
    stop("lncRNA GTF file not found: ", cfg$lnc_gtf)
  }

  # Ensure output_path is a directory (create if needed?)
  if (!dir.exists(cfg$output_path)) {
    message("Output directory does not exist. Creating: ", cfg$output_path)
    dir.create(cfg$output_path, recursive = TRUE)
  }

  # -----------------------------
  # project_name naming convention check
  # -----------------------------
  proj_name <- cfg$project_name
  if (!is.character(proj_name) || nchar(proj_name) == 0) {
    stop("'project_name' must be a non-empty character string.")
  }

  # Disallow spaces or special characters that may break downstream paths
  if (grepl("[^a-zA-Z0-9._-]", proj_name)) {
    warning(
      "Project name '", proj_name, "' contains unusual characters. ",
      "Only letters, digits, dots ('.'), underscores ('_'), and hyphens ('-') are recommended."
    )
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
# config <- load_count_config("path/to/count_config.yaml")
# print(config$project_name)
# =============================================================================
