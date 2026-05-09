# =============================================================================
# scLncR - Configuration Loader & Validator for raw FASTQ QC
# =============================================================================

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0) return(y)
  x
}

as_nonempty_chr <- function(x, field_name) {
  if (is.null(x) || !is.character(x) || length(x) != 1 || !nzchar(trimws(x))) {
    stop(sprintf("'%s' must be a non-empty string.", field_name))
  }
  trimws(x)
}

as_int_ge <- function(x, field_name, min_value = 0L, strict_gt = FALSE) {
  if (!is.numeric(x) || length(x) != 1 || is.na(x)) {
    stop(sprintf("'%s' must be numeric.", field_name))
  }
  xi <- as.integer(x)
  if (strict_gt && xi <= min_value) stop(sprintf("'%s' must be > %s.", field_name, min_value))
  if (!strict_gt && xi < min_value) stop(sprintf("'%s' must be >= %s.", field_name, min_value))
  xi
}

as_bool <- function(x, default = FALSE) {
  if (is.null(x)) return(default)
  isTRUE(x)
}

normalize_string_vector <- function(x) {
  if (is.null(x)) return(character(0))
  if (!is.character(x)) x <- as.character(unlist(x))
  unique(trimws(x[nzchar(trimws(x))]))
}

#' Load and validate raw FASTQ QC configuration
#'
#' @description
#' Reads a YAML configuration file for the optional raw FASTQ QC module.
#' The module runs FastQC and MultiQC only; it does not trim, filter, or
#' modify input FASTQ files.
#'
#' @param config_file Path to config_QC.yaml.
#' @return A validated QC configuration list.
#' @export
load_qc_config <- function(config_file) {
  if (!requireNamespace("yaml", quietly = TRUE)) {
    stop("Package 'yaml' is required but not installed.")
  }
  if (!file.exists(config_file)) stop("QC config file not found: ", config_file)
  cfg <- yaml::read_yaml(config_file)
  if (is.null(cfg) || !is.list(cfg)) stop("QC config file is empty or malformed: ", config_file)

  cfg$input_dir <- as_nonempty_chr(cfg$input_dir %||% "", "input_dir")
  cfg$output_dir <- as_nonempty_chr(cfg$output_dir %||% "", "output_dir")
  cfg$file_pattern <- as_nonempty_chr(cfg$file_pattern %||% "*.fastq.gz", "file_pattern")
  cfg$recursive <- as_bool(cfg$recursive, default = FALSE)
  cfg$threads <- as_int_ge(cfg$threads %||% 1, "threads", min_value = 1L, strict_gt = TRUE)
  cfg$dry_run <- as_bool(cfg$dry_run, default = FALSE)
  cfg$sample_limit <- as_int_ge(cfg$sample_limit %||% 0, "sample_limit", min_value = 0L)
  cfg$include_files <- normalize_string_vector(cfg$include_files)

  cfg$fastqc <- cfg$fastqc %||% list()
  cfg$fastqc$enabled <- as_bool(cfg$fastqc$enabled, default = TRUE)
  cfg$fastqc$extract <- as_bool(cfg$fastqc$extract, default = TRUE)
  cfg$fastqc$nogroup <- as_bool(cfg$fastqc$nogroup, default = FALSE)
  cfg$fastqc$contaminants <- cfg$fastqc$contaminants %||% ""
  cfg$fastqc$adapters <- cfg$fastqc$adapters %||% ""
  cfg$fastqc$extra_args <- cfg$fastqc$extra_args %||% ""

  cfg$multiqc <- cfg$multiqc %||% list()
  cfg$multiqc$enabled <- as_bool(cfg$multiqc$enabled, default = TRUE)
  cfg$multiqc$title <- cfg$multiqc$title %||% "scLncR raw FASTQ QC report"
  cfg$multiqc$filename <- cfg$multiqc$filename %||% "multiqc_report.html"
  cfg$multiqc$extra_args <- cfg$multiqc$extra_args %||% ""

  if (!dir.exists(cfg$input_dir)) stop("input_dir does not exist: ", cfg$input_dir)
  if (nzchar(cfg$fastqc$contaminants) && !file.exists(cfg$fastqc$contaminants)) {
    stop("fastqc.contaminants file not found: ", cfg$fastqc$contaminants)
  }
  if (nzchar(cfg$fastqc$adapters) && !file.exists(cfg$fastqc$adapters)) {
    stop("fastqc.adapters file not found: ", cfg$fastqc$adapters)
  }
  if (!dir.exists(cfg$output_dir)) dir.create(cfg$output_dir, recursive = TRUE, showWarnings = FALSE)

  cfg
}

#' Validate FastQC and MultiQC command availability
#'
#' @param cfg A validated QC configuration list.
#' @return Invisibly returns TRUE when required tools are available.
#' @export
validate_qc_tools <- function(cfg) {
  needed <- character(0)
  if (isTRUE(cfg$fastqc$enabled)) needed <- c(needed, "fastqc")
  if (isTRUE(cfg$multiqc$enabled)) needed <- c(needed, "multiqc")
  missing <- needed[Sys.which(needed) == ""]
  if (length(missing) == 0) return(invisible(TRUE))

  msg <- paste("Missing QC tools:", paste(missing, collapse = ", "))
  if (isTRUE(cfg$dry_run)) {
    warning(msg, ". dry_run=true, commands will only be written.")
  } else {
    stop(msg)
  }
  invisible(FALSE)
}
