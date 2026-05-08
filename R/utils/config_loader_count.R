# =============================================================================
# scLncR - Configuration Loader & Validator for count Module
# Description: Raw FASTQ-first, technology-aware count interface validation.
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
  if (!is.numeric(x) || length(x) != 1 || is.na(x)) stop(sprintf("'%s' must be numeric.", field_name))
  xi <- as.integer(x)
  if (strict_gt) {
    if (xi <= min_value) stop(sprintf("'%s' must be > %s.", field_name, min_value))
  } else {
    if (xi < min_value) stop(sprintf("'%s' must be >= %s.", field_name, min_value))
  }
  xi
}

validate_choice <- function(value, field_name, choices) {
  if (!(value %in% choices)) {
    stop(sprintf("Invalid '%s': %s. Allowed: %s", field_name, value, paste(choices, collapse = ", ")))
  }
  value
}

validate_count_tools <- function(cfg) {
  needed <- character(0)
  if (cfg$sequencing_platform == "10x" && cfg$count_engine == "cellranger") {
    needed <- c("cellranger")
  } else if (cfg$count_engine == "featurecounts") {
    needed <- c("featureCounts")
  } else if (cfg$count_engine == "starsolo") {
    needed <- c("STAR")
  }

  missing <- needed[Sys.which(needed) == ""]
  if (length(missing) == 0) return(invisible(TRUE))
  if (isTRUE(cfg$dry_run)) {
    warning("Missing external tools in dry_run mode: ", paste(missing, collapse = ", "))
  } else {
    stop("Missing required external tools: ", paste(missing, collapse = ", "))
  }
}

load_count_config <- function(config_path) {
  if (!requireNamespace("yaml", quietly = TRUE)) {
    stop("Package 'yaml' is required but not installed.")
  }
  if (!file.exists(config_path)) stop("Configuration file not found: ", config_path)
  cfg <- yaml::read_yaml(config_path)
  if (is.null(cfg) || !is.list(cfg)) stop("Config file is empty or malformed: ", config_path)

  # backward-compatible mapping from old fields
  cfg$samples_dirs <- cfg$samples_dirs %||% cfg$samples_dir %||% ""
  cfg$genome_file <- cfg$genome_file %||% cfg$genome %||% ""
  cfg$known_gtf <- cfg$known_gtf %||% cfg$gtf %||% ""
  cfg$combined_gtf <- cfg$combined_gtf %||% ""
  cfg$sequencing_platform <- tolower(cfg$sequencing_platform %||% "10x")
  cfg$count_engine <- tolower(cfg$count_engine %||% "cellranger")
  cfg$samples_info <- cfg$samples_info %||% ""
  cfg$project_name <- cfg$project_name %||% "scLncR"
  cfg$threads <- cfg$threads %||% 4
  cfg$dry_run <- isTRUE(cfg$dry_run)
  cfg$sample_limit <- cfg$sample_limit %||% 0
  cfg$include_samples <- cfg$include_samples %||% character(0)
  cfg$fastq_sample_regex <- cfg$fastq_sample_regex %||% "^(.*?)_S\\d+_L\\d+_([IR][12])_\\d+\\.fastq\\.gz$"
  cfg$r1_pattern <- cfg$r1_pattern %||% "_R1_"
  cfg$r2_pattern <- cfg$r2_pattern %||% "_R2_"
  cfg$i1_pattern <- cfg$i1_pattern %||% "_I1_"

  cfg$output_path <- as_nonempty_chr(cfg$output_path, "output_path")
  cfg$samples_dirs <- as_nonempty_chr(cfg$samples_dirs, "samples_dirs")
  cfg$genome_file <- as_nonempty_chr(cfg$genome_file, "genome_file")
  cfg$known_gtf <- as_nonempty_chr(cfg$known_gtf, "known_gtf")
  cfg$lnc_gtf <- as_nonempty_chr(cfg$lnc_gtf, "lnc_gtf")

  cfg$threads <- as_int_ge(cfg$threads, "threads", min_value = 1L, strict_gt = TRUE)
  cfg$sample_limit <- as_int_ge(cfg$sample_limit, "sample_limit", min_value = 0L)
  if (!is.character(cfg$include_samples)) cfg$include_samples <- as.character(unlist(cfg$include_samples))
  cfg$include_samples <- unique(cfg$include_samples[nzchar(cfg$include_samples)])

  cfg$sequencing_platform <- validate_choice(cfg$sequencing_platform, "sequencing_platform", c("10x", "smartseq2", "dropseq", "generic"))
  cfg$count_engine <- validate_choice(cfg$count_engine, "count_engine", c("cellranger", "featurecounts", "starsolo"))

  if (!dir.exists(cfg$samples_dirs)) stop("samples_dirs does not exist: ", cfg$samples_dirs)
  if (!file.exists(cfg$genome_file)) stop("genome_file not found: ", cfg$genome_file)
  if (!file.exists(cfg$known_gtf)) stop("known_gtf not found: ", cfg$known_gtf)
  if (!file.exists(cfg$lnc_gtf)) stop("lnc_gtf not found: ", cfg$lnc_gtf)
  if (nzchar(cfg$combined_gtf) && !file.exists(cfg$combined_gtf)) {
    warning("combined_gtf path does not exist; it will be rebuilt from known_gtf + lnc_gtf.")
    cfg$combined_gtf <- ""
  }
  if (nzchar(cfg$samples_info) && !file.exists(cfg$samples_info)) warning("samples_info not found: ", cfg$samples_info)
  if (!dir.exists(cfg$output_path)) dir.create(cfg$output_path, recursive = TRUE, showWarnings = FALSE)

  # support-level messaging
  if (!(cfg$sequencing_platform == "10x" && cfg$count_engine == "cellranger")) {
    warning("Current strong support is 10x + cellranger. Selected interface is experimental.")
  }

  validate_count_tools(cfg)
  message("Configuration loaded successfully from: ", config_path)
  message(sprintf("sequencing_platform=%s | count_engine=%s | dry_run=%s", cfg$sequencing_platform, cfg$count_engine, cfg$dry_run))
  cfg
}
