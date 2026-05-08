# =============================================================================
# scLncR - Configuration Loader & Validator for prelnc Module
# Description: Raw FASTQ-first, technology-aware prelnc config validation.
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

validate_prelnc_tools <- function(cfg) {
  needed <- c("hisat2", "hisat2-build", "samtools", "stringtie", "gffcompare", "gffread", "CPC2.py")
  missing <- needed[Sys.which(needed) == ""]
  if (length(missing) == 0) return(invisible(TRUE))
  if (isTRUE(cfg$dry_run)) {
    warning("Missing external tools in dry_run mode: ", paste(missing, collapse = ", "))
  } else {
    stop("Missing required external tools: ", paste(missing, collapse = ", "))
  }
}

load_prelnc_config <- function(config_path) {
  if (!requireNamespace("yaml", quietly = TRUE)) {
    stop("Package 'yaml' is required but not installed.")
  }
  if (!file.exists(config_path)) stop("Configuration file not found: ", config_path)
  cfg <- yaml::read_yaml(config_path)
  if (is.null(cfg) || !is.list(cfg)) stop("Config file is empty or malformed: ", config_path)

  # backward-compatible mapping
  cfg$samples_dirs <- cfg$samples_dirs %||% cfg$samples_dir %||% ""
  cfg$input_type <- tolower(cfg$input_type %||% "fastq")
  cfg$sequencing_platform <- tolower(cfg$sequencing_platform %||% "10x")
  cfg$read_layout <- tolower(cfg$read_layout %||% "auto")
  cfg$aligner_for_discovery <- tolower(cfg$aligner_for_discovery %||% cfg$aligner %||% "hisat2")
  cfg$assembler <- tolower(cfg$assembler %||% "stringtie")
  cfg$strandness <- cfg$strandness %||% "RF"
  cfg$min_transcript_length <- cfg$min_transcript_length %||% 200
  cfg$gffcompare_classes <- cfg$gffcompare_classes %||% c("i", "u", "x", "o")
  cfg$min_mapping_quality <- cfg$min_mapping_quality %||% 30
  cfg$rebuild_hisat2_index <- isTRUE(cfg$rebuild_hisat2_index)
  cfg$hisat2_index_prefix <- cfg$hisat2_index_prefix %||% ""
  cfg$keep_intermediate <- if (is.null(cfg$keep_intermediate)) TRUE else isTRUE(cfg$keep_intermediate)
  cfg$dry_run <- isTRUE(cfg$dry_run)
  cfg$sample_limit <- cfg$sample_limit %||% 0
  cfg$include_samples <- cfg$include_samples %||% character(0)
  cfg$samples_info <- cfg$samples_info %||% ""
  cfg$fastq_sample_regex <- cfg$fastq_sample_regex %||% "^(.*?)_S\\d+_L\\d+_([IR][12])_\\d+\\.fastq\\.gz$"
  cfg$r1_pattern <- cfg$r1_pattern %||% "_R1_"
  cfg$r2_pattern <- cfg$r2_pattern %||% "_R2_"
  cfg$i1_pattern <- cfg$i1_pattern %||% "_I1_"
  cfg$project_name <- cfg$project_name %||% "scLncR"
  cfg$threads <- cfg$threads %||% 4

  # required fields
  cfg$output_path <- as_nonempty_chr(cfg$output_path, "output_path")
  cfg$genome_file <- as_nonempty_chr(cfg$genome_file, "genome_file")
  cfg$gtf_file <- as_nonempty_chr(cfg$gtf_file, "gtf_file")
  cfg$lncrna_name <- as_nonempty_chr(cfg$lncrna_name, "lncrna_name")
  cfg$samples_dirs <- as_nonempty_chr(cfg$samples_dirs, "samples_dirs")

  if (grepl("[^a-zA-Z0-9-]", cfg$lncrna_name)) {
    stop("Invalid characters in 'lncrna_name'. Allowed: letters, digits, hyphen.")
  }

  cfg$threads <- as_int_ge(cfg$threads, "threads", min_value = 1L, strict_gt = TRUE)
  cfg$min_mapping_quality <- as_int_ge(cfg$min_mapping_quality, "min_mapping_quality", min_value = 0L)
  cfg$min_transcript_length <- as_int_ge(cfg$min_transcript_length, "min_transcript_length", min_value = 1L, strict_gt = TRUE)
  cfg$sample_limit <- as_int_ge(cfg$sample_limit, "sample_limit", min_value = 0L)

  if (!is.character(cfg$include_samples)) cfg$include_samples <- as.character(unlist(cfg$include_samples))
  cfg$include_samples <- unique(cfg$include_samples[nzchar(cfg$include_samples)])

  # value constraints
  if (cfg$input_type != "fastq") {
    warning("input_type='", cfg$input_type, "' is deprecated for prelnc; forcing input_type='fastq' in raw FASTQ-first workflow.")
    cfg$input_type <- "fastq"
  }
  cfg$sequencing_platform <- validate_choice(cfg$sequencing_platform, "sequencing_platform", c("10x", "smartseq2", "dropseq", "generic"))
  cfg$read_layout <- validate_choice(cfg$read_layout, "read_layout", c("auto", "single", "paired"))
  cfg$aligner_for_discovery <- validate_choice(cfg$aligner_for_discovery, "aligner_for_discovery", c("hisat2"))
  cfg$assembler <- validate_choice(cfg$assembler, "assembler", c("stringtie"))
  cfg$strandness <- validate_choice(cfg$strandness, "strandness", c("RF", "FR", "unstranded"))

  if (!is.character(cfg$gffcompare_classes) || length(cfg$gffcompare_classes) == 0) {
    stop("'gffcompare_classes' must be a non-empty character vector.")
  }
  cfg$gffcompare_classes <- unique(cfg$gffcompare_classes)

  # path checks
  if (!dir.exists(cfg$samples_dirs)) stop("samples_dirs does not exist: ", cfg$samples_dirs)
  if (!file.exists(cfg$genome_file)) stop("Reference genome not found: ", cfg$genome_file)
  if (!file.exists(cfg$gtf_file)) stop("Reference GTF not found: ", cfg$gtf_file)
  if (!dir.exists(cfg$output_path)) dir.create(cfg$output_path, recursive = TRUE, showWarnings = FALSE)
  if (nzchar(cfg$samples_info) && !file.exists(cfg$samples_info)) warning("samples_info not found: ", cfg$samples_info)

  # old BAM-first fields retained but not prioritized
  if (!is.null(cfg$bam_manifest) && nzchar(cfg$bam_manifest)) {
    warning("bam_manifest provided but prelnc is now raw FASTQ-first; bam_manifest will be ignored.")
  }
  if (!is.null(cfg$prelnc_mode) && tolower(cfg$prelnc_mode) == "bam_core") {
    warning("prelnc_mode='bam_core' is deprecated; raw FASTQ-first workflow will be used.")
  }

  validate_prelnc_tools(cfg)

  message("Configuration loaded successfully from: ", config_path)
  message(sprintf("input_type=%s | sequencing_platform=%s | dry_run=%s", cfg$input_type, cfg$sequencing_platform, cfg$dry_run))
  cfg
}
