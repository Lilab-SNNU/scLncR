suppressPackageStartupMessages({
  suppressWarnings({
    suppressMessages({
      library(optparse)
      library(yaml)
    })
  })
})

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0) return(y)
  x
}

now_string <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S")

append_line <- function(path, txt) {
  cat(txt, file = path, append = TRUE, sep = "\n")
}

sanitize_sample_id <- function(x) gsub("[^A-Za-z0-9._-]", "-", x)

log_count <- function(ctx, msg, level = "INFO") {
  line <- sprintf("[%s] [%s] %s", now_string(), level, msg)
  append_line(ctx$files$run_log, line)
  message(line)
}

add_count_warning <- function(ctx, msg) {
  ctx$warnings <- unique(c(ctx$warnings, msg))
  log_count(ctx, msg, level = "WARN")
}

record_count_command <- function(ctx, step, command, args = character(0)) {
  cmd <- if (length(args) == 0) command else paste(c(command, shQuote(args)), collapse = " ")
  append_line(ctx$files$commands_log, sprintf("[%s] [%s] %s", now_string(), step, cmd))
  ctx$command_preview <- c(ctx$command_preview, sprintf("[%s] %s", step, cmd))
  cmd
}

run_count_command <- function(ctx, step, command, args = character(0), check = TRUE, workdir = NULL) {
  cmd <- record_count_command(ctx, step, command, args)
  if (isTRUE(ctx$cfg$dry_run)) {
    log_count(ctx, sprintf("[dry_run] %s", cmd))
    return(list(ok = TRUE, status = NA_integer_, output = character(0)))
  }

  out_tmp <- tempfile("count_cmd_")
  res <- tryCatch(
    {
      oldwd <- NULL
      if (!is.null(workdir)) {
        dir.create(workdir, recursive = TRUE, showWarnings = FALSE)
        oldwd <- getwd()
        setwd(workdir)
        on.exit(setwd(oldwd), add = TRUE)
      }
      status <- suppressWarnings(system2(command, args = shQuote(args), stdout = out_tmp, stderr = out_tmp))
      list(status = as.integer(status), err = NULL)
    },
    error = function(e) list(status = 127L, err = e$message)
  )
  out <- character(0)
  if (file.exists(out_tmp)) {
    out <- readLines(out_tmp, warn = FALSE)
    unlink(out_tmp)
  }
  if (!is.null(res$err)) out <- c(out, sprintf("R ERROR: %s", res$err))
  if (length(out) > 0) append_line(ctx$files$run_log, paste(out, collapse = "\n"))

  ok <- identical(res$status, 0L)
  if (!ok) {
    log_count(ctx, sprintf("Command failed [%s] exit=%s", step, res$status), level = "ERROR")
    if (check) stop(sprintf("Command failed at step '%s': %s", step, cmd))
  }
  list(ok = ok, status = res$status, output = out)
}

init_count_context <- function(cfg) {
  dirs <- list(
    base = cfg$output_path,
    manifest = file.path(cfg$output_path, "manifest"),
    logs = file.path(cfg$output_path, "logs"),
    reference = file.path(cfg$output_path, "reference"),
    index = file.path(cfg$output_path, "index"),
    hisat2_index = file.path(cfg$output_path, "index", "hisat2"),
    alignment = file.path(cfg$output_path, "alignment"),
    sam = file.path(cfg$output_path, "alignment", "sam"),
    sorted_bam = file.path(cfg$output_path, "alignment", "sorted_bam"),
    featurecounts = file.path(cfg$output_path, "featurecounts"),
    count_results = file.path(cfg$output_path, "count_results")
  )
  for (d in unlist(dirs, use.names = FALSE)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

  ctx <- new.env(parent = emptyenv())
  ctx$cfg <- cfg
  ctx$dirs <- dirs
  ctx$files <- list(
    commands_log = file.path(dirs$logs, "commands.log"),
    run_log = file.path(dirs$logs, "count_run.log"),
    status = file.path(dirs$logs, "per_sample_status.tsv"),
    fastq_manifest = file.path(dirs$manifest, "count_raw_fastq_manifest.tsv"),
    smartseq2_manifest = file.path(dirs$manifest, "smartseq2_fastq_manifest.tsv"),
    validation_report = file.path(dirs$manifest, "count_input_validation_report.md"),
    report = file.path(dirs$base, "count_run_report.md"),
    combined_gtf = file.path(dirs$reference, "combined_mRNA_lncRNA.gtf"),
    featurecounts_output = file.path(dirs$featurecounts, "smartseq2_featureCounts.txt"),
    count_matrix_tsv = file.path(dirs$featurecounts, "smartseq2_count_matrix.tsv"),
    count_matrix_rds = file.path(dirs$featurecounts, "smartseq2_count_matrix.rds")
  )
  ctx$warnings <- character(0)
  ctx$command_preview <- character(0)
  ctx$status_rows <- data.frame(
    step = character(0),
    sample_id = character(0),
    status = character(0),
    message = character(0),
    stringsAsFactors = FALSE
  )
  writeLines(character(0), ctx$files$commands_log)
  writeLines(character(0), ctx$files$run_log)
  write.table(ctx$status_rows, ctx$files$status, sep = "\t", row.names = FALSE, quote = FALSE)
  ctx
}

add_count_status <- function(ctx, step, sample_id = "", status = "OK", message = "") {
  ctx$status_rows <- rbind(
    ctx$status_rows,
    data.frame(step = step, sample_id = sample_id, status = status, message = message, stringsAsFactors = FALSE)
  )
  write.table(ctx$status_rows, ctx$files$status, sep = "\t", row.names = FALSE, quote = FALSE)
}

parse_count_fastq_identity <- function(filename, cfg) {
  pattern <- cfg$fastq_sample_regex
  m <- regexec(pattern, filename, perl = TRUE, ignore.case = TRUE)
  reg <- regmatches(filename, m)[[1]]
  if (length(reg) >= 3) return(list(sample_id = reg[2], token = toupper(reg[3])))
  token <- "SINGLE"
  if (grepl(cfg$r2_pattern, filename, fixed = TRUE)) token <- "R2"
  if (grepl(cfg$r1_pattern, filename, fixed = TRUE)) token <- "R1"
  if (grepl(cfg$i1_pattern, filename, fixed = TRUE)) token <- "I1"
  sid <- sub("\\.(fastq|fq)(\\.gz)?$", "", filename, ignore.case = TRUE)
  sid <- sub("(_S\\d+_L\\d+_[IR][12]_\\d+)$", "", sid, perl = TRUE, ignore.case = TRUE)
  sid <- sub("(_R1_|_R2_|_I1_|_I2_).*$", "", sid, perl = TRUE, ignore.case = TRUE)
  list(sample_id = sid, token = token)
}

discover_count_fastq_manifest <- function(cfg, ctx) {
  all_files <- list.files(cfg$samples_dirs, recursive = TRUE, full.names = TRUE)
  fastq_files <- all_files[grepl("\\.(fastq|fq)(\\.gz)?$", all_files, ignore.case = TRUE)]
  if (length(fastq_files) == 0) stop("No FASTQ files found in samples_dirs: ", cfg$samples_dirs)

  rows <- lapply(fastq_files, function(f) {
    id <- parse_count_fastq_identity(basename(f), cfg)
    data.frame(sample_id = id$sample_id, token = id$token, fastq_path = normalizePath(f, winslash = "/", mustWork = FALSE), stringsAsFactors = FALSE)
  })
  df <- do.call(rbind, rows)
  sids <- sort(unique(df$sample_id))
  if (length(cfg$include_samples) > 0) sids <- intersect(sids, cfg$include_samples)
  if (cfg$sample_limit > 0) sids <- sids[seq_len(min(cfg$sample_limit, length(sids)))]
  if (length(sids) == 0) stop("No count samples left after filtering.")

  out <- data.frame(
    sample_id = sids,
    sequencing_platform = cfg$sequencing_platform,
    R1_path = NA_character_,
    R2_path = NA_character_,
    I1_path = NA_character_,
    single_fastq_path = NA_character_,
    stringsAsFactors = FALSE
  )
  for (i in seq_len(nrow(out))) {
    sub <- df[df$sample_id == out$sample_id[i], , drop = FALSE]
    pick <- function(tok) {
      vals <- sub$fastq_path[sub$token == tok]
      if (length(vals) == 0) NA_character_ else vals[1]
    }
    out$R1_path[i] <- pick("R1")
    out$R2_path[i] <- pick("R2")
    out$I1_path[i] <- pick("I1")
    out$single_fastq_path[i] <- pick("SINGLE")
  }
  write.table(out, ctx$files$fastq_manifest, sep = "\t", row.names = FALSE, quote = FALSE)
  out
}

build_combined_gtf_for_count <- function(cfg, ctx) {
  if (!is.null(cfg$combined_gtf) && nzchar(cfg$combined_gtf) && file.exists(cfg$combined_gtf)) {
    file.copy(cfg$combined_gtf, ctx$files$combined_gtf, overwrite = TRUE)
    return(ctx$files$combined_gtf)
  }
  if (!file.exists(cfg$known_gtf)) stop("known_gtf not found: ", cfg$known_gtf)
  if (!file.exists(cfg$lnc_gtf)) stop("lnc_gtf not found: ", cfg$lnc_gtf)

  writeLines(c(readLines(cfg$known_gtf, warn = FALSE), readLines(cfg$lnc_gtf, warn = FALSE)), ctx$files$combined_gtf)
  ctx$files$combined_gtf
}

build_cellranger_reference <- function(cfg, ctx, combined_gtf) {
  ref_name <- paste0(cfg$project_name, "_lncRef")
  ref_out <- file.path(ctx$dirs$reference, ref_name)
  if (dir.exists(ref_out)) {
    if (dir.exists(file.path(ref_out, "star")) && file.exists(file.path(ref_out, "genes", "genes.gtf"))) {
      log_count(ctx, sprintf("Reusing existing Cell Ranger reference: %s", ref_out))
      return(ref_out)
    }
    if (!isTRUE(cfg$dry_run)) {
      stop(
        "Cell Ranger reference output directory exists but does not look complete: ", ref_out,
        "\nRemove it or choose a new output_path before rerunning."
      )
    }
  }
  args <- c(
    "mkref",
    "--fasta", cfg$genome_file,
    "--genes", combined_gtf,
    "--genome", ref_name,
    "--nthreads", as.character(cfg$threads),
    "--output-dir", ref_out
  )
  run_count_command(ctx, "cellranger_mkref", "cellranger", args, check = TRUE)
  ref_out
}

run_cellranger_count_for_sample <- function(sample_id, cfg, ctx, transcriptome_path) {
  out_id <- sanitize_sample_id(sample_id)
  out_dir <- file.path(ctx$dirs$count_results, out_id)
  if (!isTRUE(cfg$dry_run) && dir.exists(out_dir)) {
    stop(
      "Cell Ranger count output already exists: ", out_dir,
      "\nRemove it or choose a new output_path before rerunning this sample."
    )
  }
  args <- c(
    "count",
    paste0("--id=", out_id),
    paste0("--fastqs=", cfg$samples_dirs),
    paste0("--sample=", sample_id),
    paste0("--transcriptome=", transcriptome_path),
    paste0("--localcores=", cfg$threads)
  )
  run_count_command(ctx, paste0("cellranger_count_", out_id), "cellranger", args, check = TRUE, workdir = ctx$dirs$count_results)
}

#' Discover Smart-seq2 FASTQ inputs for count analysis
#'
#' @description
#' Scans a raw FASTQ directory and builds a Smart-seq2 sample manifest for
#' paired-end or single-end gene-level quantification. Multiple lanes for the
#' same sample are retained as comma-separated FASTQ lists for HISAT2.
#'
#' @param cfg A validated count configuration list.
#' @param ctx Count workflow context containing output paths and logs.
#' @return A data.frame with sample IDs, FASTQ paths, detected read layout,
#' sequencing platform, and notes.
#' @export
discover_smartseq2_fastq_inputs <- function(cfg, ctx) {
  fastq_files <- list.files(cfg$samples_dirs, recursive = TRUE, full.names = TRUE)
  fastq_files <- fastq_files[grepl("\\.(fastq|fq)(\\.gz)?$", fastq_files, ignore.case = TRUE)]
  if (length(fastq_files) == 0) stop("No FASTQ files found in samples_dirs: ", cfg$samples_dirs)

  strip_ext <- function(x) sub("\\.(fastq|fq)(\\.gz)?$", "", x, ignore.case = TRUE)
  infer_sample <- function(filename, token) {
    x <- strip_ext(filename)
    if (token %in% c("R1", "R2")) {
      x <- sub("(_S\\d+)?(_L\\d+)?_R[12](_\\d+)?$", "", x, perl = TRUE, ignore.case = TRUE)
      x <- sub("([._-])R[12]$", "", x, perl = TRUE, ignore.case = TRUE)
      x <- sub("([._-])[12]$", "", x, perl = TRUE, ignore.case = TRUE)
    }
    x
  }
  token_for <- function(filename) {
    if (nzchar(cfg$r1_pattern) && grepl(cfg$r1_pattern, filename, fixed = TRUE)) return("R1")
    if (nzchar(cfg$r2_pattern) && grepl(cfg$r2_pattern, filename, fixed = TRUE)) return("R2")
    if (nzchar(cfg$single_pattern) && !grepl(cfg$single_pattern, filename, fixed = TRUE)) return("UNASSIGNED")
    "SINGLE"
  }

  rows <- lapply(sort(fastq_files), function(f) {
    tok <- token_for(basename(f))
    data.frame(
      sample_id = infer_sample(basename(f), tok),
      token = tok,
      fastq_path = normalizePath(f, winslash = "/", mustWork = FALSE),
      stringsAsFactors = FALSE
    )
  })
  df <- do.call(rbind, rows)
  df <- df[df$token != "UNASSIGNED", , drop = FALSE]
  if (nrow(df) == 0) stop("No Smart-seq2 FASTQ files matched r1_pattern/r2_pattern/single_pattern.")
  sids <- sort(unique(df$sample_id))
  if (length(cfg$include_samples) > 0) sids <- intersect(sids, cfg$include_samples)
  if (cfg$sample_limit > 0) sids <- sids[seq_len(min(cfg$sample_limit, length(sids)))]
  if (length(sids) == 0) stop("No Smart-seq2 samples left after filtering.")

  collapse_paths <- function(x) {
    x <- sort(unique(x[!is.na(x) & nzchar(x)]))
    if (length(x) == 0) NA_character_ else paste(x, collapse = ",")
  }
  out <- data.frame(
    sample_id = sids,
    R1_path = NA_character_,
    R2_path = NA_character_,
    single_fastq_path = NA_character_,
    read_layout_detected = NA_character_,
    sequencing_platform = cfg$sequencing_platform,
    notes = NA_character_,
    stringsAsFactors = FALSE
  )
  for (i in seq_len(nrow(out))) {
    sub <- df[df$sample_id == out$sample_id[i], , drop = FALSE]
    r1 <- collapse_paths(sub$fastq_path[sub$token == "R1"])
    r2 <- collapse_paths(sub$fastq_path[sub$token == "R2"])
    se <- collapse_paths(sub$fastq_path[sub$token == "SINGLE"])
    out$R1_path[i] <- r1
    out$R2_path[i] <- r2
    out$single_fastq_path[i] <- se
    detected <- if (!is.na(r1) || !is.na(r2)) "paired" else "single"
    if (cfg$read_layout != "auto") detected <- cfg$read_layout
    out$read_layout_detected[i] <- detected
    out$notes[i] <- "Smart-seq2 count uses HISAT2 alignment and featureCounts gene-level quantification."
  }
  write.table(out, ctx$files$smartseq2_manifest, sep = "\t", row.names = FALSE, quote = FALSE)
  out
}

#' Validate a Smart-seq2 count FASTQ manifest
#'
#' @description
#' Checks paired-end or single-end FASTQ requirements and writes a compact input
#' validation report. Mixed paired/single layouts are rejected to keep the
#' featureCounts matrix generation deterministic.
#'
#' @param manifest A Smart-seq2 FASTQ manifest from \code{discover_smartseq2_fastq_inputs()}.
#' @param cfg A validated count configuration list.
#' @param ctx Count workflow context containing output paths and logs.
#' @return Invisibly returns TRUE when validation succeeds.
#' @export
validate_smartseq2_manifest <- function(manifest, cfg, ctx) {
  if (anyDuplicated(manifest$sample_id)) stop("Smart-seq2 sample_id values must be unique.")
  layouts <- unique(manifest$read_layout_detected)
  if (length(layouts) != 1) stop("Mixed paired/single Smart-seq2 layouts are not supported in one featureCounts run.")
  layout <- layouts[1]
  if (layout == "paired") {
    bad <- manifest$sample_id[is.na(manifest$R1_path) | is.na(manifest$R2_path)]
    if (length(bad) > 0) stop("Paired Smart-seq2 samples missing R1/R2 FASTQ: ", paste(bad, collapse = ", "))
  } else if (layout == "single") {
    bad <- manifest$sample_id[is.na(manifest$single_fastq_path)]
    if (length(bad) > 0) stop("Single-end Smart-seq2 samples missing FASTQ: ", paste(bad, collapse = ", "))
  } else {
    stop("Invalid Smart-seq2 read layout: ", layout)
  }

  lines <- c(
    "# scLncR Smart-seq2 Count Input Validation",
    "",
    sprintf("- Time: %s", now_string()),
    sprintf("- samples_dirs: `%s`", cfg$samples_dirs),
    sprintf("- samples: %d", nrow(manifest)),
    sprintf("- read_layout: %s", layout),
    sprintf("- strandness: %s", cfg$strandness),
    "",
    "Validation result: PASS",
    "",
    "Smart-seq2 count uses HISAT2 and featureCounts. It does not use Cell Ranger."
  )
  writeLines(lines, ctx$files$validation_report)
  invisible(TRUE)
}

hisat2_index_files_exist <- function(prefix) {
  if (!nzchar(prefix)) return(FALSE)
  ht2 <- paste0(prefix, ".", 1:8, ".ht2")
  ht2l <- paste0(prefix, ".", 1:8, ".ht2l")
  all(file.exists(ht2)) || all(file.exists(ht2l))
}

#' Build a HISAT2 index for Smart-seq2 count if needed
#'
#' @description
#' Reuses an existing HISAT2 index when available, otherwise runs
#' \code{hisat2-build} against the configured genome FASTA. Commands are written
#' to the count commands log and skipped in dry-run mode.
#'
#' @param cfg A validated count configuration list.
#' @param ctx Count workflow context containing output paths and logs.
#' @return The HISAT2 index prefix used for alignment.
#' @export
build_hisat2_index_for_count_if_needed <- function(cfg, ctx) {
  prefix <- cfg$hisat2_index_prefix
  if (nzchar(prefix) && hisat2_index_files_exist(prefix) && !isTRUE(cfg$rebuild_hisat2_index)) {
    log_count(ctx, sprintf("Using existing HISAT2 index: %s", prefix))
    return(prefix)
  }
  prefix <- file.path(ctx$dirs$hisat2_index, "genome")
  if (hisat2_index_files_exist(prefix) && !isTRUE(cfg$rebuild_hisat2_index)) {
    log_count(ctx, sprintf("Using existing HISAT2 index: %s", prefix))
    return(prefix)
  }
  run_count_command(ctx, "hisat2_build_count_index", "hisat2-build", c("-p", as.character(cfg$threads), cfg$genome_file, prefix), check = TRUE)
  prefix
}

strandness_args_for_hisat2 <- function(strandness) {
  if (strandness %in% c("FR", "RF")) c("--rna-strandness", strandness) else character(0)
}

#' Align one Smart-seq2 sample with HISAT2 for count analysis
#'
#' @description
#' Runs HISAT2 for a paired-end or single-end Smart-seq2 sample and writes a SAM
#' file for downstream sorting and featureCounts quantification.
#'
#' @param sample_row One row from the Smart-seq2 manifest.
#' @param cfg A validated count configuration list.
#' @param ctx Count workflow context containing output paths and logs.
#' @param index_prefix HISAT2 index prefix.
#' @return The SAM file path.
#' @export
run_smartseq2_hisat2_alignment <- function(sample_row, cfg, ctx, index_prefix) {
  sid <- sanitize_sample_id(sample_row$sample_id)
  sam_path <- file.path(ctx$dirs$sam, paste0(sid, ".sam"))
  args <- c("--dta", "-x", index_prefix, "-p", as.character(cfg$threads), strandness_args_for_hisat2(cfg$strandness))
  if (sample_row$read_layout_detected == "paired") {
    args <- c(args, "-1", sample_row$R1_path, "-2", sample_row$R2_path)
  } else {
    args <- c(args, "-U", sample_row$single_fastq_path)
  }
  args <- c(args, "-S", sam_path)
  run_count_command(ctx, paste0("hisat2_smartseq2_", sid), "hisat2", args, check = TRUE)
  sam_path
}

#' Convert a Smart-seq2 SAM file to a sorted BAM file
#'
#' @description
#' Converts SAM to BAM with a minimum mapping-quality filter, coordinate-sorts
#' the BAM, and creates a BAM index. Intermediate SAM files are retained.
#'
#' @param sample_row One row from the Smart-seq2 manifest.
#' @param cfg A validated count configuration list.
#' @param ctx Count workflow context containing output paths and logs.
#' @param sam_path Input SAM path.
#' @return The sorted BAM file path.
#' @export
smartseq2_sam_to_sorted_bam <- function(sample_row, cfg, ctx, sam_path) {
  sid <- sanitize_sample_id(sample_row$sample_id)
  tmp_bam <- file.path(ctx$dirs$sorted_bam, paste0(sid, ".unsorted.bam"))
  sorted_bam <- file.path(ctx$dirs$sorted_bam, paste0(sid, ".sorted.bam"))
  run_count_command(ctx, paste0("samtools_view_", sid), "samtools", c("view", "-@", as.character(cfg$threads), "-b", "-q", as.character(cfg$min_mapping_quality), "-o", tmp_bam, sam_path), check = TRUE)
  run_count_command(ctx, paste0("samtools_sort_", sid), "samtools", c("sort", "-@", as.character(cfg$threads), "-o", sorted_bam, tmp_bam), check = TRUE)
  run_count_command(ctx, paste0("samtools_index_", sid), "samtools", c("index", sorted_bam), check = TRUE)
  if (!isTRUE(cfg$dry_run) && file.exists(tmp_bam)) unlink(tmp_bam)
  sorted_bam
}

featurecounts_strand_arg <- function(strandness) {
  if (strandness == "FR") return("1")
  if (strandness == "RF") return("2")
  "0"
}

split_extra_args <- function(x) {
  if (is.null(x) || !nzchar(trimws(x))) return(character(0))
  strsplit(trimws(x), "\\s+")[[1]]
}

#' Run featureCounts for Smart-seq2 gene-level quantification
#'
#' @description
#' Counts gene-level reads against the augmented mRNA + lncRNA GTF. Paired-end
#' runs use \code{-p --countReadPairs}; single-end runs omit paired options.
#'
#' @param bam_files Sorted BAM files for all Smart-seq2 samples.
#' @param cfg A validated count configuration list.
#' @param ctx Count workflow context containing output paths and logs.
#' @param combined_gtf Augmented GTF path.
#' @param manifest Smart-seq2 FASTQ/BAM manifest.
#' @return The featureCounts output file path.
#' @export
run_featurecounts_for_smartseq2 <- function(bam_files, cfg, ctx, combined_gtf, manifest) {
  layout <- unique(manifest$read_layout_detected)
  if (length(layout) != 1) stop("featureCounts requires one read layout per run.")
  args <- c(
    "-T", as.character(cfg$threads),
    "-t", cfg$featurecounts$feature_type,
    "-g", cfg$featurecounts$attribute_type,
    "-s", featurecounts_strand_arg(cfg$strandness),
    "-a", combined_gtf,
    "-o", ctx$files$featurecounts_output
  )
  if (layout == "paired") args <- c(args, "-p", "--countReadPairs")
  if (isTRUE(cfg$featurecounts$allow_multi_overlap)) args <- c(args, "-O")
  if (isTRUE(cfg$featurecounts$count_multi_mapping_reads)) args <- c(args, "-M")
  args <- c(args, split_extra_args(cfg$featurecounts$extra_args), bam_files)
  run_count_command(ctx, "featureCounts_smartseq2", "featureCounts", args, check = TRUE)
  ctx$files$featurecounts_output
}

#' Parse a featureCounts output table into a count matrix
#'
#' @description
#' Converts the featureCounts output table into a gene by sample count matrix,
#' replacing BAM-path column names with Smart-seq2 sample IDs.
#'
#' @param featurecounts_file featureCounts output file.
#' @param cfg A validated count configuration list.
#' @param ctx Count workflow context containing output paths and logs.
#' @param manifest Smart-seq2 manifest with sample IDs and BAM paths.
#' @return A list containing matrix paths and dimensions.
#' @export
parse_featurecounts_matrix <- function(featurecounts_file, cfg, ctx, manifest) {
  if (isTRUE(cfg$dry_run)) {
    log_count(ctx, "[dry_run] Skipping featureCounts matrix parsing.")
    return(list(tsv = ctx$files$count_matrix_tsv, rds = ctx$files$count_matrix_rds, n_genes = NA_integer_, n_samples = nrow(manifest)))
  }
  if (!file.exists(featurecounts_file)) stop("featureCounts output not found: ", featurecounts_file)
  fc <- read.delim(featurecounts_file, comment.char = "#", check.names = FALSE, stringsAsFactors = FALSE)
  if (!("Geneid" %in% colnames(fc))) stop("featureCounts output does not contain Geneid column.")
  annotation_cols <- intersect(c("Geneid", "Chr", "Start", "End", "Strand", "Length"), colnames(fc))
  count_cols <- setdiff(colnames(fc), annotation_cols)
  if (length(count_cols) == 0) stop("No count columns found in featureCounts output.")
  if (length(count_cols) != nrow(manifest)) {
    stop("featureCounts count columns (", length(count_cols), ") do not match manifest samples (", nrow(manifest), ").")
  }
  counts <- fc[, count_cols, drop = FALSE]
  colnames(counts) <- manifest$sample_id[seq_len(ncol(counts))]
  out <- data.frame(gene_id = fc$Geneid, counts, check.names = FALSE)
  write.table(out, ctx$files$count_matrix_tsv, sep = "\t", row.names = FALSE, quote = FALSE)
  mat <- as.matrix(counts)
  rownames(mat) <- fc$Geneid
  saveRDS(mat, ctx$files$count_matrix_rds)
  list(tsv = ctx$files$count_matrix_tsv, rds = ctx$files$count_matrix_rds, n_genes = nrow(mat), n_samples = ncol(mat))
}

#' Run Smart-seq2 count workflow
#'
#' @description
#' Executes the Smart-seq2 count branch: FASTQ discovery, HISAT2 alignment,
#' SAM/BAM conversion, featureCounts quantification, and count matrix parsing.
#'
#' @param cfg A validated count configuration list.
#' @param ctx Count workflow context containing output paths and logs.
#' @param combined_gtf Augmented mRNA + lncRNA GTF path.
#' @return A list with manifest, BAM files, featureCounts output, and matrix paths.
#' @export
run_smartseq2_count <- function(cfg, ctx, combined_gtf) {
  manifest <- discover_smartseq2_fastq_inputs(cfg, ctx)
  validate_smartseq2_manifest(manifest, cfg, ctx)
  add_count_status(ctx, "smartseq2_manifest", "", "OK", ctx$files$smartseq2_manifest)

  index_prefix <- build_hisat2_index_for_count_if_needed(cfg, ctx)
  add_count_status(ctx, "hisat2_index", "", if (cfg$dry_run) "DRY_RUN" else "OK", index_prefix)

  bam_files <- character(0)
  for (i in seq_len(nrow(manifest))) {
    sid <- manifest$sample_id[i]
    add_count_status(ctx, "smartseq2_alignment", sid, "RUNNING", "")
    sam_path <- run_smartseq2_hisat2_alignment(manifest[i, , drop = FALSE], cfg, ctx, index_prefix)
    bam_path <- smartseq2_sam_to_sorted_bam(manifest[i, , drop = FALSE], cfg, ctx, sam_path)
    bam_files <- c(bam_files, bam_path)
    add_count_status(ctx, "smartseq2_alignment", sid, if (cfg$dry_run) "DRY_RUN" else "OK", bam_path)
  }
  manifest$sorted_bam_path <- bam_files
  write.table(manifest, ctx$files$smartseq2_manifest, sep = "\t", row.names = FALSE, quote = FALSE)

  add_count_status(ctx, "featureCounts", "", "RUNNING", "")
  fc_file <- run_featurecounts_for_smartseq2(bam_files, cfg, ctx, combined_gtf, manifest)
  matrix_info <- parse_featurecounts_matrix(fc_file, cfg, ctx, manifest)
  add_count_status(ctx, "featureCounts", "", if (cfg$dry_run) "DRY_RUN" else "OK", fc_file)
  add_count_status(ctx, "smartseq2_count_matrix", "", if (cfg$dry_run) "DRY_RUN" else "OK", matrix_info$tsv)

  list(manifest = manifest, bam_files = bam_files, featurecounts = fc_file, matrix = matrix_info)
}

write_count_report <- function(cfg, ctx, manifest, combined_gtf, transcriptome_path = NA_character_) {
  lines <- c(
    "# scLncR count Run Report",
    "",
    sprintf("- Time: %s", now_string()),
    sprintf("- sequencing_platform: %s", cfg$sequencing_platform),
    sprintf("- count_engine: %s", cfg$count_engine),
    sprintf("- dry_run: %s", cfg$dry_run),
    "",
    "## Input Summary",
    "",
    sprintf("- samples_dirs: `%s`", cfg$samples_dirs),
    sprintf("- samples_info: `%s`", cfg$samples_info),
    sprintf("- genome_file: `%s`", cfg$genome_file),
    sprintf("- known_gtf: `%s`", cfg$known_gtf),
    sprintf("- lnc_gtf: `%s`", cfg$lnc_gtf),
    sprintf("- combined_gtf_used: `%s`", combined_gtf),
    ""
  )
  if (!is.na(transcriptome_path)) {
    lines <- c(lines, sprintf("- cellranger_transcriptome_path: `%s`", transcriptome_path), "")
  }

  lines <- c(lines, "## FASTQ Manifest Preview", "", "```text")
  lines <- c(lines, capture.output(print(utils::head(manifest, 10), row.names = FALSE)))
  if (nrow(manifest) > 10) lines <- c(lines, sprintf("... (%d more rows)", nrow(manifest) - 10))
  lines <- c(lines, "```", "")

  lines <- c(lines, "## Technology Notes", "")
  lines <- c(lines, "- 10x + cellranger is the primary supported path for lncRNA-aware quantification.")
  lines <- c(lines, "- Smart-seq2 + featurecounts uses HISAT2, samtools, and featureCounts to produce a gene x sample count matrix.")
  lines <- c(lines, "- Drop-seq interface is experimental and intended for future STARsolo/Drop-seq-tools integration.")
  lines <- c(lines, ""
  )

  if (cfg$sequencing_platform == "smartseq2") {
    lines <- c(
      lines,
      "## Smart-seq2 Outputs",
      "",
      sprintf("- FASTQ manifest: `%s`", ctx$files$smartseq2_manifest),
      sprintf("- validation_report: `%s`", ctx$files$validation_report),
      sprintf("- sorted BAM directory: `%s`", ctx$dirs$sorted_bam),
      sprintf("- featureCounts output: `%s`", ctx$files$featurecounts_output),
      sprintf("- count matrix TSV: `%s`", ctx$files$count_matrix_tsv),
      sprintf("- count matrix RDS: `%s`", ctx$files$count_matrix_rds),
      "",
      "Smart-seq2 does not use Cell Ranger. The main count matrix is generated by featureCounts from sorted BAM files and the augmented GTF.",
      ""
    )
  }

  lines <- c(lines, "## Step Status", "", "```text")
  if (nrow(ctx$status_rows) > 0) {
    lines <- c(lines, capture.output(print(ctx$status_rows, row.names = FALSE)))
  } else {
    lines <- c(lines, "No status rows.")
  }
  lines <- c(lines, "```", "")

  lines <- c(lines, "## Warnings", "")
  if (length(ctx$warnings) == 0) {
    lines <- c(lines, "- None.", "")
  } else {
    lines <- c(lines, paste0("- ", ctx$warnings), "")
  }

  lines <- c(lines, "## Command Preview", "", "```text")
  if (length(ctx$command_preview) > 0) {
    lines <- c(lines, ctx$command_preview)
  } else {
    lines <- c(lines, "No commands recorded.")
  }
  lines <- c(lines, "```", "")

  writeLines(lines, ctx$files$report)
}

print_count_help <- function() {
  cat("scLncR count - raw FASTQ-first technology-aware quantification\n")
  cat("Usage: scLncR count -c <config.yaml>\n\n")
  cat("Primary support:\n")
  cat("  - sequencing_platform=10x + count_engine=cellranger\n")
  cat("  - sequencing_platform=smartseq2 + count_engine=featurecounts\n")
  cat("Experimental interfaces:\n")
  cat("  - dropseq + starsolo/dropseq-tools\n")
}

run_count <- function(user_args, script_dir) {
  option_list <- list(
    optparse::make_option(c("-c", "--config"), type = "character", metavar = "FILE", help = "Path to count YAML config")
  )
  parser <- optparse::OptionParser(
    usage = "scLncR count [options]",
    option_list = option_list,
    description = "count: raw FASTQ-first, lncRNA-aware quantification interface"
  )
  opt <- optparse::parse_args(parser, args = user_args, print_help_and_exit = FALSE)
  if (length(user_args) == 0 || any(user_args %in% c("-h", "--help")) || is.null(opt$config)) {
    print_count_help()
    cat("\n")
    optparse::print_help(parser)
    return(invisible(NULL))
  }
  if (!file.exists(opt$config)) stop("Config file not found: ", opt$config)

  cfg <- load_count_config(opt$config)
  ctx <- init_count_context(cfg)
  copy_if_exists <- function(src, dst) if (!is.null(src) && nzchar(src) && file.exists(src)) file.copy(src, dst, overwrite = TRUE)
  copy_if_exists(opt$config, file.path(ctx$dirs$manifest, "config_Count.snapshot.yaml"))
  copy_if_exists(cfg$samples_info, file.path(ctx$dirs$manifest, "samples_info.copy.tsv"))

  log_count(ctx, "Starting count workflow.")
  combined_gtf <- build_combined_gtf_for_count(cfg, ctx)
  add_count_status(ctx, "build_combined_gtf", "", "OK", combined_gtf)

  manifest <- NULL
  transcriptome_path <- NA_character_
  if (cfg$sequencing_platform == "10x" && cfg$count_engine == "cellranger") {
    manifest <- discover_count_fastq_manifest(cfg, ctx)
    add_count_status(ctx, "cellranger_mkref", "", "RUNNING", "")
    transcriptome_path <- build_cellranger_reference(cfg, ctx, combined_gtf)
    add_count_status(ctx, "cellranger_mkref", "", if (cfg$dry_run) "DRY_RUN" else "OK", transcriptome_path)

    for (i in seq_len(nrow(manifest))) {
      sid <- manifest$sample_id[i]
      add_count_status(ctx, "cellranger_count", sid, "RUNNING", "")
      run_cellranger_count_for_sample(sid, cfg, ctx, transcriptome_path)
      add_count_status(ctx, "cellranger_count", sid, if (cfg$dry_run) "DRY_RUN" else "OK", "")
    }
  } else if (cfg$sequencing_platform == "smartseq2" && cfg$count_engine == "featurecounts") {
    smartseq2_res <- run_smartseq2_count(cfg, ctx, combined_gtf)
    manifest <- smartseq2_res$manifest
  } else {
    msg <- sprintf("Interface %s + %s is experimental and not fully implemented in this version.", cfg$sequencing_platform, cfg$count_engine)
    add_count_warning(ctx, msg)
    add_count_status(ctx, "experimental_interface", "", if (cfg$dry_run) "DRY_RUN" else "SKIPPED", msg)
    if (!cfg$dry_run) {
      stop(msg, " Use dry_run=true for interface validation/report only.")
    }
    manifest <- discover_count_fastq_manifest(cfg, ctx)
  }

  write_count_report(cfg, ctx, manifest, combined_gtf, transcriptome_path)
  log_count(ctx, "count workflow completed.")
  invisible(list(report = ctx$files$report, combined_gtf = combined_gtf))
}
