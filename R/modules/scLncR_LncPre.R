#' Null-coalescing helper
#'
#' @param x Candidate value.
#' @param y Fallback value.
#' @return \code{x} when available; otherwise \code{y}.
`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0) return(y)
  x
}

now_string <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S")

append_line <- function(path, txt) {
  cat(txt, file = path, append = TRUE, sep = "\n")
}

sanitize_sample_id <- function(x) {
  gsub("[^A-Za-z0-9._-]", "-", x)
}

#' Initialize prelnc runtime context and output directories
#'
#' @description
#' Creates a standardized output structure for raw FASTQ-first prelnc analysis
#' and initializes shared runtime objects (logs, status table, warnings).
#'
#' @param cfg A validated prelnc configuration list.
#' @return Environment containing runtime paths and mutable state.
#' @export
init_prelnc_context <- function(cfg) {
  dirs <- list(
    base = cfg$output_path,
    manifest = file.path(cfg$output_path, "manifest"),
    logs = file.path(cfg$output_path, "logs"),
    index_hisat2 = file.path(cfg$output_path, "index", "hisat2"),
    alignment_sam = file.path(cfg$output_path, "alignment", "sam"),
    alignment_bam = file.path(cfg$output_path, "alignment", "sorted_bam"),
    assembly_root = file.path(cfg$output_path, "assembly", "stringtie"),
    assembly_per_sample = file.path(cfg$output_path, "assembly", "stringtie", "per_sample_gtf"),
    gffcompare = file.path(cfg$output_path, "gffcompare"),
    cpc2 = file.path(cfg$output_path, "cpc2"),
    reference = file.path(cfg$output_path, "reference")
  )
  for (d in unlist(dirs, use.names = FALSE)) {
    dir.create(d, recursive = TRUE, showWarnings = FALSE)
  }

  ctx <- new.env(parent = emptyenv())
  ctx$cfg <- cfg
  ctx$dirs <- dirs
  ctx$files <- list(
    commands_log = file.path(dirs$logs, "commands.log"),
    run_log = file.path(dirs$logs, "prelnc_run.log"),
    per_sample_status = file.path(dirs$logs, "per_sample_status.tsv"),
    raw_fastq_manifest = file.path(dirs$manifest, "prelnc_raw_fastq_manifest.tsv"),
    validation_report = file.path(dirs$manifest, "prelnc_input_validation_report.md"),
    report = file.path(dirs$base, "prelnc_run_report.md"),
    final_gtf = file.path(dirs$base, "final_lnc.gtf"),
    final_fa = file.path(dirs$base, "final.lncRNA.fa"),
    combined_gtf = file.path(dirs$reference, "combined_mRNA_lncRNA.gtf")
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
  write.table(ctx$status_rows, file = ctx$files$per_sample_status, sep = "\t", row.names = FALSE, quote = FALSE)
  ctx
}

log_prelnc <- function(ctx, msg, level = "INFO") {
  line <- sprintf("[%s] [%s] %s", now_string(), level, msg)
  append_line(ctx$files$run_log, line)
  message(line)
}

add_warning <- function(ctx, msg) {
  ctx$warnings <- unique(c(ctx$warnings, msg))
  log_prelnc(ctx, msg, level = "WARN")
}

add_status <- function(ctx, step, sample_id = "", status = "OK", message = "") {
  ctx$status_rows <- rbind(
    ctx$status_rows,
    data.frame(step = step, sample_id = sample_id, status = status, message = message, stringsAsFactors = FALSE)
  )
  write.table(ctx$status_rows, file = ctx$files$per_sample_status, sep = "\t", row.names = FALSE, quote = FALSE)
}

record_command <- function(ctx, step, command, args = character(0)) {
  cmd_str <- if (length(args) == 0) command else paste(c(command, shQuote(args)), collapse = " ")
  append_line(ctx$files$commands_log, sprintf("[%s] [%s] %s", now_string(), step, cmd_str))
  ctx$command_preview <- c(ctx$command_preview, sprintf("[%s] %s", step, cmd_str))
  cmd_str
}

run_command <- function(ctx, step, command, args = character(0), check = TRUE) {
  cmd_str <- record_command(ctx, step, command, args)
  if (isTRUE(ctx$cfg$dry_run)) {
    log_prelnc(ctx, sprintf("[dry_run] %s", cmd_str))
    return(list(ok = TRUE, status = NA_integer_, output = character(0)))
  }

  out_tmp <- tempfile("prelnc_cmd_")
  res <- tryCatch(
    {
      status <- suppressWarnings(system2(command, args = args, stdout = out_tmp, stderr = out_tmp))
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
    log_prelnc(ctx, sprintf("Command failed [%s] exit=%s", step, res$status), level = "ERROR")
    if (check) stop(sprintf("Command failed at step '%s': %s", step, cmd_str))
  }
  list(ok = ok, status = res$status, output = out)
}

copy_if_exists <- function(src, dst) {
  if (!is.null(src) && nzchar(src) && file.exists(src)) file.copy(src, dst, overwrite = TRUE)
}

parse_fastq_identity <- function(filename, cfg) {
  pattern <- cfg$fastq_sample_regex
  m <- regexec(pattern, filename, perl = TRUE, ignore.case = TRUE)
  reg <- regmatches(filename, m)[[1]]
  if (length(reg) >= 3) {
    return(list(sample_id = reg[2], token = toupper(reg[3])))
  }

  token <- "SINGLE"
  if (grepl(cfg$r2_pattern, filename, fixed = TRUE)) token <- "R2"
  if (grepl(cfg$r1_pattern, filename, fixed = TRUE)) token <- "R1"
  if (grepl(cfg$i1_pattern, filename, fixed = TRUE)) token <- "I1"
  sample_id <- sub("\\.(fastq|fq)(\\.gz)?$", "", filename, ignore.case = TRUE)
  sample_id <- sub("(_S\\d+_L\\d+_[IR][12]_\\d+)$", "", sample_id, perl = TRUE, ignore.case = TRUE)
  sample_id <- sub("(_R1_|_R2_|_I1_|_I2_).*$", "", sample_id, perl = TRUE, ignore.case = TRUE)
  list(sample_id = sample_id, token = token)
}

#' Discover raw FASTQ inputs for scLncR prelnc analysis
#'
#' @description
#' Scans a raw FASTQ directory and builds a technology-aware sample manifest
#' for 10x, Smart-seq2, Drop-seq, or generic sequencing data.
#'
#' @param cfg A validated prelnc configuration list.
#' @return A data.frame containing sample IDs, FASTQ paths, read roles, and notes.
#' @examples
#' manifest <- discover_raw_fastq_inputs(cfg)
#' @export
discover_raw_fastq_inputs <- function(cfg) {
  if (!dir.exists(cfg$samples_dirs)) stop("samples_dirs not found: ", cfg$samples_dirs)
  all_files <- list.files(cfg$samples_dirs, recursive = TRUE, full.names = TRUE)
  fastq_files <- all_files[grepl("\\.(fastq|fq)(\\.gz)?$", all_files, ignore.case = TRUE)]
  if (length(fastq_files) == 0) stop("No FASTQ files detected in samples_dirs.")

  rows <- lapply(fastq_files, function(f) {
    id <- parse_fastq_identity(basename(f), cfg)
    data.frame(
      sample_id = id$sample_id,
      token = id$token,
      fastq_path = normalizePath(f, winslash = "/", mustWork = FALSE),
      stringsAsFactors = FALSE
    )
  })
  parsed <- do.call(rbind, rows)

  sample_ids <- unique(parsed$sample_id)
  if (length(cfg$include_samples) > 0) sample_ids <- intersect(sample_ids, cfg$include_samples)
  if (cfg$sample_limit > 0) sample_ids <- sample_ids[seq_len(min(cfg$sample_limit, length(sample_ids)))]
  parsed <- parsed[parsed$sample_id %in% sample_ids, , drop = FALSE]
  if (nrow(parsed) == 0) stop("No samples left after include_samples/sample_limit filtering.")

  out <- data.frame(
    sample_id = sort(unique(parsed$sample_id)),
    sequencing_platform = cfg$sequencing_platform,
    I1_path = NA_character_,
    R1_path = NA_character_,
    R2_path = NA_character_,
    single_fastq_path = NA_character_,
    read_layout_detected = NA_character_,
    discovery_reads_used = NA_character_,
    barcode_read = NA_character_,
    umi_read = NA_character_,
    index_read = NA_character_,
    notes = NA_character_,
    stringsAsFactors = FALSE
  )

  pick_one <- function(df, token) {
    vals <- df$fastq_path[df$token == token]
    if (length(vals) == 0) NA_character_ else vals[1]
  }

  for (i in seq_len(nrow(out))) {
    sid <- out$sample_id[i]
    sub <- parsed[parsed$sample_id == sid, , drop = FALSE]
    out$I1_path[i] <- pick_one(sub, "I1")
    out$R1_path[i] <- pick_one(sub, "R1")
    out$R2_path[i] <- pick_one(sub, "R2")
    out$single_fastq_path[i] <- pick_one(sub, "SINGLE")

    if (cfg$sequencing_platform == "10x") {
      out$read_layout_detected[i] <- "single_cDNA_R2"
      out$discovery_reads_used[i] <- out$R2_path[i]
      out$barcode_read[i] <- out$R1_path[i]
      out$umi_read[i] <- out$R1_path[i]
      out$index_read[i] <- out$I1_path[i]
      out$notes[i] <- "10x mode: R2 used for transcript evidence; R1/I1 retained for traceability."
    } else if (cfg$sequencing_platform == "dropseq") {
      layout <- if (!is.na(out$R1_path[i]) && !is.na(out$R2_path[i])) "paired" else "single"
      if (cfg$read_layout != "auto") layout <- cfg$read_layout
      out$read_layout_detected[i] <- layout
      out$discovery_reads_used[i] <- if (!is.na(out$R2_path[i])) out$R2_path[i] else if (!is.na(out$single_fastq_path[i])) out$single_fastq_path[i] else out$R1_path[i]
      out$barcode_read[i] <- out$R1_path[i]
      out$umi_read[i] <- out$R1_path[i]
      out$index_read[i] <- out$I1_path[i]
      out$notes[i] <- "Drop-seq interface is experimental; full UMI-aware quantification is planned in count module."
    } else {
      layout <- if (!is.na(out$R1_path[i]) && !is.na(out$R2_path[i])) "paired" else "single"
      if (cfg$read_layout != "auto") layout <- cfg$read_layout
      out$read_layout_detected[i] <- layout
      if (layout == "paired") {
        out$discovery_reads_used[i] <- paste(out$R1_path[i], out$R2_path[i], sep = ";")
      } else {
        out$discovery_reads_used[i] <- if (!is.na(out$single_fastq_path[i])) out$single_fastq_path[i] else if (!is.na(out$R2_path[i])) out$R2_path[i] else out$R1_path[i]
      }
      out$notes[i] <- if (cfg$sequencing_platform == "smartseq2") {
        "Smart-seq2 mode: full-length-like transcript reconstruction depends on read length, coverage, and alignment quality."
      } else {
        "Generic mode: discovery alignment follows detected/provided read layout."
      }
    }
  }
  out
}

#' Validate raw FASTQ sample manifest
#'
#' @param manifest A FASTQ manifest from \code{discover_raw_fastq_inputs()}.
#' @param cfg A validated prelnc configuration list.
#' @return List with \code{errors}, \code{warnings}, and original \code{manifest}.
#' @export
validate_raw_fastq_manifest <- function(manifest, cfg) {
  errs <- character(0)
  warns <- character(0)
  if (anyDuplicated(manifest$sample_id) > 0) errs <- c(errs, "Duplicated sample_id detected in FASTQ manifest.")

  for (i in seq_len(nrow(manifest))) {
    r <- manifest[i, , drop = FALSE]
    sid <- r$sample_id
    check_file <- function(x) is.na(x) || file.exists(x)
    for (nm in c("I1_path", "R1_path", "R2_path", "single_fastq_path")) {
      if (!check_file(r[[nm]])) errs <- c(errs, sprintf("Missing FASTQ for sample '%s': %s", sid, r[[nm]]))
    }

    if (cfg$sequencing_platform == "10x") {
      if (is.na(r$R2_path) || !nzchar(r$R2_path)) errs <- c(errs, sprintf("10x sample '%s' requires R2 FASTQ.", sid))
      if (is.na(r$R1_path)) warns <- c(warns, sprintf("10x sample '%s' missing R1 (barcode/UMI metadata).", sid))
      if (is.na(r$I1_path)) warns <- c(warns, sprintf("10x sample '%s' missing I1 (index metadata).", sid))
    } else if (cfg$sequencing_platform == "smartseq2" || cfg$sequencing_platform == "generic") {
      layout <- tolower(r$read_layout_detected)
      if (layout == "paired" && (is.na(r$R1_path) || is.na(r$R2_path))) {
        errs <- c(errs, sprintf("%s sample '%s' requires R1 and R2 for paired layout.", cfg$sequencing_platform, sid))
      }
      if (layout == "single" && is.na(r$single_fastq_path) && is.na(r$R1_path) && is.na(r$R2_path)) {
        errs <- c(errs, sprintf("%s sample '%s' has no usable single-end FASTQ.", cfg$sequencing_platform, sid))
      }
    } else if (cfg$sequencing_platform == "dropseq") {
      if (is.na(r$R2_path) && is.na(r$single_fastq_path) && is.na(r$R1_path)) {
        errs <- c(errs, sprintf("Drop-seq sample '%s' has no usable discovery read.", sid))
      }
      warns <- c(warns, "Drop-seq full UMI-aware quantification is experimental/planned; prelnc discovery outputs should be interpreted cautiously.")
    }
  }
  list(errors = unique(errs), warnings = unique(warns), manifest = manifest)
}

write_validation_report <- function(ctx, val_res) {
  lines <- c(
    "# prelnc Input Validation Report",
    "",
    sprintf("- Time: %s", now_string()),
    sprintf("- sequencing_platform: %s", ctx$cfg$sequencing_platform),
    sprintf("- dry_run: %s", ctx$cfg$dry_run),
    ""
  )
  lines <- c(lines, sprintf("- errors: %d", length(val_res$errors)))
  lines <- c(lines, sprintf("- warnings: %d", length(val_res$warnings)), "")
  if (length(val_res$errors) > 0) {
    lines <- c(lines, "## Errors", paste0("- ", val_res$errors), "")
  }
  if (length(val_res$warnings) > 0) {
    lines <- c(lines, "## Warnings", paste0("- ", val_res$warnings), "")
  }
  writeLines(lines, ctx$files$validation_report)
}

#' Build HISAT2 index if required for discovery alignment
#'
#' @param cfg A validated prelnc configuration list.
#' @param ctx Runtime context.
#' @return Index prefix path.
#' @export
build_hisat2_index_if_needed <- function(cfg, ctx) {
  prefix <- cfg$hisat2_index_prefix
  if (!nzchar(prefix)) prefix <- file.path(ctx$dirs$index_hisat2, "prelnc_genome_index")
  idx_files <- paste0(prefix, c(".1.ht2", ".2.ht2", ".3.ht2", ".4.ht2", ".5.ht2", ".6.ht2", ".7.ht2", ".8.ht2"))
  need_build <- isTRUE(cfg$rebuild_hisat2_index) || !all(file.exists(idx_files))
  if (!need_build) {
    log_prelnc(ctx, sprintf("Reusing HISAT2 index: %s", prefix))
    return(prefix)
  }
  run_command(ctx, "hisat2_build_index", "hisat2-build", c(cfg$genome_file, prefix), check = TRUE)
  prefix
}

hisat2_strand_args <- function(strandness) {
  if (strandness == "unstranded") character(0) else c("--rna-strandness", strandness)
}

stringtie_strand_args <- function(strandness) {
  if (strandness == "RF") return("--rf")
  if (strandness == "FR") return("--fr")
  character(0)
}

#' Run technology-aware discovery alignment for one sample
#'
#' @description
#' Generates HISAT2 alignment commands according to sequencing technology.
#' In 10x mode, only R2 cDNA reads are used for transcript evidence.
#'
#' @param sample_row One-row data.frame from FASTQ manifest.
#' @param cfg A validated prelnc configuration list.
#' @param ctx Runtime context.
#' @param index_prefix HISAT2 index prefix.
#' @return List with sample ID and SAM path.
#' @export
run_discovery_alignment_for_sample <- function(sample_row, cfg, ctx, index_prefix) {
  sid <- sanitize_sample_id(sample_row$sample_id)
  sam_path <- file.path(ctx$dirs$alignment_sam, paste0(sid, ".sam"))
  args <- c("--dta", "-x", index_prefix, "-p", as.character(cfg$threads), hisat2_strand_args(cfg$strandness))

  if (cfg$sequencing_platform == "10x") {
    args <- c(args, "-U", sample_row$R2_path, "-S", sam_path)
    log_prelnc(ctx, sprintf("Sample %s: 10x mode uses R2 cDNA reads for candidate transcript evidence.", sample_row$sample_id))
  } else if (cfg$sequencing_platform == "smartseq2") {
    if (tolower(sample_row$read_layout_detected) == "paired") {
      args <- c(args, "-1", sample_row$R1_path, "-2", sample_row$R2_path, "-S", sam_path)
    } else {
      single <- sample_row$single_fastq_path
      if (is.na(single)) single <- if (!is.na(sample_row$R2_path)) sample_row$R2_path else sample_row$R1_path
      args <- c(args, "-U", single, "-S", sam_path)
    }
  } else if (cfg$sequencing_platform == "dropseq") {
    add_warning(ctx, "Drop-seq discovery interface is experimental. Using available read as provisional transcript evidence.")
    use_read <- if (!is.na(sample_row$R2_path)) sample_row$R2_path else if (!is.na(sample_row$single_fastq_path)) sample_row$single_fastq_path else sample_row$R1_path
    args <- c(args, "-U", use_read, "-S", sam_path)
  } else {
    if (tolower(sample_row$read_layout_detected) == "paired") {
      args <- c(args, "-1", sample_row$R1_path, "-2", sample_row$R2_path, "-S", sam_path)
    } else {
      single <- sample_row$single_fastq_path
      if (is.na(single)) single <- if (!is.na(sample_row$R2_path)) sample_row$R2_path else sample_row$R1_path
      args <- c(args, "-U", single, "-S", sam_path)
    }
  }
  run_command(ctx, paste0("hisat2_", sid), "hisat2", args, check = TRUE)
  list(sample_id = sample_row$sample_id, sam_path = sam_path)
}

#' Convert SAM to coordinate-sorted BAM
#'
#' @param sample_row One-row list/data.frame with \code{sample_id}, \code{sam_path}.
#' @param cfg A validated prelnc configuration list.
#' @param ctx Runtime context.
#' @return List with BAM path for downstream assembly.
#' @export
sam_to_sorted_bam <- function(sample_row, cfg, ctx) {
  sid <- sanitize_sample_id(sample_row$sample_id)
  sam_path <- sample_row$sam_path
  unsorted_bam <- file.path(ctx$dirs$alignment_bam, paste0(sid, ".unsorted.bam"))
  sorted_bam <- file.path(ctx$dirs$alignment_bam, paste0(sid, ".sorted.bam"))

  run_command(
    ctx, paste0("samtools_view_", sid), "samtools",
    c("view", "-@", as.character(cfg$threads), "-bS", "-q", as.character(cfg$min_mapping_quality), "-o", unsorted_bam, sam_path),
    check = TRUE
  )
  run_command(
    ctx, paste0("samtools_sort_", sid), "samtools",
    c("sort", "-@", as.character(cfg$threads), "-o", sorted_bam, unsorted_bam),
    check = TRUE
  )
  run_command(ctx, paste0("samtools_index_", sid), "samtools", c("index", sorted_bam), check = FALSE)

  if (!isTRUE(cfg$keep_intermediate) && !isTRUE(cfg$dry_run)) {
    if (file.exists(sam_path)) file.remove(sam_path)
    if (file.exists(unsorted_bam)) file.remove(unsorted_bam)
  }

  list(sample_id = sample_row$sample_id, bam_path = sorted_bam)
}

#' Run per-sample StringTie transcript assembly
#'
#' @param sample_row One-row object with \code{sample_id} and \code{bam_path}.
#' @param cfg A validated prelnc configuration list.
#' @param ctx Runtime context.
#' @return List with sample ID and per-sample GTF path.
#' @export
run_stringtie_for_sample <- function(sample_row, cfg, ctx) {
  sid <- sanitize_sample_id(sample_row$sample_id)
  out_gtf <- file.path(ctx$dirs$assembly_per_sample, paste0(sid, ".gtf"))
  args <- c(
    "-p", as.character(cfg$threads),
    stringtie_strand_args(cfg$strandness),
    "-G", cfg$gtf_file,
    "-o", out_gtf,
    "-l", paste0(cfg$lncrna_name, "_", sid),
    sample_row$bam_path
  )
  run_command(ctx, paste0("stringtie_", sid), "stringtie", args, check = TRUE)
  list(sample_id = sample_row$sample_id, gtf_path = out_gtf)
}

#' Merge per-sample StringTie GTF files
#'
#' @param gtf_files Character vector of sample-level GTF files.
#' @param cfg A validated prelnc configuration list.
#' @param ctx Runtime context.
#' @return Path to merged GTF file.
#' @export
merge_stringtie_gtfs <- function(gtf_files, cfg, ctx) {
  list_file <- file.path(ctx$dirs$assembly_root, "gtf_list.txt")
  writeLines(gtf_files, list_file)
  merged_gtf <- file.path(ctx$dirs$assembly_root, "all.merged.gtf")
  args <- c("--merge", "-p", as.character(cfg$threads), stringtie_strand_args(cfg$strandness), "-G", cfg$gtf_file, "-o", merged_gtf, list_file)
  run_command(ctx, "stringtie_merge", "stringtie", args, check = TRUE)
  merged_gtf
}

subset_gtf_by_transcript_ids <- function(gtf_in, transcript_ids, gtf_out) {
  lines <- readLines(gtf_in, warn = FALSE)
  keep <- logical(length(lines))
  tx_re <- '.*transcript_id "([^"]+)".*'
  for (i in seq_along(lines)) {
    ln <- lines[i]
    if (startsWith(ln, "#")) next
    if (!grepl("transcript_id \"", ln, fixed = TRUE)) next
    tx <- sub(tx_re, "\\1", ln)
    keep[i] <- tx %in% transcript_ids
  }
  writeLines(lines[keep], gtf_out)
}

extract_gtf_attr <- function(line, attr_name) {
  pattern <- sprintf('.*%s "([^"]+)".*', attr_name)
  if (!grepl(sprintf('%s "', attr_name), line, fixed = TRUE)) return(NA_character_)
  sub(pattern, "\\1", line)
}

set_gtf_attr <- function(line, attr_name, value) {
  if (grepl(sprintf('%s "', attr_name), line, fixed = TRUE)) {
    sub(sprintf('%s "[^"]+"', attr_name), sprintf('%s "%s"', attr_name, value), line)
  } else {
    sub(";[[:space:]]*$", sprintf('; %s "%s";', attr_name, value), line)
  }
}

rename_gtf_transcripts <- function(gtf_in, id_map, gtf_out) {
  lines <- readLines(gtf_in, warn = FALSE)
  if (length(lines) == 0) {
    writeLines(character(0), gtf_out)
    return(invisible(FALSE))
  }

  tx_lookup <- setNames(seq_len(nrow(id_map)), id_map$old_transcript_id)
  out <- character(0)
  for (ln in lines) {
    if (startsWith(ln, "#") || !grepl("transcript_id \"", ln, fixed = TRUE)) next
    old_tx <- extract_gtf_attr(ln, "transcript_id")
    if (is.na(old_tx) || !nzchar(old_tx)) next
    idx <- unname(tx_lookup[old_tx])
    if (length(idx) == 0 || is.na(idx)) next
    ln <- set_gtf_attr(ln, "gene_id", id_map$new_gene_id[idx])
    ln <- set_gtf_attr(ln, "transcript_id", id_map$new_transcript_id[idx])
    out <- c(out, ln)
  }
  writeLines(out, gtf_out)
  invisible(length(out) > 0)
}

select_cpc2_columns <- function(cpc) {
  col_lower <- tolower(colnames(cpc))
  id_idx <- which(col_lower %in% c("#id", "id", "transcript_id", "sequence_name", "seq_id"))
  if (length(id_idx) == 0) id_idx <- 1L

  label_idx <- which(col_lower %in% c("label", "coding_label", "prediction", "class", "coding_status"))
  if (length(label_idx) == 0) {
    label_idx <- which(grepl("label|prediction|class|status", col_lower))
  }
  if (length(label_idx) == 0) stop("Cannot identify CPC2 label column. Columns: ", paste(colnames(cpc), collapse = ", "))

  prob_idx <- which(col_lower %in% c("coding_probability", "coding_prob", "probability"))
  list(id = id_idx[1], label = label_idx[1], probability = if (length(prob_idx) == 0) NA_integer_ else prob_idx[1])
}

#' Filter lncRNA candidates by gffcompare class and transcript length
#'
#' @description
#' Applies unambiguous filtering logic:
#' class_code in configured classes AND transcript length > minimum length.
#'
#' @param merged_gtf Path to merged StringTie GTF.
#' @param cfg A validated prelnc configuration list.
#' @param ctx Runtime context.
#' @return List with candidate IDs and candidate GTF path.
#' @export
filter_lncRNA_candidates <- function(merged_gtf, cfg, ctx) {
  prefix <- file.path(ctx$dirs$gffcompare, "gffcmp")
  run_command(ctx, "gffcompare", "gffcompare", c("-r", cfg$gtf_file, "-o", prefix, merged_gtf), check = TRUE)

  expected <- c(
    paste0(prefix, ".", basename(merged_gtf), ".tmap"),
    paste0(prefix, ".all.merged.gtf.tmap"),
    paste0(prefix, ".merged.gtf.tmap"),
    paste0(prefix, ".tmap")
  )
  tmap_file <- expected[file.exists(expected)][1]
  if (is.na(tmap_file) || !nzchar(tmap_file)) {
    tmap_candidates <- c(
      list.files(ctx$dirs$gffcompare, pattern = "\\.tmap$", full.names = TRUE),
      list.files(dirname(merged_gtf), pattern = "^gffcmp\\..*\\.tmap$", full.names = TRUE)
    )
    tmap_candidates <- unique(tmap_candidates[file.exists(tmap_candidates)])
    if (length(tmap_candidates) > 0) {
      mt <- file.info(tmap_candidates)$mtime
      tmap_file <- tmap_candidates[which.max(mt)]
      log_prelnc(ctx, sprintf("Using detected gffcompare tmap file: %s", tmap_file))
    } else {
      produced <- unique(c(
        file.path(ctx$dirs$gffcompare, list.files(ctx$dirs$gffcompare, full.names = FALSE)),
        file.path(dirname(merged_gtf), list.files(dirname(merged_gtf), pattern = "^gffcmp\\.", full.names = FALSE))
      ))
      stop(
        "Cannot find gffcompare tmap output. Candidate gffcompare files: ",
        paste(produced, collapse = ", ")
      )
    }
  }

  # Archive detected tmap into gffcompare directory for stable downstream paths.
  archived_tmap <- file.path(ctx$dirs$gffcompare, basename(tmap_file))
  if (normalizePath(tmap_file, winslash = "/", mustWork = FALSE) != normalizePath(archived_tmap, winslash = "/", mustWork = FALSE)) {
    file.copy(tmap_file, archived_tmap, overwrite = TRUE)
    tmap_file <- archived_tmap
  }

  annotated_gtf <- paste0(prefix, ".annotated.gtf")
  if (!file.exists(annotated_gtf)) {
    alt_annot <- file.path(dirname(merged_gtf), "gffcmp.annotated.gtf")
    if (file.exists(alt_annot)) {
      file.copy(alt_annot, annotated_gtf, overwrite = TRUE)
    } else {
      stop("Cannot find gffcompare annotated GTF output.")
    }
  }

  tmap <- read.table(tmap_file, sep = "\t", header = FALSE, stringsAsFactors = FALSE, quote = "", comment.char = "")
  if (ncol(tmap) < 10) stop("Unexpected gffcompare tmap format: ", tmap_file)
  class_code <- tmap[[3]]
  query_id <- tmap[[5]]
  tx_len <- suppressWarnings(as.numeric(tmap[[10]]))

  keep <- class_code %in% cfg$gffcompare_classes & !is.na(tx_len) & tx_len > cfg$min_transcript_length
  candidate_ids <- unique(query_id[keep & nzchar(query_id)])
  if (length(candidate_ids) == 0) stop("No transcripts passed gffcompare class/length filtering.")

  id_file <- file.path(ctx$dirs$gffcompare, "candidate_transcript_ids.txt")
  writeLines(candidate_ids, id_file)
  candidate_gtf <- file.path(ctx$dirs$gffcompare, "candidate_lnc.gtf")
  subset_gtf_by_transcript_ids(annotated_gtf, candidate_ids, candidate_gtf)

  list(candidate_ids = candidate_ids, candidate_id_file = id_file, candidate_gtf = candidate_gtf)
}

#' Filter coding potential using CPC2 and produce final lncRNA annotation
#'
#' @param candidate_gtf Candidate GTF after structural filtering.
#' @param cfg A validated prelnc configuration list.
#' @param ctx Runtime context.
#' @return List with final GTF and FASTA paths.
#' @export
filter_coding_potential_with_CPC2 <- function(candidate_gtf, cfg, ctx) {
  candidate_fa <- file.path(ctx$dirs$cpc2, "candidate_lncRNA.fa")
  run_command(ctx, "gffread_candidates", "gffread", c(candidate_gtf, "-g", cfg$genome_file, "-w", candidate_fa), check = TRUE)

  cpc_prefix <- file.path(ctx$dirs$cpc2, "CPC2.out")
  run_command(ctx, "cpc2_predict", "CPC2.py", c("-i", candidate_fa, "-o", cpc_prefix, "--ORF"), check = TRUE)

  cpc_txt <- paste0(cpc_prefix, ".txt")
  if (!file.exists(cpc_txt)) stop("CPC2 output not found: ", cpc_txt)
  cpc <- read.delim(cpc_txt, sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, comment.char = "")
  if (nrow(cpc) == 0) stop("CPC2 output is empty.")

  cpc_cols <- select_cpc2_columns(cpc)
  id_col <- colnames(cpc)[cpc_cols$id]
  label_col <- colnames(cpc)[cpc_cols$label]
  prob_col <- if (is.na(cpc_cols$probability)) NA_character_ else colnames(cpc)[cpc_cols$probability]

  tx_id <- as.character(cpc[[cpc_cols$id]])
  pred <- tolower(trimws(as.character(cpc[[cpc_cols$label]])))
  noncoding_ids <- unique(tx_id[pred %in% c("noncoding", "non-coding", "non coding") & nzchar(tx_id)])
  label_counts <- sort(table(pred), decreasing = TRUE)

  cpc_summary <- data.frame(
    metric = c(
      "total_transcripts",
      "id_column",
      "label_column",
      "coding_probability_column",
      "noncoding_count"
    ),
    value = c(nrow(cpc), id_col, label_col, prob_col %||% "", length(noncoding_ids)),
    stringsAsFactors = FALSE
  )
  label_summary <- data.frame(
    metric = paste0("label_count:", names(label_counts)),
    value = as.integer(label_counts),
    stringsAsFactors = FALSE
  )
  cpc_summary_file <- file.path(ctx$dirs$cpc2, "CPC2_filter_summary.tsv")
  write.table(rbind(cpc_summary, label_summary), cpc_summary_file, sep = "\t", row.names = FALSE, quote = FALSE)

  if (length(noncoding_ids) == 0) {
    stop(
      "No noncoding transcripts after CPC2 filtering. ",
      "Parsed label column: ", label_col, ". ",
      "Observed labels: ", paste(sprintf("%s=%s", names(label_counts), as.integer(label_counts)), collapse = ", ")
    )
  }

  noncoding_id_file <- file.path(ctx$dirs$cpc2, "noncoding_transcript_ids.txt")
  writeLines(noncoding_ids, noncoding_id_file)

  id_map <- data.frame(
    old_transcript_id = noncoding_ids,
    new_gene_id = sprintf("%s%06d", cfg$lncrna_name, seq_along(noncoding_ids)),
    new_transcript_id = sprintf("%s%06d.1", cfg$lncrna_name, seq_along(noncoding_ids)),
    stringsAsFactors = FALSE
  )
  if (!is.na(cpc_cols$probability)) {
    prob_lookup <- setNames(as.character(cpc[[cpc_cols$probability]]), tx_id)
    id_map$coding_probability <- unname(prob_lookup[id_map$old_transcript_id])
  }
  pred_lookup <- setNames(as.character(cpc[[cpc_cols$label]]), tx_id)
  id_map$cpc2_label <- unname(pred_lookup[id_map$old_transcript_id])

  id_map_file <- file.path(ctx$dirs$cpc2, "final_lnc_id_mapping.tsv")
  write.table(id_map, id_map_file, sep = "\t", row.names = FALSE, quote = FALSE)

  ok <- rename_gtf_transcripts(candidate_gtf, id_map, ctx$files$final_gtf)
  if (!isTRUE(ok)) stop("No GTF records were written after CPC2 noncoding ID mapping.")

  run_command(ctx, "gffread_final_fasta", "gffread", c(ctx$files$final_gtf, "-g", cfg$genome_file, "-w", ctx$files$final_fa), check = TRUE)

  list(
    final_gtf = ctx$files$final_gtf,
    final_fa = ctx$files$final_fa,
    noncoding_id_file = noncoding_id_file,
    id_mapping_file = id_map_file,
    cpc_summary_file = cpc_summary_file
  )
}

#' Build augmented known+lncRNA reference GTF
#'
#' @description
#' Concatenates known GTF and newly predicted lncRNA GTF into
#' \code{reference/combined_mRNA_lncRNA.gtf} for downstream count module usage.
#'
#' @param cfg A validated prelnc configuration list.
#' @param final_lnc_gtf Path to final predicted lncRNA GTF.
#' @param ctx Runtime context.
#' @return Path to combined GTF file.
#' @export
build_augmented_reference_gtf <- function(cfg, final_lnc_gtf, ctx) {
  if (!file.exists(cfg$gtf_file)) stop("Known GTF missing: ", cfg$gtf_file)
  if (!file.exists(final_lnc_gtf) && !isTRUE(cfg$dry_run)) stop("final_lnc.gtf missing: ", final_lnc_gtf)

  known_lines <- readLines(cfg$gtf_file, warn = FALSE)
  lnc_lines <- if (file.exists(final_lnc_gtf)) readLines(final_lnc_gtf, warn = FALSE) else character(0)
  writeLines(c(known_lines, lnc_lines), ctx$files$combined_gtf)

  if (file.exists(final_lnc_gtf)) {
    has_exon <- any(grepl("\texon\t", lnc_lines, fixed = TRUE))
    if (!has_exon) {
      add_warning(ctx, "final_lnc.gtf has no exon feature lines; Cell Ranger GTF compatibility may be affected.")
    }
  }
  ctx$files$combined_gtf
}

run_prelnc_core_from_bam <- function(cfg, bam_df, ctx) {
  if (isTRUE(cfg$dry_run)) {
    for (i in seq_len(nrow(bam_df))) {
      sid <- sanitize_sample_id(bam_df$sample_id[i])
      out_gtf <- file.path(ctx$dirs$assembly_per_sample, paste0(sid, ".gtf"))
      args <- c("-p", as.character(cfg$threads), stringtie_strand_args(cfg$strandness), "-G", cfg$gtf_file, "-o", out_gtf, "-l", paste0(cfg$lncrna_name, "_", sid), bam_df$bam_path[i])
      run_command(ctx, paste0("stringtie_", sid), "stringtie", args, check = FALSE)
      add_status(ctx, "stringtie", bam_df$sample_id[i], "DRY_RUN", out_gtf)
    }
    run_command(ctx, "stringtie_merge", "stringtie", c("--merge", "-p", as.character(cfg$threads), stringtie_strand_args(cfg$strandness), "-G", cfg$gtf_file, "-o", file.path(ctx$dirs$assembly_root, "all.merged.gtf"), file.path(ctx$dirs$assembly_root, "gtf_list.txt")), check = FALSE)
    run_command(ctx, "gffcompare", "gffcompare", c("-r", cfg$gtf_file, "-o", file.path(ctx$dirs$gffcompare, "gffcmp"), file.path(ctx$dirs$assembly_root, "all.merged.gtf")), check = FALSE)
    run_command(ctx, "gffread_candidates", "gffread", c(file.path(ctx$dirs$gffcompare, "candidate_lnc.gtf"), "-g", cfg$genome_file, "-w", file.path(ctx$dirs$cpc2, "candidate_lncRNA.fa")), check = FALSE)
    run_command(ctx, "cpc2_predict", "CPC2.py", c("-i", file.path(ctx$dirs$cpc2, "candidate_lncRNA.fa"), "-o", file.path(ctx$dirs$cpc2, "CPC2.out"), "--ORF"), check = FALSE)
    run_command(ctx, "gffread_final_fasta", "gffread", c(ctx$files$final_gtf, "-g", cfg$genome_file, "-w", ctx$files$final_fa), check = FALSE)
    add_status(ctx, "cpc2_filter", "", "DRY_RUN", ctx$files$final_gtf)
    combined <- build_augmented_reference_gtf(cfg, ctx$files$final_gtf, ctx)
    add_status(ctx, "build_augmented_reference_gtf", "", "DRY_RUN", combined)
    return(invisible(list(final_gtf = ctx$files$final_gtf, final_fa = ctx$files$final_fa, combined_gtf = combined)))
  }

  gtf_files <- character(0)
  for (i in seq_len(nrow(bam_df))) {
    srow <- bam_df[i, , drop = FALSE]
    add_status(ctx, "stringtie", srow$sample_id, "RUNNING", "")
    st <- run_stringtie_for_sample(srow, cfg, ctx)
    add_status(ctx, "stringtie", srow$sample_id, "OK", st$gtf_path)
    gtf_files <- c(gtf_files, st$gtf_path)
  }
  merged_gtf <- merge_stringtie_gtfs(gtf_files, cfg, ctx)
  add_status(ctx, "stringtie_merge", "", "OK", merged_gtf)
  cand <- filter_lncRNA_candidates(merged_gtf, cfg, ctx)
  add_status(ctx, "gffcompare_filter", "", "OK", cand$candidate_gtf)
  final <- filter_coding_potential_with_CPC2(cand$candidate_gtf, cfg, ctx)
  add_status(ctx, "cpc2_filter", "", "OK", final$final_gtf)
  combined <- build_augmented_reference_gtf(cfg, final$final_gtf, ctx)
  add_status(ctx, "build_augmented_reference_gtf", "", "OK", combined)
  list(final_gtf = final$final_gtf, final_fa = final$final_fa, combined_gtf = combined)
}

#' Write prelnc markdown run report
#'
#' @param cfg A validated prelnc configuration list.
#' @param manifest FASTQ manifest data.frame.
#' @param statuses Per-step status table.
#' @param ctx Runtime context.
#' @param outputs Output path list from prelnc run.
#' @return Path to report file.
#' @export
write_prelnc_run_report <- function(cfg, manifest, statuses, ctx, outputs) {
  lines <- c(
    "# scLncR prelnc Run Report",
    "",
    sprintf("- Time: %s", now_string()),
    sprintf("- input_type: %s", cfg$input_type),
    sprintf("- sequencing_platform: %s", cfg$sequencing_platform),
    sprintf("- output_path: %s", cfg$output_path),
    sprintf("- dry_run: %s", cfg$dry_run),
    "",
    "## Workflow Summary",
    "",
    "- Raw FASTQ input discovery",
    "- Technology-aware discovery alignment",
    "- StringTie transcript assembly",
    "- gffcompare structural filtering",
    "- CPC2 coding potential filtering",
    "- Augmented reference construction (known GTF + final_lnc.gtf)",
    ""
  )

  lines <- c(lines, "## Input Manifest Preview", "", "```text")
  lines <- c(lines, capture.output(print(utils::head(manifest, 10), row.names = FALSE)))
  if (nrow(manifest) > 10) lines <- c(lines, sprintf("... (%d more rows)", nrow(manifest) - 10))
  lines <- c(lines, "```", "")

  lines <- c(lines, "## Technology Notes", "")
  lines <- c(lines, "- 10x mode uses R2 cDNA reads for candidate transcript evidence; R1 and I1 are retained for traceability.")
  lines <- c(lines, "- Standard 10x 3'/5' data do not provide uniform full-length transcript coverage.")
  lines <- c(lines, "- Predicted lncRNAs should be interpreted as candidate transcript evidence and cross-validated when possible.")
  lines <- c(lines, "- Smart-seq2 is generally more suitable for transcript reconstruction but still depends on read length/coverage/alignment quality.")
  lines <- c(lines, "- Drop-seq interface is currently experimental for full UMI-aware quantification.", "")

  lines <- c(lines, "## Output Files", "")
  lines <- c(lines, sprintf("- final_lnc.gtf: `%s`", outputs$final_gtf))
  lines <- c(lines, sprintf("- final.lncRNA.fa: `%s`", outputs$final_fa))
  lines <- c(lines, sprintf("- reference/combined_mRNA_lncRNA.gtf: `%s`", outputs$combined_gtf))
  if (!is.null(outputs$id_mapping_file)) {
    lines <- c(lines, sprintf("- final lncRNA ID mapping: `%s`", outputs$id_mapping_file))
  }
  if (!is.null(outputs$cpc_summary_file)) {
    lines <- c(lines, sprintf("- CPC2 filter summary: `%s`", outputs$cpc_summary_file))
  }
  lines <- c(lines, sprintf("- raw FASTQ manifest: `%s`", ctx$files$raw_fastq_manifest))
  lines <- c(lines, sprintf("- commands log: `%s`", ctx$files$commands_log))
  lines <- c(lines, sprintf("- run log: `%s`", ctx$files$run_log), "")

  lines <- c(lines, "## Step Status", "", "```text")
  if (nrow(statuses) > 0) {
    lines <- c(lines, capture.output(print(statuses, row.names = FALSE)))
  } else {
    lines <- c(lines, "No step status rows.")
  }
  lines <- c(lines, "```", "")

  lines <- c(lines, "## Warnings", "")
  if (length(ctx$warnings) == 0) {
    lines <- c(lines, "- None.", "")
  } else {
    lines <- c(lines, paste0("- ", ctx$warnings), "")
  }

  lines <- c(lines, "## Config Snapshot", "", "```yaml")
  lines <- c(lines, strsplit(yaml::as.yaml(cfg), "\n", fixed = TRUE)[[1]])
  lines <- c(lines, "```", "")

  lines <- c(lines, "## Command Preview", "", "```text")
  if (length(ctx$command_preview) > 0) {
    lines <- c(lines, ctx$command_preview)
  } else {
    lines <- c(lines, "No commands recorded.")
  }
  lines <- c(lines, "```", "")

  writeLines(lines, ctx$files$report)
  ctx$files$report
}

check_external_tools <- function(cfg) {
  needed <- c("hisat2", "hisat2-build", "samtools", "stringtie", "gffcompare", "gffread", "CPC2.py")
  miss <- needed[Sys.which(needed) == ""]
  if (length(miss) == 0) return(invisible(TRUE))
  if (isTRUE(cfg$dry_run)) {
    warning("Missing external tools in dry_run mode: ", paste(miss, collapse = ", "))
  } else {
    stop("Missing required external tools: ", paste(miss, collapse = ", "))
  }
}

run_prelnc_fastq_first <- function(cfg, ctx) {
  manifest <- discover_raw_fastq_inputs(cfg)
  write.table(manifest, ctx$files$raw_fastq_manifest, sep = "\t", row.names = FALSE, quote = FALSE)

  val <- validate_raw_fastq_manifest(manifest, cfg)
  for (w in val$warnings) add_warning(ctx, w)
  write_validation_report(ctx, val)
  if (length(val$errors) > 0) {
    stop("FASTQ manifest validation failed:\n- ", paste(val$errors, collapse = "\n- "))
  }
  copy_if_exists(cfg$samples_info, file.path(ctx$dirs$manifest, "samples_info.copy.tsv"))

  idx <- build_hisat2_index_if_needed(cfg, ctx)
  bam_rows <- data.frame(sample_id = character(0), bam_path = character(0), stringsAsFactors = FALSE)
  for (i in seq_len(nrow(manifest))) {
    srow <- manifest[i, , drop = FALSE]
    add_status(ctx, "discovery_alignment", srow$sample_id, "RUNNING", "")
    sam <- run_discovery_alignment_for_sample(srow, cfg, ctx, idx)
    add_status(ctx, "discovery_alignment", srow$sample_id, "OK", sam$sam_path)
    add_status(ctx, "sam_to_sorted_bam", srow$sample_id, "RUNNING", "")
    bam <- sam_to_sorted_bam(sam, cfg, ctx)
    add_status(ctx, "sam_to_sorted_bam", srow$sample_id, "OK", bam$bam_path)
    bam_rows <- rbind(bam_rows, data.frame(sample_id = bam$sample_id, bam_path = bam$bam_path, stringsAsFactors = FALSE))
  }
  outputs <- run_prelnc_core_from_bam(cfg, bam_rows, ctx)
  list(manifest = manifest, outputs = outputs)
}

print_prelnc_help <- function() {
  cat("scLncR prelnc - raw FASTQ-first candidate lncRNA discovery\n")
  cat("Usage: scLncR prelnc -c <config.yaml>\n\n")
  cat("Design:\n")
  cat("  - User-facing input is raw FASTQ (technology-aware)\n")
  cat("  - Internal BAM is used as intermediate for transcript evidence\n")
  cat("  - Output final_lnc.gtf and reference/combined_mRNA_lncRNA.gtf for count module\n\n")
  cat("10x note:\n")
  cat("  - I1: sample index, R1: barcode/UMI, R2: cDNA\n")
  cat("  - prelnc discovery uses R2 for candidate transcript evidence\n")
}

#' Backward-compatible wrapper for legacy prelnc function
#'
#' @description
#' Legacy wrapper retained for compatibility. It maps old R2 FASTQ inputs to the
#' current raw FASTQ-first prelnc workflow.
#'
#' @param genome_file Reference genome FASTA.
#' @param genome_gtf Reference GTF.
#' @param sample_files Legacy FASTQ list (historically R2 files).
#' @param lncRNA_name lncRNA prefix.
#' @param threads Thread number.
#' @param output_file Legacy output argument (kept for API compatibility).
#' @param merge_name Legacy project name argument.
#' @return List with final output paths.
#' @export
get_gtf <- function(genome_file, genome_gtf, sample_files, lncRNA_name = "", threads = 4, output_file = "", merge_name = "") {
  warning("get_gtf() is deprecated; use run_prelnc() raw FASTQ-first workflow instead.")
  if (length(sample_files) == 0) stop("Legacy get_gtf() requires non-empty sample_files.")

  tmp_fastq <- tempfile("prelnc_legacy_fastq_")
  dir.create(tmp_fastq, recursive = TRUE, showWarnings = FALSE)
  for (f in sample_files) {
    if (!file.exists(f)) stop("sample file not found: ", f)
    file.link(f, file.path(tmp_fastq, basename(f)))
  }

  cfg <- list(
    samples_dirs = tmp_fastq,
    samples_info = "",
    project_name = ifelse(nzchar(merge_name), merge_name, "legacy_prelnc"),
    output_path = getwd(),
    genome_file = genome_file,
    gtf_file = genome_gtf,
    lncrna_name = ifelse(nzchar(lncRNA_name), lncRNA_name, "Lnc"),
    threads = as.integer(threads),
    input_type = "fastq",
    sequencing_platform = "10x",
    read_layout = "auto",
    fastq_sample_regex = "^(.*?)_S\\d+_L\\d+_([IR][12])_\\d+\\.(fastq|fq)(\\.gz)?$",
    r1_pattern = "_R1_",
    r2_pattern = "_R2_",
    i1_pattern = "_I1_",
    aligner_for_discovery = "hisat2",
    assembler = "stringtie",
    strandness = "RF",
    min_transcript_length = 200L,
    gffcompare_classes = c("i", "u", "x", "o"),
    min_mapping_quality = 30L,
    rebuild_hisat2_index = FALSE,
    hisat2_index_prefix = "",
    keep_intermediate = TRUE,
    dry_run = FALSE,
    sample_limit = 0L,
    include_samples = character(0)
  )
  ctx <- init_prelnc_context(cfg)
  check_external_tools(cfg)
  res <- run_prelnc_fastq_first(cfg, ctx)
  write_prelnc_run_report(cfg, res$manifest, ctx$status_rows, ctx, res$outputs)
  res$outputs
}

#' CLI entry for prelnc module
#'
#' @param user_args CLI arguments for prelnc.
#' @param script_dir scLncR installation root.
#' @return Invisibly returns output paths.
#' @export
run_prelnc <- function(user_args, script_dir) {
  option_list <- list(
    optparse::make_option(c("-c", "--config"), type = "character", metavar = "FILE", help = "Path to prelnc YAML config")
  )
  parser <- optparse::OptionParser(
    usage = "scLncR prelnc [options]",
    option_list = option_list,
    description = "prelnc: raw FASTQ-first, technology-aware candidate lncRNA discovery"
  )
  opt <- optparse::parse_args(parser, args = user_args, print_help_and_exit = FALSE)

  if (length(user_args) == 0 || any(user_args %in% c("-h", "--help")) || is.null(opt$config)) {
    print_prelnc_help()
    cat("\n")
    optparse::print_help(parser)
    return(invisible(NULL))
  }
  if (!file.exists(opt$config)) stop("Config file not found: ", opt$config)

  cfg <- load_prelnc_config(opt$config)
  ctx <- init_prelnc_context(cfg)
  copy_if_exists(opt$config, file.path(ctx$dirs$manifest, "config_LncPre.snapshot.yaml"))
  copy_if_exists(cfg$samples_info, file.path(ctx$dirs$manifest, "samples_info.copy.tsv"))

  log_prelnc(ctx, "Starting prelnc raw FASTQ-first workflow.")
  log_prelnc(ctx, sprintf("sequencing_platform=%s | input_type=%s | dry_run=%s", cfg$sequencing_platform, cfg$input_type, cfg$dry_run))
  res <- run_prelnc_fastq_first(cfg, ctx)
  write_prelnc_run_report(cfg, res$manifest, ctx$status_rows, ctx, res$outputs)
  log_prelnc(ctx, sprintf("prelnc completed. final_lnc.gtf: %s", res$outputs$final_gtf))

  invisible(list(
    final_lnc_gtf = res$outputs$final_gtf,
    final_lnc_fasta = res$outputs$final_fa,
    combined_gtf = res$outputs$combined_gtf,
    report = ctx$files$report
  ))
}

# Backward-compatibility aliases for previous helper names
discover_fastq_inputs <- discover_raw_fastq_inputs
validate_fastq_manifest <- validate_raw_fastq_manifest
run_hisat2_for_sample <- run_discovery_alignment_for_sample
