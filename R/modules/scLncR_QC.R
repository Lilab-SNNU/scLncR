suppressPackageStartupMessages({
  suppressWarnings({
    suppressMessages({
      library(optparse)
    })
  })
})

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0) return(y)
  x
}

qc_now <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S")

qc_append <- function(path, txt) {
  cat(txt, file = path, append = TRUE, sep = "\n")
}

parse_extra_args <- function(x) {
  if (is.null(x) || !nzchar(trimws(x))) return(character(0))
  scan(text = x, what = character(), quiet = TRUE)
}

#' Discover FASTQ files for raw FASTQ QC
#'
#' @description
#' Searches the configured input directory for FASTQ files and applies optional
#' include_files and sample_limit filters. This function only discovers files;
#' it does not read or modify FASTQ content.
#'
#' @param cfg A validated QC configuration list.
#' @return Character vector of FASTQ file paths.
#' @export
discover_qc_fastq_files <- function(cfg) {
  if (isTRUE(cfg$recursive)) {
    all_files <- list.files(cfg$input_dir, recursive = TRUE, full.names = TRUE)
    fastq_files <- all_files[grepl(utils::glob2rx(cfg$file_pattern), basename(all_files))]
  } else {
    fastq_files <- Sys.glob(file.path(cfg$input_dir, cfg$file_pattern))
  }
  fastq_files <- sort(normalizePath(fastq_files, winslash = "/", mustWork = FALSE))

  if (length(cfg$include_files) > 0) {
    keep <- basename(fastq_files) %in% cfg$include_files | fastq_files %in% cfg$include_files
    fastq_files <- fastq_files[keep]
  }
  if (cfg$sample_limit > 0) {
    fastq_files <- fastq_files[seq_len(min(cfg$sample_limit, length(fastq_files)))]
  }
  if (length(fastq_files) == 0) {
    stop("No FASTQ files matched QC configuration.")
  }
  fastq_files
}

init_qc_context <- function(cfg) {
  dirs <- list(
    base = cfg$output_dir,
    fastqc = file.path(cfg$output_dir, "fastqc"),
    multiqc = file.path(cfg$output_dir, "multiqc"),
    logs = file.path(cfg$output_dir, "logs")
  )
  for (d in unlist(dirs, use.names = FALSE)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

  files <- list(
    run_log = file.path(dirs$logs, "qc_run.log"),
    commands_log = file.path(dirs$logs, "commands.log"),
    fastq_file_list = file.path(dirs$logs, "fastq_file_list.txt"),
    summary = file.path(dirs$base, "qc_summary.md")
  )
  writeLines(character(0), files$run_log)
  writeLines(character(0), files$commands_log)

  list(cfg = cfg, dirs = dirs, files = files, start_time = Sys.time(), commands = character(0))
}

log_qc <- function(ctx, msg, level = "INFO") {
  line <- sprintf("[%s] [%s] %s", qc_now(), level, msg)
  qc_append(ctx$files$run_log, line)
  message(line)
}

record_qc_command <- function(ctx, step, command, args) {
  cmd <- paste(c(command, shQuote(args)), collapse = " ")
  qc_append(ctx$files$commands_log, sprintf("[%s] [%s] %s", qc_now(), step, cmd))
  ctx$commands <- c(ctx$commands, sprintf("[%s] %s", step, cmd))
  invisible(ctx)
}

run_qc_command <- function(ctx, step, command, args) {
  record_qc_command(ctx, step, command, args)
  if (isTRUE(ctx$cfg$dry_run)) {
    log_qc(ctx, sprintf("[dry_run] %s", paste(c(command, shQuote(args)), collapse = " ")))
    return(invisible(ctx))
  }

  out <- tryCatch(
    system2(command, args = args, stdout = TRUE, stderr = TRUE),
    error = function(e) {
      structure(e$message, status = 127L)
    }
  )
  if (length(out) > 0) qc_append(ctx$files$run_log, paste(out, collapse = "\n"))
  status <- attr(out, "status") %||% 0L
  if (!identical(as.integer(status), 0L)) {
    stop(sprintf("QC command failed at step '%s' with exit status %s.", step, status))
  }
  invisible(ctx)
}

build_fastqc_args <- function(cfg, fastq_files, fastqc_dir) {
  args <- c("-t", as.character(cfg$threads), "-o", fastqc_dir)
  if (isTRUE(cfg$fastqc$extract)) args <- c(args, "--extract")
  if (isTRUE(cfg$fastqc$nogroup)) args <- c(args, "--nogroup")
  if (nzchar(cfg$fastqc$contaminants)) args <- c(args, "--contaminants", cfg$fastqc$contaminants)
  if (nzchar(cfg$fastqc$adapters)) args <- c(args, "--adapters", cfg$fastqc$adapters)
  c(args, parse_extra_args(cfg$fastqc$extra_args), fastq_files)
}

build_multiqc_args <- function(cfg, fastqc_dir, multiqc_dir) {
  args <- c(fastqc_dir, "-o", multiqc_dir, "--title", cfg$multiqc$title, "--filename", cfg$multiqc$filename)
  c(args, parse_extra_args(cfg$multiqc$extra_args))
}

summarize_fastqc_summary_files <- function(fastqc_dir) {
  summary_files <- list.files(fastqc_dir, pattern = "^summary\\.txt$", recursive = TRUE, full.names = TRUE)
  if (length(summary_files) == 0) return(data.frame())

  rows <- lapply(summary_files, function(path) {
    tab <- tryCatch(read.delim(path, header = FALSE, sep = "\t", stringsAsFactors = FALSE), error = function(e) data.frame())
    if (nrow(tab) == 0 || ncol(tab) < 3) return(NULL)
    sample_name <- tab[1, 3]
    fail_modules <- tab[tab[[1]] == "FAIL", 2]
    data.frame(
      sample = sample_name,
      pass = sum(tab[[1]] == "PASS"),
      warn = sum(tab[[1]] == "WARN"),
      fail = sum(tab[[1]] == "FAIL"),
      fail_modules = if (length(fail_modules) == 0) "" else paste(fail_modules, collapse = "; "),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows[!vapply(rows, is.null, logical(1))])
}

write_qc_summary <- function(ctx, fastq_files, fastqc_summary = data.frame()) {
  cfg <- ctx$cfg
  end_time <- Sys.time()
  multiqc_report <- file.path(ctx$dirs$multiqc, cfg$multiqc$filename)

  lines <- c(
    "# scLncR Raw FASTQ QC Summary",
    "",
    sprintf("- Start time: %s", format(ctx$start_time, "%Y-%m-%d %H:%M:%S")),
    sprintf("- End time: %s", format(end_time, "%Y-%m-%d %H:%M:%S")),
    sprintf("- Input directory: `%s`", cfg$input_dir),
    sprintf("- Output directory: `%s`", cfg$output_dir),
    sprintf("- FASTQ files: %d", length(fastq_files)),
    sprintf("- dry_run: %s", cfg$dry_run),
    sprintf("- threads: %s", cfg$threads),
    sprintf("- FastQC enabled: %s", cfg$fastqc$enabled),
    sprintf("- MultiQC enabled: %s", cfg$multiqc$enabled),
    sprintf("- FastQC output directory: `%s`", ctx$dirs$fastqc),
    sprintf("- MultiQC report: `%s`", multiqc_report),
    "",
    "## Interpretation Notes",
    "",
    "- This optional module only reports raw FASTQ quality metrics.",
    "- It does not perform trimming, adapter removal, read filtering, or any modification of the original FASTQ files.",
    "- Review the MultiQC HTML report before deciding whether to continue into prelnc/count workflows.",
    ""
  )

  if (nrow(fastqc_summary) > 0) {
    lines <- c(lines, "## FastQC Summary", "", "```text")
    lines <- c(lines, capture.output(print(fastqc_summary, row.names = FALSE)))
    lines <- c(lines, "```", "")
  } else {
    lines <- c(lines, "## FastQC Summary", "", "- No extracted FastQC summary.txt files were parsed.", "")
  }

  lines <- c(lines, "## Commands", "", "```text")
  commands <- readLines(ctx$files$commands_log, warn = FALSE)
  lines <- c(lines, if (length(commands) == 0) "No commands recorded." else commands)
  lines <- c(lines, "```", "")

  writeLines(lines, ctx$files$summary)
  invisible(ctx$files$summary)
}

#' Run raw FASTQ quality control for scLncR
#'
#' @description
#' Runs FastQC on raw FASTQ files and summarizes the results using MultiQC.
#' This optional module only reports quality metrics and does not modify,
#' trim, or filter input reads.
#'
#' @param cfg A validated QC configuration list.
#' @return Invisibly returns the QC output directory.
#' @export
run_fastq_qc <- function(cfg) {
  validate_qc_tools(cfg)
  ctx <- init_qc_context(cfg)
  log_qc(ctx, "Starting raw FASTQ QC.")

  fastq_files <- discover_qc_fastq_files(cfg)
  writeLines(fastq_files, ctx$files$fastq_file_list)
  log_qc(ctx, sprintf("FASTQ files selected: %d", length(fastq_files)))

  if (isTRUE(cfg$fastqc$enabled)) {
    args <- build_fastqc_args(cfg, fastq_files, ctx$dirs$fastqc)
    ctx <- run_qc_command(ctx, "fastqc", "fastqc", args)
  } else {
    log_qc(ctx, "FastQC disabled by config.")
  }

  if (isTRUE(cfg$multiqc$enabled)) {
    args <- build_multiqc_args(cfg, ctx$dirs$fastqc, ctx$dirs$multiqc)
    ctx <- run_qc_command(ctx, "multiqc", "multiqc", args)
  } else {
    log_qc(ctx, "MultiQC disabled by config.")
  }

  fastqc_summary <- if (!isTRUE(cfg$dry_run) && isTRUE(cfg$fastqc$extract)) {
    summarize_fastqc_summary_files(ctx$dirs$fastqc)
  } else {
    data.frame()
  }
  write_qc_summary(ctx, fastq_files, fastqc_summary)
  log_qc(ctx, sprintf("QC summary written: %s", ctx$files$summary))
  invisible(cfg$output_dir)
}

print_qc_help <- function() {
  cat("scLncR qc - raw FASTQ quality control\n")
  cat("Usage: scLncR qc -c <config.yaml>\n\n")
  cat("Description:\n")
  cat("  - Runs FastQC and MultiQC on raw FASTQ files\n")
  cat("  - Does not trim, filter, or modify input FASTQ files\n")
  cat("  - Writes fastqc/, multiqc/, logs/, and qc_summary.md\n")
}

#' CLI entry for the optional raw FASTQ QC module
#'
#' @param user_args CLI arguments for qc.
#' @param script_dir scLncR installation root.
#' @return Invisibly returns the QC output directory.
#' @export
run_qc <- function(user_args, script_dir) {
  option_list <- list(
    optparse::make_option(c("-c", "--config"), type = "character", metavar = "FILE", help = "Path to QC YAML config")
  )
  parser <- optparse::OptionParser(
    usage = "scLncR qc [options]",
    option_list = option_list,
    description = "qc: optional raw FASTQ quality control using FastQC and MultiQC"
  )
  opt <- optparse::parse_args(parser, args = user_args, print_help_and_exit = FALSE)

  if (length(user_args) == 0 || any(user_args %in% c("-h", "--help")) || is.null(opt$config)) {
    print_qc_help()
    cat("\n")
    optparse::print_help(parser)
    return(invisible(NULL))
  }

  cfg <- load_qc_config(opt$config)
  run_fastq_qc(cfg)
}
