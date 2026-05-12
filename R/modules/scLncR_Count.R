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
    report = file.path(dirs$base, "count_run_report.md"),
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
  lines <- c(lines, "- Smart-seq2 interface is experimental in this version.")
  lines <- c(lines, "- Drop-seq interface is experimental and intended for future STARsolo/Drop-seq-tools integration.")
  lines <- c(lines, ""
  )

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
  cat("Experimental interfaces:\n")
  cat("  - smartseq2 + featurecounts\n")
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
  manifest <- discover_count_fastq_manifest(cfg, ctx)
  combined_gtf <- build_combined_gtf_for_count(cfg, ctx)
  add_count_status(ctx, "build_combined_gtf", "", "OK", combined_gtf)

  transcriptome_path <- NA_character_
  if (cfg$sequencing_platform == "10x" && cfg$count_engine == "cellranger") {
    add_count_status(ctx, "cellranger_mkref", "", "RUNNING", "")
    transcriptome_path <- build_cellranger_reference(cfg, ctx, combined_gtf)
    add_count_status(ctx, "cellranger_mkref", "", if (cfg$dry_run) "DRY_RUN" else "OK", transcriptome_path)

    for (i in seq_len(nrow(manifest))) {
      sid <- manifest$sample_id[i]
      add_count_status(ctx, "cellranger_count", sid, "RUNNING", "")
      run_cellranger_count_for_sample(sid, cfg, ctx, transcriptome_path)
      add_count_status(ctx, "cellranger_count", sid, if (cfg$dry_run) "DRY_RUN" else "OK", "")
    }
  } else {
    msg <- sprintf("Interface %s + %s is experimental and not fully implemented in this version.", cfg$sequencing_platform, cfg$count_engine)
    add_count_warning(ctx, msg)
    add_count_status(ctx, "experimental_interface", "", if (cfg$dry_run) "DRY_RUN" else "SKIPPED", msg)
    if (!cfg$dry_run) {
      stop(msg, " Use dry_run=true for interface validation/report only.")
    }
  }

  write_count_report(cfg, ctx, manifest, combined_gtf, transcriptome_path)
  log_count(ctx, "count workflow completed.")
  invisible(list(report = ctx$files$report, combined_gtf = combined_gtf))
}
