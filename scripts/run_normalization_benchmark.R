#!/usr/bin/env Rscript

parse_args_base <- function(args) {
  config_path <- NULL
  if (length(args) == 0 || any(args %in% c("-h", "--help"))) {
    cat("Usage:\n")
    cat("  Rscript scripts/run_normalization_benchmark.R -c <config.yaml>\n\n")
    cat("Options:\n")
    cat("  -c, --config   Path to normalization benchmark config YAML\n")
    cat("  -h, --help     Show this help message\n")
    quit(save = "no", status = 0)
  }

  i <- 1
  while (i <= length(args)) {
    token <- args[i]
    if (token %in% c("-c", "--config")) {
      if (i == length(args)) stop("Missing value for ", token)
      config_path <- args[i + 1]
      i <- i + 2
      next
    }
    i <- i + 1
  }
  if (is.null(config_path)) {
    stop("Configuration file must be provided via -c or --config")
  }
  config_path
}

args <- commandArgs(trailingOnly = TRUE)
config_path <- parse_args_base(args)

file_arg <- grep("^--file=", commandArgs(), value = TRUE)
if (length(file_arg) > 0) {
  script_path <- normalizePath(sub("^--file=", "", file_arg[1]), mustWork = FALSE)
} else {
  script_path <- normalizePath("scripts/run_normalization_benchmark.R", mustWork = FALSE)
}
repo_root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = FALSE)
package_loader <- file.path(repo_root, "R", "modules", "scLncR_load_package.R")
if (file.exists(package_loader)) {
  source(package_loader, local = TRUE)
}
module_file <- file.path(repo_root, "R", "modules", "scLncR_normalization_benchmark.R")
if (!file.exists(module_file)) {
  stop("Module file not found: ", module_file)
}
source(module_file, local = TRUE)

cfg <- read_normalization_benchmark_config(config_path)
if (exists("scLncR_load_normalization_benchmark_packages", mode = "function")) {
  scLncR_load_normalization_benchmark_packages(cfg$strategies)
}
run_normalization_benchmark(cfg)
cat("Normalization benchmark completed.\n")
