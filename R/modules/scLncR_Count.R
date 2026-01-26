suppressPackageStartupMessages({
    suppressWarnings({
        suppressMessages({
            library(optparse)
            library(stringr)
            library(gridExtra) 
            library(psych)
            library(dplyr)
            library(tidyverse)
            library(patchwork)
            library(reshape2)
            library(yaml)
        })
    })
})


run_mkgtf <- function(lncgtf, genegtf, project_name, output_path="./"){

    merge_gtf <- sprintf("cat %s %s > tmp.gtf", lncgtf, genegtf)
    mkgtf_command <- sprintf("cellranger mkgtf tmp.gtf %s/%s.gtf", output_path, project_name)
    
    system(merge_gtf)
    system(mkgtf_command)

    gtf_file <- sprintf("%s/%s.gtf", output_path, project_name)
    return(gtf_file)
}

run_mkref <- function(genome_file, gtf_file, project_name, output_path){

    mkref_command <- sprintf("cellranger mkref --fasta %s --genes %s --genome %s --output-dir %s/%s --nthreads 4", 
                             genome_file, gtf_file, project_name, output_path, project_name)
    print(mkref_command)
    system(mkref_command)
    trans_name <- sprintf("%s/%s", output_path, project_name)
    print(trans_name)

    return(trans_name)
}

run_cellranger <- function(samples, samples_dir, trans_name, cores_num=4, output_path="./"){

    for(sample in samples){
        print(sample)
        cellranger_count_command <- sprintf("cellranger count --id=%s --fastqs=%s --sample=%s --transcriptome=%s --localcores=%s --output-dir=%s/%s", 
                                            sample, samples_dir, sample, trans_name, cores_num, output_path, sample)
        print(cellranger_count_command)
        system(cellranger_count_command)
    }
}

# =============================================================================
# 
# Description: scRNA-seq express count Para Parse and Run Function
# 
# =============================================================================

run_count <- function(user_args, script_dir) {
  # Load required package
  if (!requireNamespace("optparse", quietly = TRUE)) {
    stop("Package 'optparse' is required. Please install it: install.packages('optparse')")
  }

  # Define option list
  option_list <- list(
    optparse::make_option(
      c("-c", "--config"),
      type = "character",
      help = "Path to the YAML configuration file (e.g., config.yaml)",
      metavar = "FILE"
    )
  )

  # Create parser
  parser <- optparse::OptionParser(
    usage = "scLncR count [options]",
    option_list = option_list,
    description = "count: Build reference and run CellRanger on sc/snRNA-seq data"
  )

  # Parse arguments
  opt <- optparse::parse_args(parser, args = user_args, print_help_and_exit = FALSE)

  # If config missing, show help and quit
  if (is.null(opt$config)) {
    cat("\n")
    optparse::print_help(parser)
    q()
  }

  config_file <- opt$config
  if (!file.exists(config_file)) {
    stop("Config file not found: ", config_file)
  }
  
  # Load and validate config
  # source(file.path(script_dir, "R", "utils", "config_loader_count.R"))
  cfg <- load_count_config(config_file)

  # Extract parameters
  genome_file    <- cfg$genome
  gtf_file       <- cfg$gtf
  lnc_gtf        <- cfg$lnc_gtf
  samples_dir    <- cfg$samples_dir
  output_path    <- cfg$output_path
  project_name   <- cfg$project_name
  threads        <- cfg$threads

  ####################***Work Dir Set***####################
  dir.create(output_path, recursive = TRUE)
  setwd(output_path)

  #########################***Auto-detect Samples from FASTQ***#########################
  
  # List all .fastq.gz files in samples_dir
  fastq_files <- list.files(
    path = samples_dir,
    pattern = "\\.fastq\\.gz$",
    full.names = FALSE,
    recursive = FALSE
  )
  
  if (length(fastq_files) == 0) {
    stop("No .fastq.gz files found in samples_dir: ", samples_dir)
  }
  
  # Extract sample name: everything before the first "_"
  samples <- sort(unique(sapply(strsplit(fastq_files, "_"), `[`, 1)))
  
  cat("Auto-detected samples from FASTQ filenames:\n")
  cat("  ", paste(samples, collapse = ", "), "\n\n")

  #########################***Pipeline Function Run***#########################

  ## Step 1: Merge GTF
  gtf_file_merged <- run_mkgtf(
    lncgtf       = lnc_gtf,
    genegtf      = gtf_file,
    project_name = project_name,
    output_path  = output_path
  )
  cat("####################***Function mkgtf done***####################\n")

  ## Step 2: Build Reference
  trans_name <- run_mkref(
    genome_file  = genome_file,
    gtf_file     = gtf_file_merged,
    project_name = project_name,
    output_path  = output_path
  )
  cat("####################***Function mkref done***####################\n")

  ## Step 3: Run CellRanger Count
  run_cellranger(
    samples      = samples,
    samples_dir  = samples_dir,
    trans_name   = trans_name,
    cores_num    = 4,
    output_path  = output_path
  )
  cat("####################***Function count done***####################\n")
}

