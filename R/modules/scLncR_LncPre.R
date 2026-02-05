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

# =============================================================================
# Description: LncRNA Predict main Function
# Paras: Same as Config File
# =============================================================================

get_gtf <- function(genome_file, genome_gtf, sample_files, lncRNA_name="", threads=4, output_file="", merge_name=""){

    path <- getwd()
    dir.create("./step1_res")
    system("mkdir ./step1_res/hisat2 ./step1_res/stringtie res")
    genome_index <- strsplit(genome_file, ".fa")[[1]][1]
    
    #===============***Step 1: Chain specific alignment and assembly of transcripts***===============#
    
    ## Hisat2 build index
    print("#####***Hisat2 build index***#####")
    hisat2_build_command <- sprintf("hisat2-build %s %s -q >> ./step1_res/log_file.txt 2>&1", genome_file, genome_index)
    system(hisat2_build_command)
    samples <- c()

    ## Sample trans gtf get
    for(sample_file in sample_files){
        
        infos <- strsplit(sample_file, "/")[[1]]
        num <- length(infos)
        sample <- strsplit(infos[num], "_")[[1]][1]
        sample_dir <- paste(infos[1:num-1], collapse="/")
        samples <- append(samples, sample)
        
        print(sprintf("#####***Hisat2 aligment and Stringtie assembly of %s***#####", sample))
        ## Hisat2 aligment
        hisat2_aligment_command <- sprintf("hisat2 --dta -x %s -p %s --rna-strandness RF -U %s -S %s/step1_res/hisat2/lnc.%s.sam  >> ./step1_res/log_file.txt 2>&1", genome_index, threads, sample_file, path, sample) 
        system(hisat2_aligment_command)
        
        ## Samtools deal the result from hiast2
        samtools_command <- sprintf("samtools view -@ %s -bS -q 30 ./step1_res/hisat2/lnc.%s.sam | samtools sort -@ 4 - -o ./step1_res/hisat2/lnc.%s.bam >> ./step1_res/log_file.txt 2>&1", threads, sample, sample)
        system(samtools_command)
        system(sprintf("rm ./step1_res/hisat2/lnc.%s.sam", sample))
        
        ## Stringtie assembly
        stringtie_command <- sprintf("stringtie -p %s --rf -G %s -o ./step1_res/stringtie/lnc.%s.gtf -l %s_lnc ./step1_res/hisat2/lnc.%s.bam >> ./step1_res/log_file.txt 2>&1", threads, genome_gtf, sample, sample, sample)
        system(stringtie_command)
    }

    #======================***Step 2: Filter lncRNA by the characteristics of lncRNA***======================#
    
    print("#####***LncRNA Prediction***#####")
    ## Merge all lnc_gtf
    system("ls ./step1_res/stringtie/*.gtf > ./step1_res/gtf_list.txt")
    if(merge_name==""){
        merge_name <- sample
    }else{
        merge_name <- merge_name
    }
    print("***LncRNA Info Merge***")
    stringtie_merge <- sprintf("stringtie --merge -p %s --rf -G %s -o ./step1_res/stringtie/all.merged.gtf -l %s ./step1_res/gtf_list.txt >> ./step1_res/log_file.txt 2>&1", threads, genome_gtf, lncRNA_name)
    system(stringtie_merge)
    
    print("***LncRNA Filter by Bio Characteristics***")
    ## Compare trans gtf with genome gtf
    gffcompare_command1 <- sprintf("gffcompare -r %s -o ./step1_res/stringtie/gffcmp ./step1_res/stringtie/all.merged.gtf", genome_gtf)
    system(gffcompare_command1)
    
    system("awk '$3==\"i\" || $3==\"u\" || $3==\"x\" || $3==\"o\" && $10>200 {print \"\\\"\" $5 \"\\\"\"}' ./step1_res/stringtie/gffcmp.all.merged.gtf.tmap > ./step1_res/stringtie/filter.iuxo.l200.transcript.id.txt")
    system("LC_ALL=C fgrep -f ./step1_res/stringtie/filter.iuxo.l200.transcript.id.txt ./step1_res/stringtie/gffcmp.annotated.gtf > ./step1_res/stringtie/filter.iuxo.l200.gtf")
    system(sprintf("gffread ./step1_res/stringtie/filter.iuxo.l200.gtf -g %s -w ./step1_res/stringtie/filter.iuxo.l200.transcript.fa", genome_file))
    
    #===============***Step 3: Filter lncRNA filtering lncRNA through encoding function***===============#
    
    ## Filter out transcripts with protein encoding ability
    print("***LncRNA Filter by Coding Protential***")

    CPC2_command <- "CPC2.py -i ./step1_res/stringtie/filter.iuxo.l200.transcript.fa -o ./step1_res/stringtie/CPC2.out --ORF"
    system(CPC2_command)
    system("awk '$9==\"noncoding\"{print \"\\\"\" $1 \"\\\"\"}' ./step1_res/stringtie/CPC2.out.txt > ./step1_res/stringtie/tmp.CPC2.txt")
    system("cat ./step1_res/stringtie/tmp.CPC2.txt | sort > ./step1_res/stringtie/filter.iuxo.l200.ncd.transcript.txt")
    system(sprintf("LC_ALL=C fgrep -f ./step1_res/stringtie/filter.iuxo.l200.ncd.transcript.txt ./step1_res/stringtie/filter.iuxo.l200.gtf > %s", output_file))
    system(sprintf("gffread %s -g %s -w ./res/final.lncRNA.fa", output_file, genome_file))

    res <- c(samples, output_file, sample_dir)
    return (res)
}

# =============================================================================
# 
# Description: LncRNA Predict Para Parse and Run Function
# 
# =============================================================================


run_prelnc <- function(user_args, script_dir) {
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
    usage = "scLncR prelnc [options]",
    option_list = option_list,
    description = "prelnc: Predict lncRNAs from single-cell/single-nucleus RNA-seq data"
  )

  # Parse arguments
  opt <- optparse::parse_args(parser, args = user_args, print_help_and_exit = FALSE)

  # If help requested or config missing
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
  cfg <- load_prelnc_config(config_file)

  # Extract parameters (optional: you can just pass cfg directly to functions)
  samples_dirs <- cfg$samples_dirs
  genome_file  <- cfg$genome_file
  gtf_file     <- cfg$gtf_file
  project_name <- cfg$project_name
  output_path  <- cfg$output_path
  threads      <- cfg$threads
  lncrna_name  <- cfg$lncrna_name

  ####################***Work Dir Set***####################
  dir.create(output_path, recursive = TRUE)
  setwd(output_path)

  #########################***Pipeline Function Run***#########################

  ##########***Parameters Get and Run***##########

  all_files <- list.files(samples_dirs)
  sample_files <- all_files[grep("R2_", all_files)]
  samples <- str_split_fixed(sample_files, "_S", 2)[,1]
  samples_file <- paste0(samples_dirs, sample_files)

  ## Run LncRNA Predict Function

  res_info <- get_gtf(genome_file=genome_file,
                      genome_gtf=gtf_file,
                      sample_files=samples_file,
                      lncRNA_name=lncrna_name,
                      merge_name=project_name,
                      threads=threads,
                      output_file="./final_lnc.gtf")
  print("####################***Function LncRNA Predict done***####################")
}


