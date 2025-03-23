library(getopt)


spec <- matrix(
    c("genome",  "G", 2, "character", "Genome file , format is fasta!",
      "gtf", "g", 2, "character",  "GTF annotation file, need to match your genome file version!",
      "lnc_gtf",  "l", 2, "character",  "The gtf file of lncRNA",
      "samples",  "s", 2, "character",  "The single cell seq samples.",
      "samples_dir",  "S", 2, "character",  "The dir of your samples",
      "output_path",  "o", 2, "character",  "Output dir of res",
      "project_name",  "p", 1, "character",  "Project name of your res",
      "help",   "h", 0, "logical",  "This is Help!"),
    byrow=TRUE, ncol=5)


opt <- getopt(spec=spec)

if( !is.null(opt$help) || is.null(opt$genome) || is.null(opt$gtf) ){
    # ... 这里你也可以自定义一些东放在里面
    cat(paste(getopt(spec=spec, usage = T), "\n"))
    quit()
}



print(opt$genome)
print(opt$gtf)
print(opt$samples)
samples <- strsplit(opt$samples, ",")[[1]]
print(samples)
print(opt$output_path)
print(opt$project_name)
print(opt$samples_dir)
print(opt$lnc_gtf)
