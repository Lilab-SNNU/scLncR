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

opt <- getopt(spec=spec)

if( !is.null(opt$help) || is.null(opt$genome) || is.null(opt$gtf) ){
    cat(paste(getopt(spec=spec, usage = T), "\n"))
    quit()
}


#lncgtf <- "~/work/scRNA/step1_lnc_pre/res/final.lncRNA.gtf"
#genegtf <- "~/work/scRNA/genome/Araport11_current.gtf" 
#project_name <- "test"  
#output_path <- "/home/li/anwser/ysw/work/scRNA/step2_cellranger"
#genome_file <- "~/work/scRNA/genome/TAIR10.fa"
#samples <- c("SRR10620013", "SRR10620014")
#samples_dir <- "/home/li/anwser/ysw/work/scRNA/data/test/"

print(opt)

lncgtf <- opt$lnc_gtf
genegtf <- opt$gtf
project_name <- opt$project_name
output_path <- opt$output_path
genome_file <- opt$genome
samples <- strsplit(opt$samples, ",")[[1]]
samples_dir <- opt$samples_dir

gtf_file <- run_mkgtf(lncgtf=lncgtf, 
                      genegtf=genegtf, 
                      project_name=project_name, 
                      output_path=output_path)
print("####################***Function mkgtf done***####################")

trans_name <- run_mkref(genome_file=genome_file,
                        gtf_file=gtf_file,
                        project_name=project_name,
                        output_path=output_path)
print("####################***Function mkref done***####################")

run_cellranger(samples=samples, 
               samples_dir=samples_dir, 
               trans_name=trans_name, 
               cores_num=4, 
               output_path=output_path)
print("####################***Function count done***####################")

