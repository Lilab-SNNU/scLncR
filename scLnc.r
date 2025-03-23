library(getopt)

########################################***Para Design***########################################
## step1: genome, gtf, fastq, lnc_name, software_path; return: lnc_gtf                         ##
## step2: genome, gtf, lnc_gtf, fastq, sample_dir, pro_name; return: count matrix              ##
## step3: sample_dir, sample_name, hvg_genes, nFeature_RNA_min, nFeature_RNA_max,              ##
## step4: seu_obj                                                                              ##
########################################*****************########################################

## genome, gtf, samples, 

spec <- matrix(
    c("genome",  "G", 1, "character", "Genome file , format is fasta!",
      "gtf", "g", 1, "character",  "GTF annotation file, need to match your genome file version!",
      "sample_files",  "s", 2, "character",  "The single cell seq samples. separate each sample with a comma. \
                            eg:-s SRR111_1_S1_L001_R2_001.fastq.gz,SRR22_1_S1_L001_R2_001.fastq.gz",
      "project_name",  "p", 1, "character",  "Project name of your res",
      "hvg_genes", "hnum", 2, "integer", "The hvg numbers will be chosen in data process",
      "nFeature_RNA_min", "nRmin", 2,  "integer", "The min nFeature_RNA numbers will be chosen in data filter",
      "nFeature_RNA_max", "nRmax", 2,  "integer", "The max nFeature_RNA numbers will be chosen in data filter",
      "pfam_db",  "P", 2, "character",  "The dir of pfam database",
      "CPC2_path", "c", 2, "character", "The path of CPC2",
      "output_path",  "o", 2, "character",  "Output dir of res",
      "help",   "h", 0, "logical",  "This is Help!"),
    byrow=TRUE, ncol=5)


########################################***Part1 get lnc info***########################################

get_gtf <- function(genome_file, genome_gtf, sample_files, output_file="", merge_name="", CPC2_path="", pfam_db=""){
    path <- getwd()
    dir.create("./step1_res")
    print(path)
    system("mkdir hisat2 stringtie res")
    genome_index <- strsplit(genome_file, ".fa")[[1]][1]
    
    #===============***Step 1: Chain specific alignment and assembly of transcripts***===============#
    
    ## Hisat2 build index
    hisat2_build_command <- sprintf("hisat2-build %s %s", genome_file, genome_index)
    system(hisat2_build_command)
    samples <- c()
    
    ## Sample trans gtf get
    for(sample_file in sample_files){
        
        infos <- strsplit(sample_file, "/")[[1]]
        num <- length(infos)
        sample <- strsplit(infos[num], "_")[[1]][1]
        sample_dir <- paste(infos[1:num-1], collapse="/")
        samples <- append(samples, sample)
        
        ## Hisat2 aligment
        hisat2_aligment_command <- sprintf("hisat2 --dta -x %s -p 12 --rna-strandness RF -U %s -S %s/hisat2/lnc.%s.sam", genome_index, sample_file, path, sample) 
        system(hisat2_aligment_command)
        
        ## Samtools deal the result from hiast2
        samtools_command <- sprintf("samtools view -@ 4 -bS -q 30 ./hisat2/lnc.%s.sam | samtools sort -@ 4 - -o ./hisat2/lnc.%s.bam", sample, sample)
        system(samtools_command)
        system(sprintf("rm ./hisat2/lnc.%s.sam", sample))
        
        ## Stringtie assembly
        stringtie_command <- sprintf("stringtie -p 12 --rf -G %s -o ./stringtie/lnc.%s.gtf -l %s_lnc ./hisat2/lnc.%s.bam", genome_gtf, sample, sample, sample)
        system(stringtie_command)
    }

    #======================***Step 2: Filter lncRNA by the characteristics of lncRNA***======================#
    
    ## Merge all lnc_gtf
    system("ls ./stringtie/*.gtf > gtf_list.txt")
    if(merge_name==""){
        merge_name <- sample
    }else{
        merge_name <- merge_name
    }
    stringtie_merge <- sprintf("stringtie --merge -p 12 --rf -G %s -o ./stringtie/all.merged.gtf -l %s gtf_list.txt", genome_gtf, merge_name)
    system(stringtie_merge)
    
    ## Compare trans gtf with genome gtf
    gffcompare_command1 <- sprintf("gffcompare -r %s -o ./stringtie/gffcmp ./stringtie/all.merged.gtf", genome_gtf)
    system(gffcompare_command1)
    
    ## Filter the transcripts based on the results of gffcompare, and retain only transcripts with calss codes is "i, u, x, o"
    ## and lengths greater than 200bp according to the characteristics of lncRNA
    ## i: Completely contained within the reference intron
    ## u: Located in the intergenic region
    ## x: Overlap with reference exon but opposite direction to reference transcript
    ## o: Overlap with the reference exon and align with the direction of the reference transcript
    ## Then, extract the transcript sequence from the gtf file
    
    
    system("awk '$3==\"i\" || $3==\"u\" || $3==\"x\" || $3==\"o\" && $10>200 {print \"\\\"\" $5 \"\\\"\"}' ./stringtie/gffcmp.all.merged.gtf.tmap > ./stringtie/filter.iuxo.l200.transcript.id.txt")
    system("LC_ALL=C fgrep -f ./stringtie/filter.iuxo.l200.transcript.id.txt ./stringtie/gffcmp.annotated.gtf > ./stringtie/filter.iuxo.l200.gtf")
    system(sprintf("gffread ./stringtie/filter.iuxo.l200.gtf -g %s -w ./stringtie/filter.iuxo.l200.transcript.fa", genome_file))
    
    #===============***Step 3: Filter lncRNA filtering lncRNA through encoding function***===============#
    
    ## Filter out transcripts with protein encoding ability

    if(pfam_db==""){

        ## Filter out transcripts with protein encoding ability
        # CPC2_command <- sprintf("python %s/bin/CPC2.py -i ./stringtie/filter.iuxo.l200.transcript.fa -o ./stringtie/CPC2.out --ORF", CPC2_path)
        CPC2_command <- "CPC2.py -i ./stringtie/filter.iuxo.l200.transcript.fa -o ./stringtie/CPC2.out --ORF"
        system(CPC2_command)
        system("awk '$9==\"noncoding\"{print \"\\\"\" $1 \"\\\"\"}' ./stringtie/CPC2.out.txt > ./stringtie/tmp.CPC2.txt")
        system("cat ./stringtie/tmp.CPC2.txt | sort > ./stringtie/filter.iuxo.l200.ncd.transcript.txt")
        system(sprintf("LC_ALL=C fgrep -f ./stringtie/filter.iuxo.l200.ncd.transcript.txt ./stringtie/filter.iuxo.l200.gtf > %s", output_file))
        system(sprintf("gffread %s -g %s -w ./res/final.lncRNA.fa", output_file, genome_file))
    }else{
        ## Filter out transcripts with pfam database
        ## Translate the transcript into six possible protein sequences
        ## Set E-value<1e-5 as the threshold of pfam_scan
        CPC2_command <- sprintf("python %s/bin/CPC2.py -i ./stringtie/filter.iuxo.l200.transcript.fa -o ./stringtie/CPC2.out --ORF", CPC2_path)
        system(CPC2_command)
        system("awk '$9==\"noncoding\"{print \"\\\"\" $1 \"\\\"\"}' ./stringtie/CPC2.out.txt > ./stringtie/tmp.CPC2.txt")
        system("cat ./stringtie/tmp.CPC2.txt | sort > ./stringtie/filter.iuxo.l200.ncd.transcript.txt")
        system("LC_ALL=C fgrep -f ./stringtie/filter.iuxo.l200.ncd.transcript.txt ./stringtie/filter.iuxo.l200.gtf > ./stringtie/filter.iuxo.l200.ncd.gtf")
        system(sprintf("gffread ./stringtie/filter.iuxo.l200.ncd.gtf -g %s -w ./stringtie/filter.iuxo.l200.ncd.transcript.fa", genome_file))
        transeq_command <- "transeq ./stringtie/filter.iuxo.l200.ncd.transcript.fa ./stringtie/filter.iuxo.l200.ncd.pep6.fa -frame=6"
        system(transeq_command)
        pfam_scan_command <- sprintf("pfam_scan.pl -cpu 12 -fasta ./stringtie/filter.iuxo.l200.ncd.pep6.fa -dir %s -outfile ./stringtie/pfam.out", pfam_db)
        system(pfam_scan_command)
        system("sed \"/#/d\" ./stringtie/pfam.out | sed '/^$/d'  | awk '$13 < 1e-5 {print $1}' | sed \"s/_/\\t/g\" | awk '{print \"\\\"\" $1 \"\\\"\"}' > ./stringtie/pfam.coding.txt")
        system(sprintf("LC_ALL=C fgrep -v -f ./stringtie/pfam.coding.txt ./stringtie/filter.iuxo.l200.ncd.gtf > %s", output_file))
        system(sprintf("gffread %s -g %s -w ./res/final.lncRNA.fa", output_file, genome_file))
    }
    res <- c(samples, output_file, sample_dir)
    return (res)
}

#################################################***Part2 Single Cell Sequence Aligment***##################################################

## genome, gtf, lnc.gtf(Part 1), sample(Part 1), project_name = sample, sample_dir = fastq[:-1]

run_mkgtf <- function(lncgtf, genegtf, project_name, output_path="./"){

    merge_gtf <- sprintf("cat %s %s > tmp.gtf", lncgtf, genegtf)
    mkgtf_command <- sprintf("cellranger mkgtf tmp.gtf %s/%s.gtf", output_path, project_name)
    
    system(merge_gtf)
    system(mkgtf_command)

    gtf_file <- sprintf("%s/%s.gtf", output_path, project_name)
    return(gtf_file)
}

run_mkref <- function(genome_file, gtf_file, project_name, output_path){

    mkref_command <- sprintf("cellranger mkref --fasta %s --genes %s --genome %s --output-dir %s/%s --nthreads 12", 
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
        cellranger_count_command <- sprintf("cellranger count --id=%s --fastqs=%s --sample=%s --transcriptome=%s --localcores=%s %s/%s", 
                                            sample, samples_dir, sample, trans_name, cores_num, output_path, sample)
        print(cellranger_count_command)
        system(cellranger_count_command)
    }
    return (output_path)
}

opt <- getopt(spec=spec)

if( !is.null(opt$help) || is.null(opt$genome) || is.null(opt$gtf) ){
    cat(paste(getopt(spec=spec, usage = T), "\n"))
    quit()
}


# lncgtf <- opt$lnc_gtf
# genegtf <- opt$gtf
# project_name <- opt$project_name
# output_path <- opt$output_path
# genome_file <- opt$genome
# samples <- strsplit(opt$samples, ",")[[1]]
# samples_dir <- opt$samples_dir

if(is.null(opt$mean)){ opt$pfam_db =""}
res_info <- get_gtf(genome_file=opt$genome,
                    genome_gtf=opt$gtf,
                    sample_files=strsplit(opt$sample_files, ",")[[1]],
                    merge_name=opt$project_name,
                    pfam_db=opt$pfam_db,
                    CPC2_path=opt$CPC2_path,
                    output_file="./lnc_gtf")
n <- length(res_info)
m <- n-2
samples <- res_info[1:m] 
lnc_gtf <-  res_info[n-1]
samples_dir <- res_info[n]

project_name <- opt$project_name
output_path <- opt$output_path
gtf_file <- run_mkgtf(lncgtf=lnc_gtf, 
                      genegtf=opt$gtf, 
                      project_name=project_name, 
                      output_path=output_path)
print("####################***Function mkgtf done***####################")

trans_name <- run_mkref(genome_file=opt$genome,
                        gtf_file=gtf_file,
                        project_name=project_name,
                        output_path=output_path)
print("####################***Function mkref done***####################")

run_cellranger(samples=samples, 
               samples_dir=samples_dir, 
               trans_name=trans_name, 
               cores_num=24, 
               output_path=output_path)
print("####################***Function count done***####################")
