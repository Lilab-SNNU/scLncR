library(getopt)

spec <- matrix(
    c("genome",  "G", 2, "character", "Genome file , format is fasta!",
      "gtf", "g", 2, "character",  "GTF annotation file, need to match your genome file version!",
      "samples",  "s", 2, "character",  "The single cell seq samples. separate each sample with a comma. \
                            eg:-s SRR111_1_S1_L001_R2_001.fastq.gz,SRR22_1_S1_L001_R2_001.fastq.gz",
      "output_file",  "o", 2, "character",  "Output file path",
      "lncrna_name",  "n", 1, "character",  "If you have multiple samples, you can set the name prefix for your lncrnas, which will become the prefix to all lncrna in result, \
                            eg: -n test, the 'transcript_id' in final gtf result will be 'test.1.1', the transcript name in final fasta result will be >test.1.1. \
                            If you do not set the paramater, it will be choose the final sample name as the lncrna name prefix.",
      "CPC2_path",  "c", 2, "character",  "The path of software CPC2",
      "pfam_db",  "p", 2, "character",  "The dir of pfam database",
      "help",   "h", 0, "logical",  "This is Help!"),
    byrow=TRUE, ncol=5)


get_gtf <- function(genome_file, genome_gtf, sample_files, output_file, merge_name="", CPC2_path="", pfam_db=""){
    path <- getwd()
    print(path)
    system("mkdir hisat2 stringtie res")
    genome_index <- strsplit(genome_file, ".fa")[[1]][1]
    
    ###############***Step 1: Chain specific alignment and assembly of transcripts***###############
    
    ## Hisat2 build index
    hisat2_build_command <- sprintf("hisat2-build %s %s", genome_file, genome_index)
    system(hisat2_build_command)
    
    ## Sample trans gtf get
    for(sample_file in sample_files){
        
        infos <- strsplit(sample_file, "/")[[1]]
        num <- length(infos)
        sample <- strsplit(infos[num], "_")[[1]][1]
        
        ## Hisat2 aligment
        hisat2_aligment_command <- sprintf("hisat2 --dta -x %s -p 4 --rna-strandness RF -U %s -S %s/hisat2/lnc.%s.sam", genome_index, sample_file, path, sample) 
        system(hisat2_aligment_command)
        
        ## Samtools deal the result from hiast2
        samtools_command <- sprintf("samtools view -@ 4 -bS -q 30 ./hisat2/lnc.%s.sam | samtools sort -@ 4 - -o ./hisat2/lnc.%s.bam", sample, sample)
        system(samtools_command)
        system(sprintf("rm ./hisat2/lnc.%s.sam", sample))
        
        ## Stringtie assembly
        stringtie_command <- sprintf("stringtie -p 4 --rf -G %s -o ./stringtie/lnc.%s.gtf -l %s_lnc ./hisat2/lnc.%s.bam", genome_gtf, sample, sample, sample)
        system(stringtie_command)
    }
    
    ###############***Step 2: Filter lncRNA by the characteristics of lncRNA***###############
    
    ## Merge all lnc_gtf
    system("ls ./stringtie/*.gtf > gtf_list.txt")
    if(merge_name==""){
        merge_name <- sample
    }else{
        merge_name <- merge_name
    }
    stringtie_merge <- sprintf("stringtie --merge -p 4 --rf -G %s -o ./stringtie/all.merged.gtf -l %s gtf_list.txt", genome_gtf, merge_name)
    system(stringtie_merge)
    
    ## Compare trans gtf with genome gtf
    gffcompare_command1 <- sprintf("gffcompare -r %s -o ./stringtie/gffcmp ./stringtie/all.merged.gtf", genome_gtf)
    system(gffcompare_command1)
    
    ## Filter the transcripts based on the results of gffcompare, and retain only transcripts with calss codes is "i, u, x, o"
    ## and lengths greater than 200bp according to the characteristics of lncRNA
    # i: Completely contained within the reference intron
    # u: Located in the intergenic region
    # x: Overlap with reference exon but opposite direction to reference transcript
    # o: Overlap with the reference exon and align with the direction of the reference transcript
    ## Then, extract the transcript sequence from the gtf file
    system("awk '$3==\"i\" || $3==\"u\" || $3==\"x\" || $3==\"o\" && $10>200 {print \"\\\"\" $5 \"\\\"\"}' ./stringtie/gffcmp.all.merged.gtf.tmap > ./stringtie/filter.iuxo.l200.transcript.id.txt")
    system("LC_ALL=C fgrep -f ./stringtie/filter.iuxo.l200.transcript.id.txt ./stringtie/gffcmp.annotated.gtf > ./stringtie/filter.iuxo.l200.gtf")
    system(sprintf("gffread ./stringtie/filter.iuxo.l200.gtf -g %s -w ./stringtie/filter.iuxo.l200.transcript.fa", genome_file))
    
    ###############***Step 3: Filter lncRNA filtering lncRNA through encoding function***###############
    
    ## Filter out transcripts with protein encoding ability
    CPC2_command <- sprintf("python %s/bin/CPC2.py -i ./stringtie/filter.iuxo.l200.transcript.fa -o ./stringtie/CPC2.out --ORF", CPC2_path)
    system(CPC2_command)
    system("awk '$9==\"noncoding\"{print \"\\\"\" $1 \"\\\"\"}' ./stringtie/CPC2.out.txt > ./stringtie/tmp.CPC2.txt")
    system("cat ./stringtie/tmp.CPC2.txt | sort > ./stringtie/filter.iuxo.l200.ncd.transcript.txt")
    system("LC_ALL=C fgrep -f ./stringtie/filter.iuxo.l200.ncd.transcript.txt ./stringtie/filter.iuxo.l200.gtf > ./stringtie/filter.iuxo.l200.ncd.gtf")
    system(sprintf("gffread ./stringtie/filter.iuxo.l200.ncd.gtf -g %s -w ./stringtie/filter.iuxo.l200.ncd.transcript.fa", genome_file))
    
    ## Filter out transcripts with pfam database
    ## Translate the transcript into six possible protein sequences
    transeq_command <- "transeq ./stringtie/filter.iuxo.l200.ncd.transcript.fa ./stringtie/filter.iuxo.l200.ncd.pep6.fa -frame=6"
    system(transeq_command)
    
    ## Set E-value<1e-5 as the threshold of pfam_scan
    pfam_scan_command <- sprintf("pfam_scan.pl -cpu 4 -fasta ./stringtie/filter.iuxo.l200.ncd.pep6.fa -dir %s -outfile ./stringtie/pfam.out", pfam_db)
    system(pfam_scan_command)
    system("sed \"/#/d\" ./stringtie/pfam.out | sed '/^$/d'  | awk '$13 < 1e-5 {print $1}' | sed \"s/_/\\t/g\" | awk '{print \"\\\"\" $1 \"\\\"\"}' > ./stringtie/pfam.coding.txt")
    system(sprintf("LC_ALL=C fgrep -v -f ./stringtie/pfam.coding.txt ./stringtie/filter.iuxo.l200.ncd.gtf > %s", output_file))
    system(sprintf("gffread %s -g %s -w ./res/final.lncRNA.fa", output_file, genome_file))
}

# setwd("~/work/scRNA/step1_lnc_pre/")
# genome_file <- "~/work/scRNA/genome/TAIR10.fa"
# genome_gtf <- "~/work/scRNA/genome/Araport11_current.gtf"
# samples <- c("~/work/scRNA/data/SRR10620013_S1_L001_R2_001.fastq.gz", "~/work/scRNA/data/SRR10620014_S1_L001_R2_001.fastq.gz")
# output_file <- "~/work/scRNA/step1_lnc_pre/res/final.lncRNA.gtf"
# CPC2_path <- "/home/li/anwser/ysw/software/CPC2_standalone-1.0.1"
# pfam_db <- "~/work/scRNA/step1_lnc_pre/pfam"

opt <- getopt(spec=spec)

if( !is.null(opt$help) || is.null(opt$genome) || is.null(opt$gtf) ){
    cat(paste(getopt(spec=spec, usage = T), "\n"))
    quit()
}

genome_file <- opt$genome
genome_gtf <- opt$gtf
samples <- strsplit(opt$samples, ",")[[1]]
output_file <- opt$output_file
lncrna_name <- opt$lncrna_name
CPC2_path <- opt$CPC2_path
pfam_db <- opt$pfam_db

get_gtf(genome_file=genome_file, 
        genome_gtf=genome_gtf, 
        sample_files=samples, 
        output_file=output_file, 
        merge_name=lncrna_name, 
        CPC2_path=CPC2_path, 
        pfam_db=pfam_db)

