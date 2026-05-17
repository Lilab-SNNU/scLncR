# scLncR

Primary project documentation is maintained in [ReadMe.md](./ReadMe.md).

Current design summary:
- raw FASTQ-first, technology-aware workflow;
- optional raw FASTQ QC with FastQC and MultiQC before prelnc/count;
- prelnc performs candidate lncRNA discovery from raw FASTQ (with internal BAM intermediates);
- count performs lncRNA-aware quantification from raw FASTQ + augmented reference.

## Raw FASTQ quality control

Run optional QC before prelnc/count:

```bash
scLncR qc -c R/confings/config_QC.yaml
```

Equivalent shell wrapper:

```bash
bash scripts/run_fastq_qc.sh -c R/confings/config_QC.yaml
```

This module only runs FastQC and MultiQC. It does not trim, filter, remove adapters, or modify raw FASTQ files. Outputs are written under `output_dir/fastqc`, `output_dir/multiqc`, `output_dir/logs`, and `output_dir/qc_summary.md`.

## Smart-seq2 Count Workflow

Smart-seq2 is a full-length-like single-cell RNA-seq technology and should not be quantified with Cell Ranger. For `sequencing_platform: "smartseq2"` and `count_engine: "featurecounts"`, scLncR runs:

```text
Raw Smart-seq2 FASTQ
→ HISAT2 alignment
→ samtools sort/index
→ featureCounts with combined_mRNA_lncRNA.gtf
→ gene x sample raw count matrix
```

Example:

```bash
scLncR count -c R/confings/config_Count.yaml
```

Minimal Smart-seq2 configuration:

```yaml
sequencing_platform: "smartseq2"
count_engine: "featurecounts"
samples_dirs: "/path/to/smartseq2_fastq"
combined_gtf: "/path/to/step_LncPre/reference/combined_mRNA_lncRNA.gtf"
read_layout: "paired"
strandness: "unstranded"
```

The main outputs are `featurecounts/smartseq2_featureCounts.txt`, `featurecounts/smartseq2_count_matrix.tsv`, and `featurecounts/smartseq2_count_matrix.rds`. The count matrix rows are `gene_id` values from the augmented GTF and columns are sample IDs.

## Smart-seq2 DataProcess And Benchmark

Smart-seq2 count output can be passed to `dataProcess` through the featureCounts matrix input mode:

```yaml
input_format: "featurecounts_matrix"
sequencing_platform: "smartseq2"
counts_matrix: "/path/to/featurecounts/smartseq2_count_matrix.tsv"
samples_info: "/path/to/samples_info.txt"
anno_method: "none"
lnc_name: "AthLnc"
```

Run:

```bash
scLncR dataProcess -c R/confings/examples/config_dataProcess.smartseq2.yaml
Rscript scripts/run_normalization_benchmark.R -c R/confings/examples/config_normalization_benchmark.smartseq2.yaml
```

The Smart-seq2 example configs are stored in `R/confings/examples/`. This path preserves the existing 10x `filtered_feature_bc_matrix` input mode while adding a separate featureCounts matrix input path for Smart-seq2.
