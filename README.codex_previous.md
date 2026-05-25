# scLncR

scLncR is a raw FASTQ-first, technology-aware workflow for single-cell lncRNA candidate discovery and downstream lncRNA analysis.

The tested workflows currently include:

- 10x-style scRNA/snRNA-seq FASTQ processing.
- Smart-seq2 FASTQ processing with featureCounts quantification.
- Optional raw FASTQ QC, candidate lncRNA prediction, lncRNA-aware counting, Seurat preprocessing, normalization benchmarking, and downstream function analysis.

![scLncR workflow](figures/scLncR_workflow.png)

## 1. Installation

### 1.1 Clone scLncR

```bash
git clone https://github.com/Lilab-SNNU/scLncR.git
cd scLncR
```

### 1.2 Create the conda environment

The recommended installation file is:

```bash
conda env create -f scLncR.yaml
conda activate scLncR
```

For reproducibility, exported environment files are also provided:

```bash
# Minimal manually requested packages
conda env create -f envs/scLncR_environment_from_history.yaml

# Full exported environment, more reproducible but less portable
conda env create -f envs/scLncR_environment_full.yaml

# Same-platform exact reproduction
conda create --name scLncR --file envs/scLncR_conda_explicit_linux64.txt
```

If the environment already exists and only GUI packages are missing, install them in the active environment:

```bash
conda activate scLncR
conda install -c conda-forge r-dt r-shiny r-shinyjs r-shinydashboard r-shinywidgets r-shinyfiles r-fs
```

### 1.3 Make the CLI available

The `scLncR` launcher sets `SC_LNCR_HOME` automatically before calling `R/main.R`. Users normally only need to make the launcher executable and add the repository directory to `PATH`.

```bash
chmod +x scLncR
echo "export PATH=\"$(pwd):\$PATH\"" >> ~/.bashrc
source ~/.bashrc
conda activate scLncR
scLncR -h
```

If you run `R/main.R` directly instead of using the launcher, set `SC_LNCR_HOME` manually:

```bash
export SC_LNCR_HOME=/path/to/scLncR
Rscript "$SC_LNCR_HOME/R/main.R" --help
```

## 2. Dependencies By Module

Most users should install the conda environment above. The list below explains which tools are required by each workflow step.

| Module | Main dependencies |
| --- | --- |
| QC | FastQC, MultiQC |
| prelnc | HISAT2, samtools, StringTie, gffcompare, gffread, CPC2 |
| 10x count | Cell Ranger |
| Smart-seq2 count | HISAT2, samtools, featureCounts/subread |
| dataProcess | Seurat, dplyr/tidyverse, ggplot2, SingleR optional, scMayoMap optional |
| normalization benchmark | Seurat, ggplot2, dplyr, reshape2, yaml; scran optional |
| function | Seurat, monocle, VGAM, hdWGCNA, WGCNA, UCell, Biobase |
| Shiny GUI | shiny, shinydashboard, shinyjs, shinyWidgets, shinyFiles, DT, fs |

Cell Ranger and CPC2 may require separate installation depending on your system and license/distribution policy. Make sure each executable is available in `PATH` before running the corresponding module.

```bash
which fastqc
which multiqc
which hisat2
which stringtie
which gffcompare
which gffread
which CPC2.py
which cellranger
which featureCounts
```

## 3. Workflow Overview

```text
Raw FASTQ
-> optional FASTQ QC
-> technology-aware lncRNA candidate discovery
-> final_lnc.gtf and combined_mRNA_lncRNA.gtf
-> technology-aware quantification
-> mRNA + lncRNA count matrix
-> dataProcess
-> normalization benchmark
-> downstream function analysis
```

The main design is raw FASTQ-first. BAM files generated during prelnc or count are internal intermediates, not the recommended user-facing input.

## 4. Quick Start: 10x Workflow

Edit the default configs in `R/confings/` first, especially paths to FASTQ, genome FASTA, GTF, output directories, and lncRNA prefix.

```bash
conda activate scLncR

# Step 0: raw FASTQ QC, optional but recommended
scLncR qc -c R/confings/config_QC.yaml

# Step 1: candidate lncRNA discovery
scLncR prelnc -c R/confings/config_LncPre.yaml

# Step 2: 10x lncRNA-aware quantification
scLncR count -c R/confings/config_Count.yaml

# Step 3: preprocessing and optional annotation
scLncR dataProcess -c R/confings/config_dataProcess.yaml

# Step 4: normalization benchmark
Rscript scripts/run_normalization_benchmark.R \
  -c R/confings/config_normalization_benchmark.yaml

# Step 5: downstream function analysis
scLncR function -c R/confings/config_function.yaml
```

Important 10x notes:

- `I1` is the sample index read.
- `R1` contains cell barcode and UMI information.
- `R2` is the cDNA/insert read.
- prelnc uses R2 reads as candidate transcript evidence.
- 10x lncRNA prediction should be interpreted as candidate evidence, not full-length transcript reconstruction.
- count uses raw FASTQ again with Cell Ranger and the augmented reference.

## 5. Quick Start: Smart-seq2 Workflow

Smart-seq2 uses example configs under `R/confings/examples/`.

```bash
conda activate scLncR

scLncR qc -c R/confings/examples/config_QC.smartseq2.yaml
scLncR prelnc -c R/confings/examples/config_LncPre.smartseq2.yaml
scLncR count -c R/confings/examples/config_Count.smartseq2.yaml
scLncR dataProcess -c R/confings/examples/config_dataProcess.smartseq2.yaml

Rscript scripts/run_normalization_benchmark.R \
  -c R/confings/examples/config_normalization_benchmark.smartseq2.yaml
```

Smart-seq2 notes:

- Smart-seq2 is full-length-like and should not use Cell Ranger.
- scLncR uses HISAT2, samtools, and featureCounts for Smart-seq2 gene-level quantification.
- The main count matrix is `featurecounts/smartseq2_count_matrix.tsv`.
- `dataProcess` reads Smart-seq2 counts with `input_format: "featurecounts_matrix"`.

## 6. Configuration Files

Main configs:

| Config | Purpose |
| --- | --- |
| `R/confings/config_QC.yaml` | FastQC/MultiQC raw FASTQ QC |
| `R/confings/config_LncPre.yaml` | candidate lncRNA discovery and augmented GTF construction |
| `R/confings/config_Count.yaml` | 10x Cell Ranger or Smart-seq2 featureCounts quantification |
| `R/confings/config_dataProcess.yaml` | Seurat preprocessing and optional annotation |
| `R/confings/config_normalization_benchmark.yaml` | normalization stability benchmark |
| `R/confings/config_function.yaml` | snRNA/scRNA enrichment, Monocle2, WGCNA |

Smart-seq2 examples:

```text
R/confings/examples/config_QC.smartseq2.yaml
R/confings/examples/config_LncPre.smartseq2.yaml
R/confings/examples/config_Count.smartseq2.yaml
R/confings/examples/config_dataProcess.smartseq2.yaml
R/confings/examples/config_normalization_benchmark.smartseq2.yaml
```

Only a few fields usually need editing before a run:

- FASTQ directory: `input_dir` or `samples_dirs`
- sample metadata: `samples_info`
- genome FASTA: `genome_file`
- reference GTF: `gtf_file` or `known_gtf`
- lncRNA prefix: `lncrna_name`, `lnc_name`, or `lncRNA_name`
- output directory: `output_dir` or `output_path`
- platform settings: `sequencing_platform`, `read_layout`, `count_engine`

## 7. Main Outputs

| Step | Main output |
| --- | --- |
| QC | `qc_summary.md`, FastQC reports, MultiQC report |
| prelnc | `final_lnc.gtf`, `final.lncRNA.fa`, `reference/combined_mRNA_lncRNA.gtf` |
| count 10x | Cell Ranger `filtered_feature_bc_matrix` |
| count Smart-seq2 | `featurecounts/smartseq2_count_matrix.tsv` |
| dataProcess | `data_preprocess/preprocessed_result.rds`, `data_annotation/anno_result.rds` |
| normalization benchmark | `normalization_benchmark_report.md`, marker/stability tables and figures |
| function | snRNA/scRNA enrichment, Monocle2 trajectory, WGCNA outputs |

## 8. Shiny GUI

The GUI is useful for configuring and launching workflows interactively.

```bash
conda activate scLncR
scLncR shiny
```

By default, the GUI is served at:

```text
http://localhost:3838
```

On a server, open the displayed host/port through your browser or SSH tunnel according to your server policy.

You can also launch it from RStudio:

```r
setwd("/path/to/scLncR/shiny_app")
shiny::runApp()
```

## 9. Normalization Benchmark

scLncR includes a reviewer-oriented benchmark to compare normalization strategies and assess whether separate lncRNA/mRNA normalization changes lncRNA-associated signals.

```bash
Rscript scripts/run_normalization_benchmark.R \
  -c R/confings/config_normalization_benchmark.yaml
```

Use a Seurat object from `dataProcess` as input:

```text
data_preprocess/preprocessed_result.rds
data_annotation/anno_result.rds
```

The benchmark evaluates marker counts, lncRNA marker overlap, lncRNA signal distribution, HVG overlap, and group-level expression correlations. It is intended to assess stability, not to prove that one normalization strategy is universally superior.

## 10. Downstream Function Analysis

```bash
scLncR function -c R/confings/config_function.yaml
```

Available module keys:

- `location`: snRNA/scRNA expression enrichment analysis.
- `monocle2`: trajectory analysis.
- `wgcna`: co-expression network analysis.

The `location` key is kept for backward compatibility. Its output should be interpreted as snRNA-enriched or scRNA-enriched expression, not as direct nuclear/cytoplasmic localization evidence.

## 11. Troubleshooting

Check basic CLI availability:

```bash
conda activate scLncR
which scLncR
scLncR -h
```

If `scLncR` is not found, add the repository directory to `PATH`:

```bash
cd /path/to/scLncR
chmod +x scLncR
echo "export PATH=\"$(pwd):\$PATH\"" >> ~/.bashrc
source ~/.bashrc
```

Check module help:

```bash
scLncR qc --help
scLncR prelnc --help
scLncR count --help
scLncR dataProcess --help
scLncR function --help
scLncR shiny --help
```

Common issues:

- Missing external tools: install the tool and ensure it is in `PATH`.
- Missing Cell Ranger: only affects 10x count.
- Missing CPC2: affects prelnc coding-potential filtering.
- Missing `DT` or Shiny packages: affects GUI only.
- Invalid FASTQ naming: check `fastq_sample_regex`, `r1_pattern`, `r2_pattern`, and `i1_pattern`.
- Missing Cell Ranger matrix: check the count output path and `filtered_feature_bc_matrix`.
- Wrong benchmark input: use a Seurat RDS from `dataProcess`.

## 12. Limitations

- 10x lncRNA prediction provides candidate transcript evidence, not definitive full-length reconstruction.
- lncRNAs are often lowly expressed and sparse in single-cell data.
- snRNA/scRNA enrichment results should not be used alone as proof of nuclear/cytoplasmic localization.
- Drop-seq and generic modes are experimental/planned interfaces in this version.
- Biological conclusions should be treated as hypothesis-generating unless supported by independent evidence.

## 13. Contact

For questions or issues, contact the maintainers or open an issue at:

```text
https://github.com/Lilab-SNNU/scLncR/issues
```
