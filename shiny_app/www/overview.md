# scLncR Shiny GUI

## Overview

scLncR is a comprehensive pipeline for single-cell long non-coding RNA analysis. This Shiny GUI provides a user-friendly interface to configure and run the entire pipeline.

## Modules

### 1. LncRNA Prediction
- Raw FASTQ-first candidate lncRNA discovery (technology-aware)
- 10x/Smart-seq2/Drop-seq/generic interfaces with explicit input-read-role parsing
- Outputs: `final_lnc.gtf`, `final.lncRNA.fa`, `reference/combined_mRNA_lncRNA.gtf`, manifests, logs, and run report

### 2. Expression Counting
- lncRNA-aware quantification from raw FASTQ + augmented reference
- Primary supported path: 10x + Cell Ranger
- Experimental interfaces: Smart-seq2 and Drop-seq

### 3. Data Processing
- Quality control, normalization, and preprocessing
- Cell type annotation using SingleR or scMM
- Clustering and visualization

### 4. Function Analysis
- snRNA/scRNA expression enrichment analysis (differential expression between snRNA-seq and scRNA-seq groups)
- Trajectory analysis with Monocle2
- Co-expression network analysis with WGCNA

Note: enrichment results indicate relative expression patterns (snRNA-enriched or scRNA-enriched) and are not standalone evidence of true subcellular localization.

## Quick Start

1. Configure each module with your data paths
2. Validate the configurations
3. Generate configuration files
4. Run the pipeline
5. Monitor progress and view results

## File Requirements

### Input Files
- Raw FASTQ directory and optional sample metadata table
- Reference genome (FASTA format)
- Gene annotation (GTF format)

### Output Structure
