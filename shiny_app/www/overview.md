# scLncR Shiny GUI

## Overview

scLncR is a comprehensive pipeline for single-cell long non-coding RNA analysis. This Shiny GUI provides a user-friendly interface to configure and run the entire pipeline.

## Modules

### 1. LncRNA Prediction
- Predicts lncRNAs from single-cell/single-nucleus RNA-seq data
- Requires: Sample directories, reference genome, annotation GTF
- Outputs: Predicted lncRNA GTF file

### 2. Expression Counting
- Generates expression count matrices
- Uses CellRanger for alignment and counting
- Requires: lncRNA GTF from previous step

### 3. Data Processing
- Quality control, normalization, and preprocessing
- Cell type annotation using SingleR or scMM
- Clustering and visualization

### 4. Function Analysis
- Location analysis (differential expression)
- Trajectory analysis with Monocle2
- Co-expression network analysis with WGCNA

## Quick Start

1. Configure each module with your data paths
2. Validate the configurations
3. Generate configuration files
4. Run the pipeline
5. Monitor progress and view results

## File Requirements

### Input Files
- FASTQ/BAM files for each sample
- Reference genome (FASTA format)
- Gene annotation (GTF format)

### Output Structure