# scLncR Pipeline Overview

The scLncR pipeline is designed for comprehensive analysis of long non-coding RNAs (lncRNAs) from single-cell RNA sequencing data. The pipeline consists of four main modules:

## Pipeline Modules

1. **LncRNA Prediction** - Identify and annotate lncRNAs from single-nucleus RNA-seq data
2. **Expression Counting** - Generate expression count matrices from aligned BAM files
3. **Data Processing** - Quality control, normalization, and preprocessing of expression data
4. **Functional Analysis** - Downstream analysis to explore lncRNA function

## Workflow

The standard workflow proceeds sequentially through these modules, but each module can also be run independently if you have intermediate data.

## Features

- **Complete Pipeline Management**: Configure and run all four scLncR modules
- **Visual Parameter Configuration**: Intuitive interface for all YAML parameters
- **File Browser Integration**: Built-in file/directory selection
- **Real-time Monitoring**: Live job monitoring and log display
- **Configuration Export**: Save and load configuration files
- **Project Management**: Organize analyses by project

## Quick Start

1. **Install R dependencies**:
   ```bash
   Rscript -e "install.packages(c('shiny', 'shinydashboard', 'shinyjs', 'shinyWidgets', 'shinyFiles', 'yaml', 'DT', 'fs', 'processx'))"