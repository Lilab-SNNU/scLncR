# scLncR
A pipeline to predict and analysis lncRNA from scRNA-seq.

---
### Download
```shell
$ git clone https://github.com/Lilab-SNNU/scLncR.git
```
---

### Install
#### Enveriment 

You can configure scLncR environment using conda.
```shell
$ cd scLncR
$ conda create -f scLncR.yaml 
$ conda activate scLncR
```

If you do not have conda, you can also manually install the required R packages listed below.

**R package list**
- seurat(v4.3.0.1): https://satijalab.org/seurat/
- seuratobject(v4.1.3): https://satijalab.org/seurat/
- monocle2(v2.30.0): https://cole-trapnell-lab.github.io/monocle-release/docs/
- stringr(v1.5.1): https://cran.r-project.org/web/packages/stringr/index.html
- shiny(version 1.8.1.1): https://shiny.posit.co/
- shinyjs(v2.1.0): https://cran.r-project.org/web/packages/shinyjs/index.html
- tidyverse(v2.0.0): https://www.tidyverse.org/
- dplyr(v1.1.4): https://dplyr.tidyverse.org/
- psych(v2.4.3): https://www.rdocumentation.org/packages/psych
- pheatmap(v1.0.12): https://cran.r-project.org/web/packages/pheatmap/index.html
- ggsci, version 3.1.0(https://cran.r-project.org/web/packages/ggsci/index.html)
- ggplot2, version 3.5.1(https://ggplot2.tidyverse.org/)
- patchwork, version 1.2.0(https://patchwork.data-imaginist.com/)
- reshape2, version 1.4.4(https://cran.r-project.org/web/packages/reshape2/index.html)


**Aligment Software**
- hisat2,  version  2.2.1(https://daehwankimlab.github.io/hisat2/)
- stringtie,  version  2.2.3(https://ccb.jhu.edu/software/stringtie/)
- gffcompare, version 0.12.6(https://ccb.jhu.edu/software/stringtie/gffcompare.shtml)
- gffread, version 0.9.12(https://ccb.jhu.edu/software/stringtie/gff.shtml#gffread)
- CPC2, version 1.0.1(https://github.com/gao-lab/CPC2_standalone/releases/tag/v1.0.1)
---

### Usage
#### - You can run all scLncR pipeline by one script
If you are deploying on a server, use SSH protocol for terminal deployment, please follow the steps below:
1. Fill the config file (https://github.com/Lilab-SNNU/scLncR/config.yaml), input full path of each required file.
2. Run scLncR
``` shell
$ Rscripts /opt/scLncR/scLncR_main.r -c /opt/scLncR/config.yaml
```

#### - Or you can run all scLncR pipeline by one script
If you are using it locally, and there is a GUI graphical interface, you can use shiny to run scLncR.

Run the Shiny.R file in terminal:
```shell
$ Rscripts /opt/scLncR/Shiny.R
```
  And, you will see the shiny GUI graphical interface in your browser. 

- The basic parameters are required, you need give the path of each required file.

![alt text](2025-05-13_17-01.png)

- The optional modules are optional, you can choose the modules you want to analysis in the pipeline and set the parameters.

![alt text](2025-05-13_17-12.png)

- Then, click the "Run" button to start the pipeline.

If the GUI interface does not automatically pop up, please check your browser settings or open the IP address which return by the terminal command in your browser.

---
### How to fill in the config.yaml file or Shiny Parameters? 
When opening the config file in text format, or use shiny there are some lines that need to be filled in, they are:
- ***Requireds***
  - **samples_dirs**: The directory where the scRNA-seq data is stored.
  - **project_name**: The name of the project, which will be used as the name of the output directory.
  - **output_path**: The path where the output files will be stored.
  - **genome_file**: You need to fill in the absolute path of the corresponding species genome here (**not the relative path!**)
  - **gtf_file**: You need to fill in the absolute path of the corresponding annotation file of species genome here (**not the relative path!**)
  - **lncrna_name**: The name of the lncRNA, which will be used in lncRNA prediction and select from scRNA-seq count matrix.
  - **threads**: The number of threads running, only the number in the input int format is supported here, for example: 1 or 2 or 3 or 4 or 5..., default is 4

- ***Optionals***

  - **Data process**: The data processing module
    - **marge**: Merge multiple samples or not, default is "Yes".
    - **threads**: The number of threads running, only the number in the input int format is supported here, for example: 1 or 2 or 3 or 4 or 5..., default is 4.
    - **nFeature_RNA**: The minimum number of features for each cell, only the number in the input int format is supported here, for example: 1000 or 2000 or 3000 or 4000 or 5000..., default is 2000.
    - **resolution**: Value of the resolution parameter, use a value above (below) 1.0 if you want to obtain a larger (smaller) number of communities. Default is 0.5.
    - **dims**: Dimensions of reduction to use as input.
  
  - **Cell Type Annotation**
    - **anno_method**: Method used in cell annotation, the methods "SingleR" and "scMayoMap" can be choosed 
    - **anno_ref**: The reference dataset used in cell annotation, the reference datasets, only used when annotation method "SingleR" is selected.
    - **anno_marker_file**: The marker gene files collected from the marker gene database, only used when annotation method "scMayoMap" is selected.
  
  - **Pserdotime Analysis**
    - **qval**: Gene filter threshold used in differential analysis process, default is 0.05.
    - **reduceDimensionMethod**: ReduceDimension method used in WGCNA analysis process, you can choose from "DDRTree", "ICA", "tSNE", "SimplePPT", "L1-graph", "SGL-tree", default is "DDRTree".

  - **WGCNA Analysis**
