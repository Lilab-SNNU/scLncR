# scLncR
A pipeline to predict and analysis lncRNA from scRNA-seq.

![alt text](<figures/scLncR_workflow.png>)

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
$ export PATH="/path/to/scLncR:$PATH"
```
The R package of hdWGCNA need install independent in terminal or Rstudio.
```R
# install Bioconductor
install.packages("BiocManager")
BiocManager::install()

# install hdWGCNA from GitHub
devtools::install_github('smorabit/hdWGCNA', ref='dev')
```

If you do not have conda, you can also manually install the required R packages listed below.

**R package list**
- seurat, version 4.3.0.1(https://satijalab.org/seurat/)
- seuratobject, version 4.1.3(https://satijalab.org/seurat/)
- monocle2, version 2.30.0(https://cole-trapnell-lab.github.io/monocle-release/docs/)
- hdWGCNA, version 0.4.4(https://smorabit.github.io/hdWGCNA/)
- stringr, version 1.5.1(https://cran.r-project.org/web/packages/stringr/index.html)
- shiny, version 1.8.1.1(https://shiny.posit.co/)
- shinyjs, version 2.1.0(https://cran.r-project.org/web/packages/shinyjs/index.html)
- tidyverse, version 2.0.0(https://www.tidyverse.org/)
- dplyr, version 1.1.4(https://dplyr.tidyverse.org/)
- psych, version 2.4.3(https://www.rdocumentation.org/packages/psych)
- pheatmap, version 1.0.12(https://cran.r-project.org/web/packages/pheatmap/index.html)
- ggsci, version 3.1.0(https://cran.r-project.org/web/packages/ggsci/index.html)
- ggplot2, version 3.5.1(https://ggplot2.tidyverse.org/)
- patchwork, version 1.2.0(https://patchwork.data-imaginist.com/)
- reshape2, version 1.4.4(https://cran.r-project.org/web/packages/reshape2/index.html)
- yaml, version 2.3.10(https://cran.r-project.org/web/packages/yaml/index.html)
- optparse, version 1.7.5(https://cran.r-project.org/web/packages/optparse/index.html)


**Bioinformatics Software**
- hisat2,  version  2.2.1(https://daehwankimlab.github.io/hisat2/)
- stringtie,  version  2.2.3(https://ccb.jhu.edu/software/stringtie/)
- gffcompare, version 0.12.6(https://ccb.jhu.edu/software/stringtie/gffcompare.shtml)
- gffread, version 0.9.12(https://ccb.jhu.edu/software/stringtie/gff.shtml#gffread)
- CPC2, version 1.0.1(https://github.com/gao-lab/CPC2_standalone/releases/tag/v1.0.1)
---

### Usage
scLncR supports two modes of usage: command-line operation and a Shiny GUI interface.

### Command-line usage:
***scLncR Main Program***
```shell
$ scLncR -h
scLncR v0.1.0 - Single-cell lncRNA Discovery Pipeline
Usage: scLncR <command> [options]

Available commands:
  prelnc        Predict and annotate lncRNAs from single-nucleus RNA-seq data 
  count         Get scRNA-seq expression count matrix 
  dataProcess   ScRNA-seq expression count preprocess and annotation 
  function      DownStream analysis to explore lncRNA function 
  shiny         Launch Shiny GUI

Examples:
  scLncR prelnc -c config.yaml
  scLncR prelnc --help

Options:
  -h, --help    Show this global help message

Note: Each command has its own --help. For example:
      scLncR prelnc --help
```
scLncR allows each subroutine to accept corresponding parameters via YAML configuration files. Users can modify the settings according to the parameter specifications provided in each YAML file and then execute the program sequentially using the following command. Alternatively, they may run specific modules independently based on their needs.

```shell
$ scLncR prelnc -c scLncR/R/confings/config_LncPre.yaml
$ scLncR count -c scLncR/R/confings/config_Count.yaml
$ scLncR dataProcess -c scLncR/R/confings/config_dataProcess.yaml
$ scLncR function -c scLncR/R/confings/config_function.yaml
```
### Run scLncR in graphical user interface(GUI) 
#### Method 1 :Open the shiny in Rstudio
```text
After opening RStudio, click the menu:
File → New Project → Existing Directory → Select the scLncR directory

In RStudio's file explorer:
1. Navigate to the scLncR/shiny_app directory
2. Double-click to open the app.R file

Method A: Using the Run Button
In the upper right corner of the app.R editor, click the "Run App" button.
Or click the drop-down arrow next to it and select the run port.

Method B: Using Keyboard Shortcuts
In the app.R file, press Ctrl + Shift + Enter (Windows/Linux)
Or Cmd + Shift + Enter (Mac)

Method C: Running in the Console
Ensure the working directory is correct:
setwd("/home/data/scLncR/scLncR/shiny_app")
shiny::runApp()
```
#### Method 2 :Open the shiny in command

```shell
$ scLCnR shiny
$ # Then Open the port print in shell,such as 
Setting up environment...
Found scLncR at: opt/scLncR/scLncR
Installation directory: opt/scLncR
Launching Shiny app from: opt/scLncR/shiny_app
Server will be available at: http://localhost:3838
Press Ctrl+C to stop the server


Listening on http://0.0.0.0:3838
```
Then, open the IP in your browser.

You can get the GUI:

![alt text](figures/Shiny_GUI.png)


---
## Contact us

If you encounter any problems while using scLncR, please send an email (glli@snnu.edu.cn) or submit the issues on GitHub (https://github.com/Lilab-SNNU/scLcnR/issues) and we will resolve it as soon as possible.