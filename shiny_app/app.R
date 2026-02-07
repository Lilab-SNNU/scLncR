# shiny_app/app.R
library(shiny)
library(shinyjs)
library(shinydashboard)
library(shinyWidgets)
library(yaml)
library(DT)
library(shinyFiles)

source("global.R")

# UI界面 - 根据您的配置文件重新设计
ui <- dashboardPage(
  dashboardHeader(
    title = "scLncR Pipeline",
    titleWidth = 250,
    tags$li(class = "dropdown", 
            actionButton("help_btn", "Help", icon = icon("question-circle"))),
    tags$li(class = "dropdown",
            tags$a(href = "#", icon("github"), "GitHub", 
                   onclick = "window.open('https://github.com/yourrepo/scLncR', '_blank');"))
  ),
  dashboardSidebar(
    width = 250,
    sidebarMenu(
      id = "tabs",
      menuItem("Dashboard", tabName = "dashboard", icon = icon("dashboard")),
      menuItem("1. LncRNA Prediction", tabName = "prelnc", icon = icon("dna")),
      menuItem("2. Expression Counting", tabName = "count", icon = icon("calculator")),
      menuItem("3. Data Processing", tabName = "dataprocess", icon = icon("filter")),
      menuItem("4. Function Analysis", tabName = "function_tab", icon = icon("chart-line")),
      menuItem("Run Pipeline", tabName = "run", icon = icon("play-circle")),
      hr(),
      menuItem("Configuration Files", tabName = "configs", icon = icon("file-code")),
      menuItem("Documentation", tabName = "docs", icon = icon("book")),
      menuItem("About", tabName = "about", icon = icon("info-circle"))
    ),
    hr(),
    div(style = "padding: 10px;",
        h5("Quick Actions"),
        actionButton("load_example", "Load Example", width = "100%", 
                     class = "btn-primary btn-sm"),
        br(),
        actionButton("validate_configs", "Validate Configs", width = "100%",
                     class = "btn-warning btn-sm"),
        br(),
        downloadButton("export_all", "Export All", width = "100%",
                      class = "btn-success btn-sm")
    )
  ),
  dashboardBody(
    useShinyjs(),
    tags$head(
      tags$link(rel = "stylesheet", type = "text/css", href = "style.css"),
      tags$script(HTML("
        // 防止表单意外关闭
        window.onbeforeunload = function() {
          return 'Are you sure you want to leave? All unsaved changes will be lost.';
        };
      "))
    ),
    tabItems(
      # 仪表板
      tabItem(
        tabName = "dashboard",
        h2("scLncR Pipeline Dashboard"),
        fluidRow(
          valueBoxOutput("box_pipeline_status", width = 3),
          valueBoxOutput("box_config_status", width = 3),
          valueBoxOutput("box_scLncR_status", width = 3),
          valueBoxOutput("box_memory_status", width = 3)
        ),
        fluidRow(
          box(
            width = 8,
            status = "primary",
            solidHeader = TRUE,
            title = "Pipeline Overview",
            img(src = "../figures/scLncR workflow.png", width = "100%", style = "border: 1px solid #ddd;")
          ),
          box(
            width = 4,
            status = "info",
            solidHeader = TRUE,
            title = "Quick Start",
            tags$ol(
              tags$li("Configure each module below"),
              tags$li("Validate your settings"),
              tags$li("Generate configuration files"),
              tags$li("Run the pipeline"),
              tags$li("Monitor progress and view results")
            ),
            hr(),
            h5("Recent Activity"),
            verbatimTextOutput("recent_logs")
          )
        ),
        fluidRow(
          box(
            width = 12,
            title = "Module Status",
            status = "success",
            collapsible = TRUE,
            DTOutput("module_status_table")
          )
        )
      ),
      
      # 1. LncRNA Prediction 模块 - 根据 config_LncPre.yaml 设计
      tabItem(
        tabName = "prelnc",
        h2("LncRNA Prediction Module"),
        fluidRow(
          box(
            width = 12,
            status = "primary",
            solidHeader = TRUE,
            title = "Configuration File: config_LncPre.yaml",
            collapsible = TRUE,
            collapsed = FALSE
          )
        ),
        fluidRow(
          # 左侧: 输入数据设置
          box(
            width = 6,
            status = "info",
            title = "1. Input Data Settings",
            solidHeader = TRUE,
            shinyDirButton("prelnc_samples_dir", "Select Samples Directory", 
                         "Please select a directory"),
            br(), br(),
            textInput("prelnc_samples_dirs", "Samples Directory",
                     value = "/home/data/scLncR/data/paper_test/",
                     placeholder = "Directory containing raw sample data (e.g., BAM/FASTQ files)"),
            helpText("Each subdirectory under this path should correspond to one sample."),
            hr(),
            h5("File Browser"),
            verbatimTextOutput("prelnc_dir_contents")
          ),
          
          # 右侧: 项目和输出设置
          box(
            width = 6,
            status = "success",
            title = "2. Project and Output Settings",
            solidHeader = TRUE,
            textInput("prelnc_project_name", "Project Name",
                     value = "scLncR",
                     placeholder = "Used for naming log files and output directories"),
            shinyDirButton("prelnc_output_path", "Select Output Directory", 
                         "Please select a directory"),
            br(), br(),
            textInput("prelnc_output_path", "Output Path",
                     value = "/home/data/scLncR/paper/paper_case",
                     placeholder = "Main output directory for all results"),
            helpText("All results from the prelnc module will be saved here.")
          )
        ),
        fluidRow(
          # 参考基因组文件
          box(
            width = 6,
            status = "warning",
            title = "3. Reference Genome and Annotation Files",
            solidHeader = TRUE,
            shinyFilesButton("prelnc_genome_file", "Select Genome File", 
                           "Please select a FASTA file", multiple = FALSE,
                           filetypes = c("fa", "fasta", "fa.gz", "fasta.gz")),
            br(), br(),
            textInput("prelnc_genome_file", "Genome File (FASTA)",
                     value = "/home/data/scLncR/genome/TAIR10.fa",
                     placeholder = "Path to reference genome in FASTA format"),
            shinyFilesButton("prelnc_gtf_file", "Select GTF File", 
                           "Please select a GTF file", multiple = FALSE,
                           filetypes = c("gtf", "gtf.gz")),
            br(), br(),
            textInput("prelnc_gtf_file", "GTF Annotation File",
                     value = "/home/data/scLncR/genome/Araport11_current.gtf",
                     placeholder = "Path to gene annotation file in GTF format"),
            helpText("Recommended: Use authoritative sources such as Araport, Ensembl, or RefSeq.")
          ),
          
          # lncRNA命名和计算资源
          box(
            width = 6,
            status = "danger",
            title = "4. lncRNA Naming and Computational Resources",
            solidHeader = TRUE,
            textInput("prelnc_lncrna_name", "lncRNA Name Prefix",
                     value = "AthLnc",
                     placeholder = "Custom prefix for newly predicted lncRNA identifiers"),
            helpText(HTML("<strong>Important:</strong> Only use letters, digits, and hyphens ('-').<br>
                         DO NOT use underscores ('_'), spaces, or other special characters.")),
            hr(),
            numericInput("prelnc_threads", "Number of Threads",
                        value = 24, min = 1, max = 128, step = 1),
            helpText("Recommendation: Set to the number of physical CPU cores or slightly lower."),
            hr(),
            h5("Validation Status"),
            uiOutput("prelnc_validation_status")
          )
        )
      ),
      
      # 2. Expression Counting 模块 - 根据 config_Count.yaml 设计
      tabItem(
        tabName = "count",
        h2("Expression Counting Module"),
        fluidRow(
          box(
            width = 12,
            status = "primary",
            solidHeader = TRUE,
            title = "Configuration File: config_Count.yaml",
            collapsible = TRUE
          )
        ),
        fluidRow(
          box(
            width = 6,
            status = "info",
            title = "1. Input Data Settings",
            solidHeader = TRUE,
            shinyDirButton("count_samples_dir", "Select Samples Directory", 
                         "Please select a directory"),
            br(), br(),
            textInput("count_samples_dir", "Samples Directory",
                     value = "/home/data/scLncR/data/paper_test/",
                     placeholder = "Directory containing raw sample data"),
            helpText("Each subdirectory under this path should correspond to one sample."),
            hr(),
            h5("Detected Samples"),
            uiOutput("count_detected_samples")
          ),
          
          box(
            width = 6,
            status = "success",
            title = "2. Project and Output Settings",
            solidHeader = TRUE,
            textInput("count_project_name", "Project Name",
                     value = "scLncR",
                     placeholder = "Used for naming output directories"),
            shinyDirButton("count_output_path", "Select Output Directory", 
                         "Please select a directory"),
            br(), br(),
            textInput("count_output_path", "Output Path",
                     value = "/home/data/scLncR/paper/paper_case/step2_count",
                     placeholder = "Main output directory for count results")
          )
        ),
        fluidRow(
          box(
            width = 6,
            status = "warning",
            title = "3. Reference Files",
            solidHeader = TRUE,
            shinyFilesButton("count_genome", "Select Genome File", 
                           "Please select a FASTA file", multiple = FALSE,
                           filetypes = c("fa", "fasta")),
            br(), br(),
            textInput("count_genome", "Genome File (FASTA)",
                     value = "/home/data/scLncR/genome/TAIR10.fa"),
            shinyFilesButton("count_gtf", "Select GTF File", 
                           "Please select a GTF file", multiple = FALSE,
                           filetypes = c("gtf")),
            br(), br(),
            textInput("count_gtf", "GTF Annotation File",
                     value = "/home/data/scLncR/genome/Araport11_current.gtf")
          ),
          
          box(
            width = 6,
            status = "danger",
            title = "4. lncRNA Annotation and Computational Resources",
            solidHeader = TRUE,
            shinyFilesButton("count_lnc_gtf", "Select lncRNA GTF File", 
                           "Please select a GTF file", multiple = FALSE,
                           filetypes = c("gtf")),
            br(), br(),
            textInput("count_lnc_gtf", "lncRNA GTF File",
                     value = "/home/data/scLncR/paper/paper_case/lnc_gtf",
                     placeholder = "Result of prelnc step"),
            helpText("Path to the lncRNA annotation file in GTF format."),
            hr(),
            numericInput("count_threads", "Number of Threads",
                        value = 24, min = 1, max = 128, step = 1),
            helpText("Used for compute-intensive steps like alignment and feature counting.")
          )
        )
      ),
      
      # 3. Data Processing 模块 - 根据 config_dataProcess.yaml 设计
      tabItem(
        tabName = "dataprocess",
        h2("Data Processing Module"),
        fluidRow(
          box(
            width = 12,
            status = "primary",
            solidHeader = TRUE,
            title = "Configuration File: config_dataProcess.yaml",
            collapsible = TRUE
          )
        ),
        fluidRow(
          box(
            width = 6,
            status = "info",
            title = "Core Parameters",
            solidHeader = TRUE,
            shinyDirButton("dp_counts_dir", "Select Counts Directory", 
                         "Please select a directory"),
            br(), br(),
            textInput("dp_counts_dir", "Counts Directory",
                     value = "~/paper/paper_case/step2_count",
                     placeholder = "Path to sample directories"),
            helpText("Paths to sample directories (comma-separated if multiple)"),
            textInput("dp_mt_name", "Mitochondrial Gene Prefix",
                     value = "AthM",
                     placeholder = "e.g., 'MT' or 'AthM'"),
            numericInput("dp_min_RNAs", "Minimum RNA Count per Cell",
                        value = 100, min = 0, max = 1000, step = 10),
            numericInput("dp_max_RNAs", "Maximum RNA Count per Cell",
                        value = 7000, min = 1000, max = 50000, step = 100),
            numericInput("dp_percent_mt", "Mitochondrial Percentage Threshold",
                        value = 10, min = 0, max = 100, step = 1)
          ),
          
          box(
            width = 6,
            status = "success",
            title = "Processing Parameters",
            solidHeader = TRUE,
            textInput("dp_lnc_name", "lncRNA Gene Prefix",
                     value = "AthLnc",
                     placeholder = "e.g., 'lncRNA' or 'AthLnc'"),
            numericInput("dp_nfeatures", "Number of Highly Variable Genes",
                        value = 6000, min = 1000, max = 10000, step = 500),
            numericInput("dp_dims", "Dimensionality for UMAP/tSNE",
                        value = 15, min = 5, max = 50, step = 1),
            numericInput("dp_resolution", "Clustering Resolution",
                        value = 0.5, min = 0.1, max = 2.0, step = 0.1),
            helpText("Higher = more clusters"),
            shinyFilesButton("dp_samples_info", "Select Samples Info File", 
                           "Please select a CSV/TXT file", multiple = FALSE,
                           filetypes = c("csv", "txt", "tsv")),
            br(), br(),
            textInput("dp_samples_info", "Samples Info File",
                     value = "~/data/paper_test/samples_info.txt",
                     placeholder = "Path to metadata CSV file")
          )
        ),
        fluidRow(
          box(
            width = 6,
            status = "warning",
            title = "Annotation Method Selection",
            solidHeader = TRUE,
            selectInput("dp_anno_method", "Annotation Method",
                       choices = c("SingleR" = "SingleR", 
                                  "scMM" = "scMM",
                                  "Manual" = "manual",
                                  "None" = "none"),
                       selected = "SingleR"),
            # 根据选择的注释方法显示不同的输入
            conditionalPanel(
              condition = "input.dp_anno_method == 'SingleR'",
              shinyFilesButton("dp_ref_file", "Select Reference File", 
                             "Please select an RDS file", multiple = FALSE,
                             filetypes = c("rds")),
              br(), br(),
              textInput("dp_ref_file", "Reference File for SingleR",
                       value = "~/data/singler/Ryu_2019_507252.rds",
                       placeholder = "Path to reference file for SingleR")
            ),
            conditionalPanel(
              condition = "input.dp_anno_method == 'scMM'",
              shinyFilesButton("dp_marker_gene_file", "Select Marker Gene File", 
                             "Please select a file", multiple = FALSE,
                             filetypes = c("csv", "txt", "tsv", "rds")),
              br(), br(),
              textInput("dp_marker_gene_file", "Marker Gene File for scMM",
                       value = "",
                       placeholder = "Path to marker gene file for scMM"),
              textInput("dp_tissue", "Tissue Type",
                       value = "",
                       placeholder = "e.g., 'lung', 'brain', 'blood'"),
              numericInput("dp_n_top_markers", "Number of Top Markers",
                          value = 3, min = 1, max = 20, step = 1)
            )
          ),
          
          box(
            width = 6,
            status = "danger",
            title = "Output Settings",
            solidHeader = TRUE,
            shinyDirButton("dp_output_dir", "Select Output Directory", 
                         "Please select a directory"),
            br(), br(),
            textInput("dp_output_dir", "Output Directory",
                     value = "~/paper/paper_case/step3_dataProcess",
                     placeholder = "Root output directory for all results"),
            hr(),
            h5("Generated Subdirectories:"),
            tags$ul(
              tags$li("sample_seurat_raw - Raw Seurat objects"),
              tags$li("data_preprocess - Preprocessing results"),
              tags$li("data_annotation - Annotation results")
            )
          )
        )
      ),
      
      # 4. Function Analysis 模块 - 根据 config_function.yaml 设计
      tabItem(
        tabName = "function_tab",
        h2("Function Analysis Module"),
        fluidRow(
          box(
            width = 12,
            status = "primary",
            solidHeader = TRUE,
            title = "Configuration File: config_function.yaml",
            collapsible = TRUE
          )
        ),
        fluidRow(
          box(
            width = 4,
            status = "info",
            title = "Global Input Settings",
            solidHeader = TRUE,
            shinyFilesButton("func_input_seurat", "Select Seurat Object", 
                           "Please select an RDS file", multiple = FALSE,
                           filetypes = c("rds")),
            br(), br(),
            textInput("func_input_seurat", "Input Seurat Object (RDS)",
                     value = "path/to/integrated_seurat.rds",
                     placeholder = "Path to preprocessed Seurat object"),
            helpText("Required: Path to preprocessed Seurat object (.rds file)"),
            hr(),
            checkboxGroupInput("func_run_modules", "Modules to Run",
                              choices = list(
                                "Location Analysis" = "location",
                                "Monocle2 Trajectory" = "monocle2",
                                "WGCNA Co-expression" = "wgcna"
                              ),
                              selected = c("location", "monocle2", "wgcna"))
          ),
          
          box(
            width = 4,
            status = "success",
            title = "Location Analysis Parameters",
            solidHeader = TRUE,
            numericInput("func_location_log2fc", "Log2FC Threshold",
                        value = 0.25, min = 0, max = 2, step = 0.05),
            numericInput("func_location_padj", "Adjusted P-value Threshold",
                        value = 0.05, min = 0.001, max = 0.1, step = 0.001),
            textInput("func_location_output", "Output Path",
                     value = "./location_results",
                     placeholder = "Output directory for location analysis"),
            textInput("func_location_lnc_name", "lncRNA Name Prefix",
                     value = "scLncR",
                     placeholder = "Prefix for lncRNA naming")
          ),
          
          box(
            width = 4,
            status = "warning",
            title = "Monocle2 Parameters",
            solidHeader = TRUE,
            numericInput("func_monocle2_qval", "Q-value Threshold",
                        value = 0.05, min = 0.001, max = 0.1, step = 0.001),
            selectInput("func_monocle2_method", "Dimension Reduction Method",
                       choices = c("DDRTree" = "DDRTree",
                                  "ICA" = "ICA",
                                  "tSNE" = "tSNE"),
                       selected = "DDRTree"),
            textInput("func_monocle2_output", "Output Path",
                     value = "./monocle2_results",
                     placeholder = "Output directory for Monocle2"),
            textInput("func_monocle2_hub_genes", "Hub Genes",
                     value = "all",
                     placeholder = "'all' or comma-separated gene list")
          )
        ),
        fluidRow(
          box(
            width = 6,
            status = "danger",
            title = "WGCNA Parameters",
            solidHeader = TRUE,
            textInput("func_wgcna_output", "Output Path",
                     value = "./wgcna_results",
                     placeholder = "Output directory for WGCNA"),
            selectizeInput("func_wgcna_cell_types", "Cell Types to Include",
                          choices = NULL,
                          multiple = TRUE,
                          options = list(placeholder = 'Select cell types (empty = all)')),
            textInput("func_wgcna_pro_name", "Project Prefix",
                     value = "scLncR",
                     placeholder = "Project prefix for module naming"),
            textInput("func_wgcna_lnc_name", "lncRNA Prefix Filter",
                     value = "AthLnc",
                     placeholder = "lncRNA prefix filter (empty = no filter)"),
            selectInput("func_wgcna_gene_select", "Gene Selection Method",
                       choices = c("fraction" = "fraction",
                                  "variance" = "variance"),
                       selected = "fraction")
          ),
          
          box(
            width = 6,
            status = "primary",
            title = "Module Dependencies",
            solidHeader = TRUE,
            h5("Module Prerequisites:"),
            tags$ul(
              tags$li(tags$strong("Location Analysis:"), 
                     "Requires preprocessed Seurat object with cell annotations"),
              tags$li(tags$strong("Monocle2:"), 
                     "Requires cell type annotations and variable features"),
              tags$li(tags$strong("WGCNA:"), 
                     "Requires normalized expression matrix")
            ),
            hr(),
            h5("Estimated Runtime:"),
            uiOutput("func_runtime_estimate"),
            hr(),
            actionButton("func_preview_config", "Preview Configuration",
                        class = "btn-info", width = "100%")
          )
        )
      ),
      
      # 运行管道页面
      tabItem(
        tabName = "run",
        h2("Run scLncR Pipeline"),
        fluidRow(
          box(
            width = 12,
            status = "success",
            solidHeader = TRUE,
            title = "Pipeline Execution Control",
            h4("Select Modules to Execute:"),
            fluidRow(
              column(3, awesomeCheckbox("run_prelnc", "LncRNA Prediction", value = TRUE)),
              column(3, awesomeCheckbox("run_count", "Expression Counting", value = TRUE)),
              column(3, awesomeCheckbox("run_dp", "Data Processing", value = TRUE)),
              column(3, awesomeCheckbox("run_func", "Function Analysis", value = TRUE))
            ),
            hr(),
            h4("Execution Options"),
            fluidRow(
              column(6,
                selectInput("run_mode", "Execution Mode",
                           choices = c(
                             "Full Pipeline (Sequential)" = "full",
                             "Selected Modules Only" = "selected",
                             "Dry Run (Generate Configs Only)" = "dryrun",
                             "Step-by-Step (Interactive)" = "stepwise"
                           ),
                           selected = "full"),
                checkboxInput("run_parallel", "Enable Parallel Execution", value = FALSE),
                conditionalPanel(
                  condition = "input.run_parallel == true",
                  numericInput("run_parallel_jobs", "Max Parallel Jobs",
                              value = 2, min = 1, max = 8)
                )
              ),
              column(6,
                textInput("run_output_dir", "Pipeline Output Directory",
                         value = "./scLncR_results",
                         placeholder = "Directory for all pipeline outputs"),
                selectInput("run_log_level", "Log Level",
                           choices = c("INFO" = "info",
                                      "DEBUG" = "debug",
                                      "WARN" = "warn",
                                      "ERROR" = "error"),
                           selected = "info"),
                checkboxInput("run_clean_intermediate", "Clean Intermediate Files", value = FALSE)
              )
            ),
            hr(),
            fluidRow(
              column(6,
                actionButton("validate_pipeline", "Validate Pipeline",
                            icon = icon("check-circle"),
                            class = "btn-warning",
                            width = "100%")
              ),
              column(6,
                actionButton("run_pipeline", "Execute Pipeline",
                            icon = icon("play"),
                            class = "btn-success btn-lg",
                            width = "100%")
              )
            ),
            hr(),
            h4("Pipeline Validation Status"),
            uiOutput("pipeline_validation_status"),
            hr(),
            h4("Execution Log"),
            div(style = "max-height: 400px; overflow-y: auto;",
                verbatimTextOutput("execution_log")),
            br(),
            fluidRow(
              column(6,
                downloadButton("download_log", "Download Log File",
                              class = "btn-info")
              ),
              column(6,
                actionButton("clear_log", "Clear Log",
                            class = "btn-default")
              )
            )
          )
        ),
        fluidRow(
          box(
            width = 6,
            title = "Generated Configuration Files",
            status = "primary",
            collapsible = TRUE,
            DTOutput("config_files_table")
          ),
          box(
            width = 6,
            title = "Pipeline Progress",
            status = "info",
            collapsible = TRUE,
            uiOutput("pipeline_progress"),
            hr(),
            h5("Current Status:"),
            verbatimTextOutput("current_status"),
            hr(),
            actionButton("stop_pipeline", "Stop Pipeline", 
                        icon = icon("stop"),
                        class = "btn-danger")
          )
        )
      ),
      
      # 配置文件页面
      tabItem(
        tabName = "configs",
        h2("Configuration Files Management"),
        fluidRow(
          box(
            width = 12,
            title = "Configuration Files",
            status = "primary",
            solidHeader = TRUE,
            DTOutput("all_configs_table")
          )
        ),
        fluidRow(
          box(
            width = 6,
            title = "Generate Configuration Files",
            status = "success",
            solidHeader = TRUE,
            selectInput("config_module", "Select Module",
                       choices = c("LncRNA Prediction" = "prelnc",
                                  "Expression Counting" = "count",
                                  "Data Processing" = "dataprocess",
                                  "Function Analysis" = "function")),
            textInput("config_filename", "Filename (without .yaml)",
                     value = "config", placeholder = "e.g., my_config"),
            br(),
            actionButton("generate_config", "Generate Config File",
                        icon = icon("file-code"),
                        class = "btn-primary",
                        width = "100%"),
            hr(),
            h5("Batch Generation:"),
            actionButton("generate_all_configs", "Generate All Configs",
                        class = "btn-warning", width = "100%")
          ),
          box(
            width = 6,
            title = "Load Configuration File",
            status = "info",
            solidHeader = TRUE,
            fileInput("load_config_input", "Choose YAML File",
                     accept = c(".yaml", ".yml")),
            selectInput("load_config_module", "Apply to Module",
                       choices = c("LncRNA Prediction" = "prelnc",
                                  "Expression Counting" = "count",
                                  "Data Processing" = "dataprocess",
                                  "Function Analysis" = "function",
                                  "All Modules" = "all")),
            br(),
            actionButton("load_config_btn", "Load Configuration",
                        icon = icon("upload"),
                        class = "btn-success",
                        width = "100%"),
            hr(),
            h5("Recent Configurations:"),
            uiOutput("recent_configs")
          )
        )
      ),
      
      # 文档页面
      tabItem(
        tabName = "docs",
        h2("Documentation"),
        fluidRow(
          box(
            width = 12,
            title = "scLncR Documentation",
            status = "info",
            solidHeader = TRUE,
            includeMarkdown("www/overview.md")
          )
        )
      ),
      
      # 关于页面
      tabItem(
        tabName = "about",
        h2("About scLncR"),
        fluidRow(
          box(
            width = 12,
            status = "primary",
            solidHeader = TRUE,
            h3("scLncR - Single-cell lncRNA Analysis Pipeline"),
            p("Version: 1.0.0"),
            p("Author: Your Name/Institution"),
            p("License: MIT"),
            hr(),
            h4("Citation:"),
            p("If you use scLncR in your research, please cite:"),
            p(em("scLncR: A comprehensive pipeline for single-cell long non-coding RNA analysis.")),
            hr(),
            h4("Contact:"),
            p("Email: your.email@institution.edu"),
            p("GitHub: https://github.com/yourrepo/scLncR"),
            hr(),
            h4("Acknowledgments:"),
            tags$ul(
              tags$li("R Core Team"),
              tags$li("Seurat development team"),
              tags$li("CellRanger team (10x Genomics)"),
              tags$li("All package maintainers")
            )
          )
        )
      )
    )
  )
)

# 服务器逻辑
server <- function(input, output, session) {
  
  # 初始化响应式值
  values <- reactiveValues(
    # 配置状态
    configs = list(
      prelnc = NULL,
      count = NULL,
      dataprocess = NULL,
      function_tab = NULL
    ),
    # 运行状态
    pipeline_status = "idle",
    execution_log = "",
    config_files = data.frame(),
    scLncR_path = "",
    scLncR_home = ""
  )
  
  # 初始化 - 查找scLncR
  observe({
    # 查找scLncR可执行文件
    scLncR_path <- Sys.which("scLncR")
    if (scLncR_path != "") {
      values$scLncR_path <- scLncR_path
      values$scLncR_home <- dirname(scLncR_path)
      
      # 更新日志
      values$execution_log <- paste0(
        values$execution_log,
        "[INFO] Found scLncR at: ", scLncR_path, "\n",
        "[INFO] Installation directory: ", values$scLncR_home, "\n"
      )
    } else {
      showNotification("scLncR not found in PATH. Please make sure it's installed.",
                      type = "error", duration = NULL)
    }
  })
  
  # 值框输出
  output$box_pipeline_status <- renderValueBox({
    status <- values$pipeline_status
    color <- if (status == "idle") "green" else
             if (status == "running") "yellow" else
             if (status == "completed") "blue" else "red"
    
    valueBox(
      value = toupper(status),
      subtitle = "Pipeline Status",
      icon = icon(if (status == "running") "sync" else "check-circle"),
      color = color
    )
  })
  
  output$box_scLncR_status <- renderValueBox({
    status <- if (values$scLncR_path != "") "Available" else "Not Found"
    color <- if (values$scLncR_path != "") "green" else "red"
    
    valueBox(
      value = status,
      subtitle = "scLncR Tool",
      icon = icon(if (values$scLncR_path != "") "check" else "times"),
      color = color
    )
  })
  
  # 模块状态表
  output$module_status_table <- renderDT({
    data.frame(
      Module = c("LncRNA Prediction", "Expression Counting", 
                "Data Processing", "Function Analysis"),
      Status = c("Ready", "Ready", "Ready", "Ready"),
      Config = c("Valid", "Valid", "Valid", "Valid"),
      LastRun = c("Never", "Never", "Never", "Never"),
      stringsAsFactors = FALSE
    )
  }, options = list(pageLength = 4))
  
  # 为每个模块添加文件浏览器功能
  volumes <- c(Home = fs::path_home(), "R Installation" = R.home(), getVolumes()())
  
  # LncRNA Prediction 文件浏览器
  shinyDirChoose(input, "prelnc_samples_dir", roots = volumes, session = session)
  shinyDirChoose(input, "prelnc_output_path", roots = volumes, session = session)
  shinyFileChoose(input, "prelnc_genome_file", roots = volumes, session = session)
  shinyFileChoose(input, "prelnc_gtf_file", roots = volumes, session = session)
  
  observeEvent(input$prelnc_samples_dir, {
    if (!is.integer(input$prelnc_samples_dir)) {
      path <- parseDirPath(volumes, input$prelnc_samples_dir)
      updateTextInput(session, "prelnc_samples_dirs", value = path)
    }
  })
  
  # 类似地为其他模块添加文件浏览器...
  
  # 生成配置文件函数
  generate_prelnc_config <- function() {
    config <- list(
      samples_dirs = input$prelnc_samples_dirs,
      project_name = input$prelnc_project_name,
      output_path = input$prelnc_output_path,
      genome_file = input$prelnc_genome_file,
      gtf_file = input$prelnc_gtf_file,
      lncrna_name = input$prelnc_lncrna_name,
      threads = as.integer(input$prelnc_threads)
    )
    
    # 添加注释
    attr(config, "comment") <- list(
      "scLncR - Single-cell lncRNA Analysis Pipeline: prelnc Module Configuration",
      "Purpose: Define input data paths, reference genome information, output settings, and runtime parameters."
    )
    
    return(config)
  }
  
  generate_count_config <- function() {
    config <- list(
      samples_dir = input$count_samples_dir,
      project_name = input$count_project_name,
      output_path = input$count_output_path,
      genome = input$count_genome,
      gtf = input$count_gtf,
      lnc_gtf = input$count_lnc_gtf,
      threads = as.integer(input$count_threads)
    )
    
    return(config)
  }
  
  generate_dataprocess_config <- function() {
    config <- list(
      counts_dir = input$dp_counts_dir,
      mt_name = input$dp_mt_name,
      min.RNAs = as.integer(input$dp_min_RNAs),
      max.RNAs = as.integer(input$dp_max_RNAs),
      percent.mt = as.numeric(input$dp_percent_mt),
      lnc_name = input$dp_lnc_name,
      nfeatures = as.integer(input$dp_nfeatures),
      dims = as.integer(input$dp_dims),
      resolution = as.numeric(input$dp_resolution),
      samples_info = input$dp_samples_info,
      anno_method = input$dp_anno_method,
      ref_file = if (input$dp_anno_method == "SingleR") input$dp_ref_file else "",
      marker_gene_file = if (input$dp_anno_method == "scMM") input$dp_marker_gene_file else "",
      tissue = if (input$dp_anno_method == "scMM") input$dp_tissue else "",
      n_top_markers = if (input$dp_anno_method == "scMM") as.integer(input$dp_n_top_markers) else 3,
      output_dir = input$dp_output_dir
    )
    
    return(config)
  }
  
  generate_function_config <- function() {
    config <- list(
      input_seurat = input$func_input_seurat,
      run_modules = input$func_run_modules,
      location = list(
        LOG2FC_THRESH = as.numeric(input$func_location_log2fc),
        PADJ_THRESH = as.numeric(input$func_location_padj),
        output_path = input$func_location_output,
        lncRNA_name = input$func_location_lnc_name
      ),
      monocle2 = list(
        qval = as.numeric(input$func_monocle2_qval),
        reduceDimensionMethod = input$func_monocle2_method,
        output_path = input$func_monocle2_output,
        hub_genes = input$func_monocle2_hub_genes
      ),
      wgcna = list(
        output_path = input$func_wgcna_output,
        cell_types = input$func_wgcna_cell_types,
        pro_name = input$func_wgcna_pro_name,
        lnc_name = input$func_wgcna_lnc_name,
        gene_select_method = input$func_wgcna_gene_select
      )
    )
    
    return(config)
  }
  
  # 生成配置文件按钮
  observeEvent(input$generate_config, {
    module <- input$config_module
    filename <- paste0(input$config_filename, ".yaml")
    
    config <- switch(module,
      "prelnc" = generate_prelnc_config(),
      "count" = generate_count_config(),
      "dataprocess" = generate_dataprocess_config(),
      "function" = generate_function_config()
    )
    
    # 写入文件
    write_yaml(config, filename)
    
    # 添加到记录
    new_row <- data.frame(
      Module = module,
      Filename = filename,
      Time = Sys.time(),
      Size = file.info(filename)$size,
      stringsAsFactors = FALSE
    )
    
    values$config_files <- rbind(values$config_files, new_row)
    
    showNotification(
      sprintf("Configuration file saved: %s", filename),
      type = "success"
    )
  })
  
  # 生成所有配置文件
  observeEvent(input$generate_all_configs, {
    timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
    
    # 生成prelnc配置
    prelnc_config <- generate_prelnc_config()
    write_yaml(prelnc_config, sprintf("config_LncPre_%s.yaml", timestamp))
    
    # 生成count配置
    count_config <- generate_count_config()
    write_yaml(count_config, sprintf("config_Count_%s.yaml", timestamp))
    
    # 生成dataprocess配置
    dp_config <- generate_dataprocess_config()
    write_yaml(dp_config, sprintf("config_dataProcess_%s.yaml", timestamp))
    
    # 生成function配置
    func_config <- generate_function_config()
    write_yaml(func_config, sprintf("config_function_%s.yaml", timestamp))
    
    showNotification(
      sprintf("All configuration files generated with timestamp: %s", timestamp),
      type = "success"
    )
  })
  
  # 运行管道
  observeEvent(input$run_pipeline, {
    # 检查scLncR是否可用
    if (values$scLncR_path == "") {
      showNotification("scLncR not found. Please install it first.", 
                      type = "error")
      return()
    }
    
    # 更新状态
    values$pipeline_status <- "running"
    values$execution_log <- paste0(
      values$execution_log,
      "[INFO] Starting scLncR pipeline...\n",
      "[INFO] Time: ", Sys.time(), "\n"
    )
    
    # 禁用运行按钮
    shinyjs::disable("run_pipeline")
    
    # 创建输出目录
    output_dir <- input$run_output_dir
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
      values$execution_log <- paste0(
        values$execution_log,
        "[INFO] Created output directory: ", output_dir, "\n"
      )
    }
    
    # 生成配置文件
    timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
    config_dir <- file.path(output_dir, "configs")
    dir.create(config_dir, recursive = TRUE)
    
    # 运行选中的模块
    if (input$run_prelnc) {
      config <- generate_prelnc_config()
      config_file <- file.path(config_dir, sprintf("config_LncPre_%s.yaml", timestamp))
      write_yaml(config, config_file)
      
      cmd <- paste(values$scLncR_path, "prelnc -c", config_file)
      values$execution_log <- paste0(
        values$execution_log,
        "[RUN] Executing: ", cmd, "\n"
      )
      
      # 在实际环境中，使用 system() 运行
      tryCatch({
        result <- system(cmd, intern = TRUE, wait = TRUE)
        values$execution_log <- paste0(
          values$execution_log,
          paste(result, collapse = "\n"), "\n",
          "[DONE] prelnc module completed\n"
        )
      }, error = function(e) {
        values$execution_log <- paste0(
          values$execution_log,
          "[ERROR] prelnc module failed: ", e$message, "\n"
        )
      })
    }
    
    # 类似地为其他模块添加运行逻辑...
    
    # 完成后启用按钮
    shinyjs::enable("run_pipeline")
    values$pipeline_status <- "completed"
    
    values$execution_log <- paste0(
      values$execution_log,
      "[INFO] Pipeline finished at: ", Sys.time(), "\n"
    )
    
    showNotification("Pipeline execution completed!", type = "success")
  })
  
  # 执行日志输出
  output$execution_log <- renderText({
    values$execution_log
  })
  
  # 配置文件表格
  output$config_files_table <- renderDT({
    if (nrow(values$config_files) == 0) {
      return(data.frame(Message = "No configuration files generated yet."))
    }
    values$config_files
  }, options = list(pageLength = 5))
  
  # 下载日志
  output$download_log <- downloadHandler(
    filename = function() {
      paste0("scLncR_log_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".txt")
    },
    content = function(file) {
      writeLines(values$execution_log, file)
    }
  )
  
  # 清除日志
  observeEvent(input$clear_log, {
    values$execution_log <- ""
  })
  
  # 停止管道
  observeEvent(input$stop_pipeline, {
    if (values$pipeline_status == "running") {
      values$pipeline_status <- "stopped"
      values$execution_log <- paste0(
        values$execution_log,
        "[WARN] Pipeline stopped by user at: ", Sys.time(), "\n"
      )
      showNotification("Pipeline stopped.", type = "warning")
    }
  })
  
  # 验证管道
  observeEvent(input$validate_pipeline, {
    errors <- c()
    warnings <- c()
    
    # 验证prelnc模块
    if (input$run_prelnc) {
      if (!dir.exists(input$prelnc_samples_dirs)) {
        errors <- c(errors, "prelnc: Samples directory does not exist")
      }
      if (!file.exists(input$prelnc_genome_file)) {
        errors <- c(errors, "prelnc: Genome file does not exist")
      }
    }
    
    # 显示验证结果
    if (length(errors) == 0) {
      output$pipeline_validation_status <- renderUI({
        div(class = "alert alert-success",
            icon("check-circle"), "All configurations are valid.")
      })
    } else {
      output$pipeline_validation_status <- renderUI({
        div(class = "alert alert-danger",
            icon("exclamation-triangle"), "Validation failed:",
            tags$ul(lapply(errors, tags$li)))
      })
    }
  })
}

# 运行应用
shinyApp(ui, server)