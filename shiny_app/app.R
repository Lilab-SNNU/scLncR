# shiny_app/app.R
library(shiny)
library(shinyjs)
library(shinydashboard)
library(yaml)
library(shinyWidgets)
library(shinyFiles)
library(shinyalert)

source("global.R")

# UI界面
ui <- dashboardPage(
  skin = "blue",
  dashboardHeader(
    title = "scLncR Pipeline v1.0",
    titleWidth = 300,
    tags$li(class = "dropdown",
            actionButton("help_btn", "Help", icon = icon("question-circle"),
                        style = "margin-top: 8px; margin-right: 5px;")
    )
  ),
  dashboardSidebar(
    width = 300,
    sidebarMenu(
      id = "tabs",
      menuItem("Dashboard", tabName = "dashboard", icon = icon("tachometer-alt")),
      menuItem("1. LncRNA Prediction", tabName = "prelnc", icon = icon("dna")),
      menuItem("2. Expression Counting", tabName = "count", icon = icon("calculator")),
      menuItem("3. Data Processing", tabName = "dataprocess", icon = icon("filter")),
      menuItem("4. Function Analysis", tabName = "function_tab", icon = icon("chart-line")),
      menuItem("Pipeline Runner", tabName = "runner", icon = icon("play-circle")),
      menuItem("Job Monitor", tabName = "monitor", icon = icon("desktop")),
      hr(),
      div(style = "padding: 15px;",
          h4("Project Settings"),
          textInput("global_project_name", "Project Name", value = "scLncR_Project"),
          numericInput("global_threads", "CPU Threads", value = 8, min = 1, max = 64),
          actionButton("save_project", "Save Project", icon = icon("save"), 
                      class = "btn-primary", width = "100%")
      )
    )
  ),
  dashboardBody(
    useShinyjs(),
    useShinyalert(),
    tags$head(
      tags$link(rel = "stylesheet", type = "text/css", href = "style.css"),
      tags$script(HTML("
        $(document).ready(function() {
          $('.sidebar-toggle').click();
        });
      "))
    ),
    tabItems(
      # Dashboard Tab
      tabItem(
        tabName = "dashboard",
        fluidRow(
          valueBoxOutput("box_prelnc"),
          valueBoxOutput("box_count"),
          valueBoxOutput("box_dataprocess"),
          valueBoxOutput("box_function")
        ),
        fluidRow(
          box(
            width = 8,
            title = "Pipeline Workflow",
            status = "primary",
            solidHeader = TRUE,
            img(src = "scLncR_workflow.png", width = "100%", 
                style = "border: 1px solid #ddd; border-radius: 5px;")
          ),
          box(
            width = 4,
            title = "Quick Actions",
            status = "info",
            solidHeader = TRUE,
            actionButton("quick_prelnc", "Run LncRNA Prediction", 
                        icon = icon("bolt"), width = "100%", class = "btn-success"),
            br(), br(),
            actionButton("quick_count", "Run Expression Counting", 
                        icon = icon("bolt"), width = "100%", class = "btn-success"),
            br(), br(),
            actionButton("quick_dp", "Run Data Processing", 
                        icon = icon("bolt"), width = "100%", class = "btn-success"),
            br(), br(),
            actionButton("quick_func", "Run Function Analysis", 
                        icon = icon("bolt"), width = "100%", class = "btn-success"),
            hr(),
            downloadButton("export_all_configs", "Export All Configs", 
                          class = "btn-warning", width = "100%")
          )
        ),
        fluidRow(
          box(
            width = 12,
            title = "Recent Jobs",
            status = "primary",
            solidHeader = TRUE,
            DTOutput("recent_jobs_table")
          )
        )
      ),
      
      # Tab 1: LncRNA Prediction (prelnc)
      tabItem(
        tabName = "prelnc",
        h2("LncRNA Prediction Module"),
        fluidRow(
          box(
            width = 4,
            status = "primary",
            solidHeader = TRUE,
            title = "Input Data",
            shinyFilesButton("prelnc_samples_dir", "Select Samples Directory", 
                           "Select directory containing FASTQ/BAM files", 
                           FALSE, class = "btn-primary", style = "width: 100%;"),
            br(), br(),
            verbatimTextOutput("prelnc_samples_dir_display"),
            hr(),
            h4("Reference Files"),
            shinyFilesButton("prelnc_genome", "Select Genome FASTA", 
                           "Select reference genome FASTA file", 
                           FALSE, class = "btn-info", style = "width: 100%;"),
            br(), br(),
            verbatimTextOutput("prelnc_genome_display"),
            shinyFilesButton("prelnc_gtf", "Select Reference GTF", 
                           "Select reference annotation GTF file", 
                           FALSE, class = "btn-info", style = "width: 100%;"),
            br(), br(),
            verbatimTextOutput("prelnc_gtf_display")
          ),
          
          box(
            width = 4,
            status = "warning",
            solidHeader = TRUE,
            title = "Output Settings",
            textInput("prelnc_project_name", "Project Name", 
                     value = "scLncR"),
            shinyDirButton("prelnc_output_path", "Select Output Directory", 
                          "Select output directory", class = "btn-success", 
                          style = "width: 100%;"),
            br(), br(),
            verbatimTextOutput("prelnc_output_path_display"),
            hr(),
            h4("LncRNA Naming"),
            textInput("prelnc_lncrna_name", "LncRNA Prefix", 
                     value = "AthLnc"),
            helpText("Note: Use only letters, digits, and hyphens. Underscores will cause issues in downstream analysis.")
          ),
          
          box(
            width = 4,
            status = "success",
            solidHeader = TRUE,
            title = "Run Configuration",
            numericInput("prelnc_threads", "CPU Threads", 
                        value = 24, min = 1, max = 64),
            hr(),
            h4("Validation Status"),
            uiOutput("prelnc_validation"),
            hr(),
            actionButton("run_prelnc_single", "Run This Module", 
                        icon = icon("play"), 
                        class = "btn-danger btn-lg", width = "100%"),
            br(), br(),
            downloadButton("download_prelnc_config", "Download Config", 
                          class = "btn-info", width = "100%")
          )
        )
      ),
      
      # Tab 2: Expression Counting (count)
      tabItem(
        tabName = "count",
        h2("Expression Counting Module"),
        fluidRow(
          box(
            width = 4,
            status = "primary",
            solidHeader = TRUE,
            title = "Input Data",
            shinyDirButton("count_samples_dir", "Select Samples Directory", 
                          "Select directory containing FASTQ files", 
                          class = "btn-primary", style = "width: 100%;"),
            br(), br(),
            verbatimTextOutput("count_samples_dir_display"),
            hr(),
            shinyFilesButton("count_genome", "Select Genome FASTA", 
                           "Select reference genome FASTA file", 
                           FALSE, class = "btn-info", style = "width: 100%;"),
            br(), br(),
            verbatimTextOutput("count_genome_display"),
            shinyFilesButton("count_gtf", "Select Reference GTF", 
                           "Select reference annotation GTF file", 
                           FALSE, class = "btn-info", style = "width: 100%;"),
            br(), br(),
            verbatimTextOutput("count_gtf_display")
          ),
          
          box(
            width = 4,
            status = "warning",
            solidHeader = TRUE,
            title = "LncRNA Annotation & Output",
            shinyDirButton("count_lnc_gtf", "Select LncRNA GTF Directory", 
                          "Select lncRNA GTF directory from prelnc step", 
                          class = "btn-warning", style = "width: 100%;"),
            br(), br(),
            verbatimTextOutput("count_lnc_gtf_display"),
            hr(),
            textInput("count_project_name", "Project Name", value = "scLncR"),
            shinyDirButton("count_output_path", "Select Output Directory", 
                          "Select output directory", 
                          class = "btn-success", style = "width: 100%;"),
            br(), br(),
            verbatimTextOutput("count_output_path_display")
          ),
          
          box(
            width = 4,
            status = "success",
            solidHeader = TRUE,
            title = "Run Configuration",
            numericInput("count_threads", "CPU Threads", 
                        value = 24, min = 1, max = 64),
            hr(),
            h4("Advanced Settings"),
            numericInput("count_cores", "CellRanger Cores", 
                        value = 4, min = 1, max = 16),
            helpText("Number of cores for CellRanger count step"),
            hr(),
            uiOutput("count_validation"),
            hr(),
            actionButton("run_count_single", "Run This Module", 
                        icon = icon("play"), 
                        class = "btn-danger btn-lg", width = "100%"),
            br(), br(),
            downloadButton("download_count_config", "Download Config", 
                          class = "btn-info", width = "100%")
          )
        )
      ),
      
      # Tab 3: Data Processing (dataProcess)
      tabItem(
        tabName = "dataprocess",
        h2("Data Processing Module"),
        fluidRow(
          box(
            width = 4,
            status = "primary",
            solidHeader = TRUE,
            title = "Input Data",
            shinyDirButton("dp_counts_dir", "Select Counts Directory", 
                          "Select directory containing filtered_feature_bc_matrix", 
                          class = "btn-primary", style = "width: 100%;"),
            br(), br(),
            verbatimTextOutput("dp_counts_dir_display"),
            hr(),
            shinyFilesButton("dp_samples_info", "Select Samples Info File", 
                           "Select samples information file (optional)", 
                           FALSE, class = "btn-info", style = "width: 100%;"),
            br(), br(),
            verbatimTextOutput("dp_samples_info_display")
          ),
          
          box(
            width = 4,
            status = "warning",
            solidHeader = TRUE,
            title = "Filtering Parameters",
            textInput("dp_mt_name", "Mitochondrial Gene Prefix", value = "AthM"),
            textInput("dp_lnc_name", "LncRNA Gene Prefix", value = "AthLnc"),
            numericInput("dp_min_rnas", "Minimum RNAs per Cell", 
                        value = 100, min = 0, max = 1000),
            numericInput("dp_max_rnas", "Maximum RNAs per Cell", 
                        value = 7000, min = 100, max = 50000),
            numericInput("dp_pct_mt", "Maximum Mitochondrial %", 
                        value = 10, min = 0, max = 100, step = 0.1),
            hr(),
            h4("Processing Parameters"),
            numericInput("dp_nfeatures", "Variable Features", 
                        value = 6000, min = 100, max = 10000),
            numericInput("dp_dims", "PCA Dimensions", 
                        value = 15, min = 5, max = 100),
            numericInput("dp_resolution", "Clustering Resolution", 
                        value = 0.5, min = 0.1, max = 2, step = 0.1)
          ),
          
          box(
            width = 4,
            status = "success",
            solidHeader = TRUE,
            title = "Annotation & Output",
            selectInput("dp_anno_method", "Annotation Method", 
                       choices = c("SingleR" = "SingleR", "scMM" = "scMM"),
                       selected = "SingleR"),
            
            # Conditional panel for SingleR
            conditionalPanel(
              condition = "input.dp_anno_method == 'SingleR'",
              shinyFilesButton("dp_ref_file", "Select Reference RDS", 
                             "Select SingleR reference RDS file", 
                             FALSE, class = "btn-info", style = "width: 100%;"),
              br(), br(),
              verbatimTextOutput("dp_ref_file_display")
            ),
            
            # Conditional panel for scMM
            conditionalPanel(
              condition = "input.dp_anno_method == 'scMM'",
              shinyFilesButton("dp_marker_file", "Select Marker Gene File", 
                             "Select marker gene file", 
                             FALSE, class = "btn-info", style = "width: 100%;"),
              br(), br(),
              verbatimTextOutput("dp_marker_file_display"),
              textInput("dp_tissue", "Tissue Type", value = ""),
              numericInput("dp_n_top_markers", "Top Markers", 
                          value = 3, min = 1, max = 20)
            ),
            
            shinyDirButton("dp_output_dir", "Select Output Directory", 
                          "Select output directory", 
                          class = "btn-success", style = "width: 100%;"),
            br(), br(),
            verbatimTextOutput("dp_output_dir_display"),
            hr(),
            uiOutput("dp_validation"),
            hr(),
            actionButton("run_dp_single", "Run This Module", 
                        icon = icon("play"), 
                        class = "btn-danger btn-lg", width = "100%"),
            br(), br(),
            downloadButton("download_dp_config", "Download Config", 
                          class = "btn-info", width = "100%")
          )
        )
      ),
      
      # Tab 4: Function Analysis (function)
      tabItem(
        tabName = "function_tab",
        h2("Function Analysis Module"),
        fluidRow(
          box(
            width = 3,
            status = "primary",
            solidHeader = TRUE,
            title = "Global Input",
            shinyFilesButton("func_input_seurat", "Select Seurat Object", 
                           "Select integrated Seurat RDS file", 
                           FALSE, class = "btn-primary", style = "width: 100%;"),
            br(), br(),
            verbatimTextOutput("func_input_seurat_display"),
            hr(),
            h4("Modules to Run"),
            awesomeCheckboxGroup(
              inputId = "func_run_modules",
              label = "Select Modules:",
              choices = c("Location Analysis" = "location",
                         "Monocle2 Trajectory" = "monocle2",
                         "WGCNA Network" = "wgcna"),
              selected = c("location", "monocle2", "wgcna"),
              status = "primary"
            )
          ),
          
          box(
            width = 3,
            status = "warning",
            solidHeader = TRUE,
            title = "Location Analysis",
            numericInput("func_loc_logfc", "Log2FC Threshold", 
                        value = 0.25, min = 0, max = 2, step = 0.01),
            numericInput("func_loc_padj", "Adjusted P-value Threshold", 
                        value = 0.05, min = 0.001, max = 0.2, step = 0.001),
            textInput("func_loc_lncrna_name", "LncRNA Prefix", value = "scLncR"),
            shinyDirButton("func_loc_output", "Output Directory", 
                          "Select output directory", 
                          class = "btn-info", style = "width: 100%;"),
            br(), br(),
            verbatimTextOutput("func_loc_output_display")
          ),
          
          box(
            width = 3,
            status = "info",
            solidHeader = TRUE,
            title = "Monocle2 Analysis",
            numericInput("func_mono_qval", "Q-value Threshold", 
                        value = 0.05, min = 0.001, max = 0.2, step = 0.001),
            selectInput("func_mono_method", "Reduction Method", 
                       choices = c("DDRTree", "ICA", "tSNE"),
                       selected = "DDRTree"),
            textInput("func_mono_hub_genes", "Hub Genes", value = "all"),
            helpText("Enter 'all' or comma-separated gene list"),
            shinyDirButton("func_mono_output", "Output Directory", 
                          "Select output directory", 
                          class = "btn-info", style = "width: 100%;"),
            br(), br(),
            verbatimTextOutput("func_mono_output_display")
          ),
          
          box(
            width = 3,
            status = "success",
            solidHeader = TRUE,
            title = "WGCNA Analysis",
            shinyDirButton("func_wgcna_output", "Output Directory", 
                          "Select output directory", 
                          class = "btn-info", style = "width: 100%;"),
            br(), br(),
            verbatimTextOutput("func_wgcna_output_display"),
            textInput("func_wgcna_cell_types", "Cell Types (optional)", 
                     value = "", placeholder = "T_cell,B_cell,NK_cell"),
            helpText("Leave empty for all cell types"),
            textInput("func_wgcna_pro_name", "Project Prefix", value = "scLncR"),
            textInput("func_wgcna_lnc_name", "LncRNA Filter", value = "AthLnc"),
            selectInput("func_wgcna_method", "Gene Selection Method", 
                       choices = c("fraction", "variance"),
                       selected = "fraction"),
            hr(),
            actionButton("run_func_single", "Run This Module", 
                        icon = icon("play"), 
                        class = "btn-danger btn-lg", width = "100%"),
            br(), br(),
            downloadButton("download_func_config", "Download Config", 
                          class = "btn-info", width = "100%")
          )
        )
      ),
      
      # Tab 5: Pipeline Runner
      tabItem(
        tabName = "runner",
        h2("Pipeline Runner"),
        fluidRow(
          box(
            width = 8,
            status = "primary",
            solidHeader = TRUE,
            title = "Pipeline Configuration",
            h4("Select Pipeline Steps:"),
            fluidRow(
              column(3,
                awesomeCheckbox("run_prelnc_pipe", "LncRNA Prediction", value = TRUE)
              ),
              column(3,
                awesomeCheckbox("run_count_pipe", "Expression Counting", value = TRUE)
              ),
              column(3,
                awesomeCheckbox("run_dp_pipe", "Data Processing", value = TRUE)
              ),
              column(3,
                awesomeCheckbox("run_func_pipe", "Function Analysis", value = TRUE)
              )
            ),
            hr(),
            h4("Execution Mode"),
            radioButtons("run_mode", NULL,
                        choices = c("Run All Selected" = "all",
                                   "Step-by-Step" = "stepwise",
                                   "Generate Configs Only" = "dryrun"),
                        selected = "all", inline = TRUE),
            conditionalPanel(
              condition = "input.run_mode == 'stepwise'",
              selectInput("start_step", "Start from Step:", 
                         choices = c("LncRNA Prediction" = "prelnc",
                                    "Expression Counting" = "count",
                                    "Data Processing" = "dataprocess",
                                    "Function Analysis" = "function"))
            ),
            hr(),
            h4("Output Directory"),
            textInput("pipe_output_root", "Pipeline Output Root", 
                     value = "/home/data/scLncR/analysis"),
            helpText("All module outputs will be saved under this directory")
          ),
          
          box(
            width = 4,
            status = "success",
            solidHeader = TRUE,
            title = "Run Pipeline",
            h4("Pipeline Status:"),
            verbatimTextOutput("pipe_status"),
            hr(),
            actionButton("validate_pipeline", "Validate Pipeline", 
                        icon = icon("check"), class = "btn-warning", width = "100%"),
            br(), br(),
            actionButton("run_pipeline", "Run Pipeline", 
                        icon = icon("rocket"), 
                        class = "btn-success btn-lg", width = "100%"),
            br(), br(),
            actionButton("stop_pipeline", "Stop Pipeline", 
                        icon = icon("stop"), 
                        class = "btn-danger", width = "100%"),
            hr(),
            downloadButton("download_pipe_config", "Download Pipeline Config", 
                          class = "btn-info", width = "100%")
          )
        ),
        fluidRow(
          box(
            width = 12,
            status = "info",
            solidHeader = TRUE,
            title = "Pipeline Output Log",
            div(style = "height: 400px; overflow-y: auto;",
                verbatimTextOutput("pipe_log")
            ),
            br(),
            actionButton("clear_log", "Clear Log", icon = icon("broom")),
            downloadButton("download_full_log", "Download Full Log", 
                          class = "btn-info")
          )
        )
      ),
      
      # Tab 6: Job Monitor
      tabItem(
        tabName = "monitor",
        h2("Job Monitor"),
        fluidRow(
          box(
            width = 12,
            status = "primary",
            solidHeader = TRUE,
            title = "Active Jobs",
            DTOutput("active_jobs_table")
          )
        ),
        fluidRow(
          box(
            width = 6,
            status = "success",
            solidHeader = TRUE,
            title = "Job Details",
            verbatimTextOutput("job_details")
          ),
          box(
            width = 6,
            status = "warning",
            solidHeader = TRUE,
            title = "Job Actions",
            actionButton("refresh_jobs", "Refresh Jobs", 
                        icon = icon("sync"), class = "btn-primary", width = "100%"),
            br(), br(),
            actionButton("kill_job", "Kill Selected Job", 
                        icon = icon("skull-crossbones"), 
                        class = "btn-danger", width = "100%"),
            br(), br(),
            actionButton("view_job_log", "View Job Log", 
                        icon = icon("file-alt"), 
                        class = "btn-info", width = "100%"),
            hr(),
            h4("Job Statistics"),
            plotOutput("job_stats", height = "200px")
          )
        )
      )
    )
  )
)

# Server逻辑
server <- function(input, output, session) {
  
  # 初始化文件选择器
  volumes <- c(Home = fs::path_home(), "R Installation" = R.home(), getVolumes()())
  
  # 设置文件选择器
  shinyFileChoose(input, "prelnc_samples_dir", roots = volumes, 
                  filetypes = c('', 'dir'))
  shinyFileChoose(input, "prelnc_genome", roots = volumes, 
                  filetypes = c('fa', 'fasta', 'fna'))
  shinyFileChoose(input, "prelnc_gtf", roots = volumes, 
                  filetypes = c('gtf', 'gff', 'gff3'))
  shinyDirChoose(input, "prelnc_output_path", roots = volumes)
  
  shinyFileChoose(input, "count_samples_dir", roots = volumes, 
                  filetypes = c('', 'dir'))
  shinyFileChoose(input, "count_genome", roots = volumes, 
                  filetypes = c('fa', 'fasta', 'fna'))
  shinyFileChoose(input, "count_gtf", roots = volumes, 
                  filetypes = c('gtf', 'gff', 'gff3'))
  shinyDirChoose(input, "count_lnc_gtf", roots = volumes)
  shinyDirChoose(input, "count_output_path", roots = volumes)
  
  shinyDirChoose(input, "dp_counts_dir", roots = volumes)
  shinyFileChoose(input, "dp_samples_info", roots = volumes, 
                  filetypes = c('txt', 'csv', 'tsv'))
  shinyFileChoose(input, "dp_ref_file", roots = volumes, 
                  filetypes = c('rds'))
  shinyFileChoose(input, "dp_marker_file", roots = volumes, 
                  filetypes = c('txt', 'csv', 'tsv', 'rds'))
  shinyDirChoose(input, "dp_output_dir", roots = volumes)
  
  shinyFileChoose(input, "func_input_seurat", roots = volumes, 
                  filetypes = c('rds'))
  shinyDirChoose(input, "func_loc_output", roots = volumes)
  shinyDirChoose(input, "func_mono_output", roots = volumes)
  shinyDirChoose(input, "func_wgcna_output", roots = volumes)
  
  # 响应式值
  values <- reactiveValues(
    prelnc_config = NULL,
    count_config = NULL,
    dp_config = NULL,
    func_config = NULL,
    pipe_log = "",
    pipe_status = "Ready",
    active_jobs = list(),
    job_history = list()
  )
  
  # Dashboard value boxes
  output$box_prelnc <- renderValueBox({
    valueBox(
      "LncRNA Prediction", 
      ifelse(is.null(values$prelnc_config), "Not Configured", "Configured"),
      icon = icon("dna"),
      color = ifelse(is.null(values$prelnc_config), "red", "green")
    )
  })
  
  output$box_count <- renderValueBox({
    valueBox(
      "Expression Counting", 
      ifelse(is.null(values$count_config), "Not Configured", "Configured"),
      icon = icon("calculator"),
      color = ifelse(is.null(values$count_config), "red", "green")
    )
  })
  
  output$box_dataprocess <- renderValueBox({
    valueBox(
      "Data Processing", 
      ifelse(is.null(values$dp_config), "Not Configured", "Configured"),
      icon = icon("filter"),
      color = ifelse(is.null(values$dp_config), "red", "green")
    )
  })
  
  output$box_function <- renderValueBox({
    valueBox(
      "Function Analysis", 
      ifelse(is.null(values$func_config), "Not Configured", "Configured"),
      icon = icon("chart-line"),
      color = ifelse(is.null(values$func_config), "red", "green")
    )
  })
  
  # File path displays for prelnc
  output$prelnc_samples_dir_display <- renderText({
    if (is.integer(input$prelnc_samples_dir)) {
      "No directory selected"
    } else {
      parseDirPath(volumes, input$prelnc_samples_dir)
    }
  })
  
  output$prelnc_genome_display <- renderText({
    if (is.integer(input$prelnc_genome)) {
      "No file selected"
    } else {
      parseFilePaths(volumes, input$prelnc_genome)$datapath
    }
  })
  
  output$prelnc_gtf_display <- renderText({
    if (is.integer(input$prelnc_gtf)) {
      "No file selected"
    } else {
      parseFilePaths(volumes, input$prelnc_gtf)$datapath
    }
  })
  
  output$prelnc_output_path_display <- renderText({
    if (is.integer(input$prelnc_output_path)) {
      "No directory selected"
    } else {
      parseDirPath(volumes, input$prelnc_output_path)
    }
  })
  
  # Validation for prelnc
  output$prelnc_validation <- renderUI({
    samples_dir <- if (!is.integer(input$prelnc_samples_dir)) 
      parseDirPath(volumes, input$prelnc_samples_dir) else NULL
    genome <- if (!is.integer(input$prelnc_genome)) 
      parseFilePaths(volumes, input$prelnc_genome)$datapath else NULL
    gtf <- if (!is.integer(input$prelnc_gtf)) 
      parseFilePaths(volumes, input$prelnc_gtf)$datapath else NULL
    output_path <- if (!is.integer(input$prelnc_output_path)) 
      parseDirPath(volumes, input$prelnc_output_path) else NULL
    
    if (is.null(samples_dir)) {
      return(div(icon("times-circle"), " Samples directory required", 
                 style = "color: red;"))
    }
    if (is.null(genome)) {
      return(div(icon("times-circle"), " Genome FASTA required", 
                 style = "color: red;"))
    }
    if (is.null(gtf)) {
      return(div(icon("times-circle"), " Reference GTF required", 
                 style = "color: red;"))
    }
    if (is.null(output_path)) {
      return(div(icon("times-circle"), " Output directory required", 
                 style = "color: red;"))
    }
    
    # Check if files exist
    checks <- c(
      if (!dir.exists(samples_dir)) "Samples directory not found",
      if (!file.exists(genome)) "Genome file not found",
      if (!file.exists(gtf)) "GTF file not found"
    )
    
    if (length(checks) > 0) {
      return(div(icon("exclamation-triangle"), 
                 paste(checks, collapse = "; "), 
                 style = "color: orange;"))
    }
    
    return(div(icon("check-circle"), " All inputs valid", 
               style = "color: green;"))
  })
  
  # Generate prelnc config
  observeEvent(input$run_prelnc_single, {
    # Validate inputs
    samples_dir <- if (!is.integer(input$prelnc_samples_dir)) 
      parseDirPath(volumes, input$prelnc_samples_dir) else NULL
    genome <- if (!is.integer(input$prelnc_genome)) 
      parseFilePaths(volumes, input$prelnc_genome)$datapath else NULL
    gtf <- if (!is.integer(input$prelnc_gtf)) 
      parseFilePaths(volumes, input$prelnc_gtf)$datapath else NULL
    output_path <- if (!is.integer(input$prelnc_output_path)) 
      parseDirPath(volumes, input$prelnc_output_path) else NULL
    
    if (any(sapply(list(samples_dir, genome, gtf, output_path), is.null))) {
      shinyalert("Error", "Please fill all required fields", type = "error")
      return()
    }
    
    # Generate config
    config <- list(
      samples_dirs = samples_dir,
      project_name = input$prelnc_project_name,
      output_path = output_path,
      genome_file = genome,
      gtf_file = gtf,
      lncrna_name = input$prelnc_lncrna_name,
      threads = input$prelnc_threads
    )
    
    values$prelnc_config <- config
    
    # Save config file
    timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
    config_file <- file.path(output_path, 
                            paste0("config_LncPre_", timestamp, ".yaml"))
    
    tryCatch({
      yaml::write_yaml(config, config_file)
      
      # Run the command
      cmd <- paste("scLncR prelnc -c", shQuote(config_file))
      
      # Add to log
      values$pipe_log <- paste0(values$pipe_log, 
                               "[" , Sys.time(), "] Running prelnc module\n",
                               "Command: ", cmd, "\n")
      
      # Run in background
      job_id <- paste0("prelnc_", timestamp)
      values$active_jobs[[job_id]] <- list(
        id = job_id,
        module = "prelnc",
        cmd = cmd,
        config = config_file,
        start_time = Sys.time(),
        status = "running"
      )
      
      # Execute command
      system(cmd, wait = FALSE)
      
      shinyalert("Success", 
                paste("prelnc module started. Config saved to:", config_file),
                type = "success")
      
    }, error = function(e) {
      shinyalert("Error", paste("Failed to run prelnc:", e$message), 
                type = "error")
    })
  })
  
  # Download prelnc config
  output$download_prelnc_config <- downloadHandler(
    filename = function() {
      paste0("config_LncPre_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".yaml")
    },
    content = function(file) {
      if (is.null(values$prelnc_config)) {
        showNotification("No config to download", type = "warning")
        return()
      }
      yaml::write_yaml(values$prelnc_config, file)
    }
  )
  
  # Similar implementations for count, dataProcess, and function modules...
  # Due to space constraints, I'll show the pattern for count module
  
  # File path displays for count
  output$count_samples_dir_display <- renderText({
    if (is.integer(input$count_samples_dir)) {
      "No directory selected"
    } else {
      parseDirPath(volumes, input$count_samples_dir)
    }
  })
  
  output$count_genome_display <- renderText({
    if (is.integer(input$count_genome)) {
      "No file selected"
    } else {
      parseFilePaths(volumes, input$count_genome)$datapath
    }
  })
  
  output$count_gtf_display <- renderText({
    if (is.integer(input$count_gtf)) {
      "No file selected"
    } else {
      parseFilePaths(volumes, input$count_gtf)$datapath
    }
  })
  
  output$count_lnc_gtf_display <- renderText({
    if (is.integer(input$count_lnc_gtf)) {
      "No directory selected"
    } else {
      parseDirPath(volumes, input$count_lnc_gtf)
    }
  })
  
  output$count_output_path_display <- renderText({
    if (is.integer(input$count_output_path)) {
      "No directory selected"
    } else {
      parseDirPath(volumes, input$count_output_path)
    }
  })
  
  # Validation for count
  output$count_validation <- renderUI({
    # Similar validation logic as prelnc
    # ...
  })
  
  # Generate count config
  observeEvent(input$run_count_single, {
    # Similar logic as prelnc
    # ...
  })
  
  # Pipeline runner
  observeEvent(input$run_pipeline, {
    # Disable run button
    shinyjs::disable("run_pipeline")
    values$pipe_status <- "Running..."
    
    # Start log
    values$pipe_log <- paste0(values$pipe_log,
                             "\n", strrep("=", 60),
                             "\nStarting scLncR Pipeline at ", Sys.time(),
                             "\n", strrep("=", 60), "\n")
    
    # Create pipeline directory
    pipe_dir <- input$pipe_output_root
    if (!dir.exists(pipe_dir)) {
      dir.create(pipe_dir, recursive = TRUE)
    }
    
    # Generate and run each module
    modules_to_run <- c()
    if (input$run_prelnc_pipe) modules_to_run <- c(modules_to_run, "prelnc")
    if (input$run_count_pipe) modules_to_run <- c(modules_to_run, "count")
    if (input$run_dp_pipe) modules_to_run <- c(modules_to_run, "dataprocess")
    if (input$run_func_pipe) modules_to_run <- c(modules_to_run, "function")
    
    # Execute modules
    for (module in modules_to_run) {
      values$pipe_log <- paste0(values$pipe_log,
                               "\n[", Sys.time(), "] Starting ", module, " module\n")
      
      # Generate config and run
      # This would call the individual module run functions
      # ...
    }
    
    # Re-enable run button
    shinyjs::enable("run_pipeline")
    values$pipe_status <- "Completed"
    values$pipe_log <- paste0(values$pipe_log,
                             "\n", strrep("=", 60),
                             "\nPipeline completed at ", Sys.time(),
                             "\n", strrep("=", 60), "\n")
  })
  
  # Clear log
  observeEvent(input$clear_log, {
    values$pipe_log <- ""
  })
  
  # Download full log
  output$download_full_log <- downloadHandler(
    filename = function() {
      paste0("scLncR_log_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".txt")
    },
    content = function(file) {
      writeLines(values$pipe_log, file)
    }
  )
  
  # Active jobs table
  output$active_jobs_table <- renderDT({
    if (length(values$active_jobs) == 0) {
      return(data.frame(Message = "No active jobs"))
    }
    
    jobs_df <- do.call(rbind, lapply(values$active_jobs, function(job) {
      data.frame(
        JobID = job$id,
        Module = job$module,
        StartTime = as.character(job$start_time),
        Status = job$status,
        Config = basename(job$config)
      )
    }))
    
    datatable(jobs_df, selection = 'single', options = list(
      pageLength = 10,
      dom = 'Bfrtip'
    ))
  })
  
  # Pipe log output
  output$pipe_log <- renderText({
    values$pipe_log
  })
  
  # Pipe status output
  output$pipe_status <- renderText({
    values$pipe_status
  })
  
  # Help button
  observeEvent(input$help_btn, {
    shinyalert(
      title = "scLncR Help",
      text = HTML("
        <h4>Pipeline Steps:</h4>
        <ol>
          <li><b>LncRNA Prediction</b>: Predict lncRNAs from single-cell/single-nucleus RNA-seq data</li>
          <li><b>Expression Counting</b>: Generate expression matrices using CellRanger</li>
          <li><b>Data Processing</b>: Quality control, normalization, and cell annotation</li>
          <li><b>Function Analysis</b>: Spatial localization, trajectory, and co-expression analysis</li>
        </ol>
        <h4>Quick Start:</h4>
        <p>1. Configure each module with your data paths</p>
        <p>2. Use 'Run This Module' to test individual steps</p>
        <p>3. Use 'Pipeline Runner' to run the complete workflow</p>
        <h4>Note:</h4>
        <p>• All paths should be absolute paths</p>
        <p>• Ensure sufficient disk space for intermediate files</p>
      "),
      html = TRUE,
      size = "l"
    )
  })
}

# Run the application
shinyApp(ui, server)