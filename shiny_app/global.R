# shiny_app/global.R
# 全局设置和辅助函数

# 加载必要的库
library(shiny)
library(shinydashboard)
library(shinyjs)
library(shinyWidgets)
library(shinyFiles)
library(yaml)
library(DT)
library(fs)
library(processx)
library(future)
library(promises)
library(bs4Dash)

# 检查必要的包
check_required_packages <- function() {
  required_packages <- c(
    "shiny", "shinydashboard", "shinyjs", "shinyWidgets", 
    "shinyFiles", "yaml", "DT", "fs", "processx"
  )
  
  missing_packages <- required_packages[!sapply(required_packages, requireNamespace, quietly = TRUE)]
  
  if (length(missing_packages) > 0) {
    warning("The following packages are missing: ", 
            paste(missing_packages, collapse = ", "))
    return(FALSE)
  }
  return(TRUE)
}

# 初始化文件系统
init_filesystem <- function() {
  # 创建必要的目录
  dirs_to_create <- c(
    "configs",
    "logs",
    "results",
    "temp"
  )
  
  for (dir in dirs_to_create) {
    if (!dir.exists(dir)) {
      dir.create(dir, recursive = TRUE, showWarnings = FALSE)
    }
  }
}

# 生成模块配置
generate_prelnc_config <- function(params) {
  config <- list(
    samples_dirs = params$samples_dirs,
    project_name = params$project_name,
    output_path = params$output_path,
    genome_file = params$genome_file,
    gtf_file = params$gtf_file,
    lncrna_name = params$lncrna_name,
    threads = params$threads
  )
  return(config)
}

generate_count_config <- function(params) {
  config <- list(
    samples_dir = params$samples_dir,
    project_name = params$project_name,
    output_path = params$output_path,
    genome = params$genome,
    gtf = params$gtf,
    lnc_gtf = params$lnc_gtf,
    threads = params$threads
  )
  return(config)
}

generate_dataprocess_config <- function(params) {
  config <- list(
    counts_dir = params$counts_dir,
    mt_name = params$mt_name,
    min.RNAs = params$min.RNAs,
    max.RNAs = params$max.RNAs,
    percent.mt = params$percent.mt,
    lnc_name = params$lnc_name,
    nfeatures = params$nfeatures,
    dims = params$dims,
    resolution = params$resolution,
    samples_info = params$samples_info,
    anno_method = params$anno_method,
    ref_file = params$ref_file,
    marker_gene_file = params$marker_gene_file,
    tissue = params$tissue,
    n_top_markers = params$n_top_markers,
    output_dir = params$output_dir
  )
  return(config)
}

generate_function_config <- function(params) {
  config <- list(
    input_seurat = params$input_seurat,
    run_modules = params$run_modules,
    location = list(
      LOG2FC_THRESH = params$loc_logfc_thresh,
      PADJ_THRESH = params$loc_padj_thresh,
      output_path = params$loc_output_path,
      lncRNA_name = params$loc_lncrna_name
    ),
    monocle2 = list(
      qval = params$mono_qval,
      reduceDimensionMethod = params$mono_method,
      output_path = params$mono_output_path,
      hub_genes = params$mono_hub_genes
    ),
    wgcna = list(
      output_path = params$wgcna_output_path,
      cell_types = if (params$wgcna_cell_types == "") list() else 
        strsplit(params$wgcna_cell_types, ",")[[1]],
      pro_name = params$wgcna_pro_name,
      lnc_name = params$wgcna_lnc_name,
      gene_select_method = params$wgcna_method
    )
  )
  return(config)
}

# 验证文件路径
validate_path <- function(path, must_exist = TRUE, is_dir = FALSE) {
  if (is.null(path) || path == "") {
    return(FALSE)
  }
  
  if (must_exist) {
    if (is_dir) {
      return(dir.exists(path))
    } else {
      return(file.exists(path))
    }
  }
  return(TRUE)
}

# 运行命令
run_command <- function(cmd, log_file = NULL, background = FALSE) {
  if (is.null(log_file)) {
    log_file <- tempfile("scLncR_log_", fileext = ".txt")
  }
  
  if (background) {
    # 在后台运行
    proc <- processx::process$new(
      command = "bash",
      args = c("-c", cmd),
      stdout = log_file,
      stderr = log_file
    )
    return(proc)
  } else {
    # 在前台运行并捕获输出
    result <- system(cmd, intern = TRUE, wait = TRUE)
    writeLines(result, log_file)
    return(list(
      success = attr(result, "status") == 0,
      output = result,
      log_file = log_file
    ))
  }
}

# 检查scLncR是否可用
check_scLncR <- function() {
  cmd <- "which scLncR"
  result <- try(system(cmd, intern = TRUE, ignore.stderr = TRUE), silent = TRUE)
  
  if (length(result) == 0 || result == "") {
    return(FALSE)
  }
  return(TRUE)
}

# 获取scLncR版本
get_scLncR_version <- function() {
  if (!check_scLncR()) {
    return("Not installed")
  }
  
  cmd <- "scLncR --version 2>&1 || scLncR -h 2>&1 | head -1"
  result <- try(system(cmd, intern = TRUE), silent = TRUE)
  
  if (length(result) > 0) {
    return(result[1])
  }
  return("Unknown")
}

# 解析样本信息
parse_samples_info <- function(file_path) {
  if (!file.exists(file_path)) {
    return(NULL)
  }
  
  # 尝试读取不同的格式
  if (grepl("\\.csv$", file_path)) {
    return(read.csv(file_path, stringsAsFactors = FALSE))
  } else if (grepl("\\.tsv$", file_path)) {
    return(read.delim(file_path, stringsAsFactors = FALSE))
  } else if (grepl("\\.txt$", file_path)) {
    return(read.table(file_path, header = TRUE, stringsAsFactors = FALSE))
  }
  return(NULL)
}

# 创建项目结构
create_project_structure <- function(project_name, output_root) {
  project_dirs <- c(
    file.path(output_root, project_name, "00_rawdata"),
    file.path(output_root, project_name, "01_prelnc"),
    file.path(output_root, project_name, "02_count"),
    file.path(output_root, project_name, "03_dataprocess"),
    file.path(output_root, project_name, "04_function"),
    file.path(output_root, project_name, "05_reports"),
    file.path(output_root, project_name, "configs"),
    file.path(output_root, project_name, "logs")
  )
  
  for (dir in project_dirs) {
    if (!dir.exists(dir)) {
      dir.create(dir, recursive = TRUE, showWarnings = FALSE)
    }
  }
  
  return(project_dirs)
}

# 初始化应用
.onAttach <- function(libname, pkgname) {
  # 检查必要的包
  check_required_packages()
  
  # 初始化文件系统
  init_filesystem()
  
  # 检查scLncR
  if (!check_scLncR()) {
    packageStartupMessage(
      "WARNING: scLncR command-line tool not found in PATH.\n",
      "Make sure it's installed and accessible for pipeline execution."
    )
  } else {
    version <- get_scLncR_version()
    packageStartupMessage(paste("scLncR version:", version))
  }
}