# shiny_app/global.R
# 全局设置和工具函数

resolve_scLncR_home <- function() {
  env_home <- Sys.getenv("SC_LNCR_HOME", unset = "")
  if (nzchar(env_home) && dir.exists(env_home)) return(normalizePath(env_home))
  normalizePath(file.path(getwd(), ".."), mustWork = FALSE)
}

scLncR_home <- resolve_scLncR_home()
package_loader <- file.path(scLncR_home, "R", "modules", "scLncR_load_package.R")
if (!file.exists(package_loader)) {
  stop("Package loader not found: ", package_loader)
}
source(package_loader, local = TRUE)
scLncR_load_shiny_packages()
.sclncr_shiny_global_loaded <- TRUE

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0) return(y)
  x
}

safe_include_markdown <- function(path) {
  if (requireNamespace("markdown", quietly = TRUE)) {
    return(includeMarkdown(path))
  }
  if (!file.exists(path)) {
    return(tags$div(class = "alert alert-warning", paste("Documentation file not found:", path)))
  }
  tags$pre(
    style = "white-space: pre-wrap; background: #fff; border: 0; font-family: inherit;",
    paste(readLines(path, warn = FALSE), collapse = "\n")
  )
}

# 获取scLncR安装信息
get_scLncR_info <- function() {
  scLncR_path <- Sys.which("scLncR")
  
  if (scLncR_path == "") {
    return(list(
      found = FALSE,
      path = "",
      home = "",
      message = "scLncR not found in PATH. Please install it first."
    ))
  }
  
  scLncR_home <- dirname(scLncR_path)
  
  # 检查目录结构
  expected_dirs <- c("R", "shiny_app", "figures")
  dir_exists <- sapply(file.path(scLncR_home, expected_dirs), dir.exists)
  
  if (!all(dir_exists)) {
    return(list(
      found = TRUE,
      path = scLncR_path,
      home = scLncR_home,
      message = "scLncR found but directory structure appears incomplete.",
      valid = FALSE
    ))
  }
  
  return(list(
    found = TRUE,
    path = scLncR_path,
    home = scLncR_home,
    message = "scLncR found and directory structure is valid.",
    valid = TRUE
  ))
}

# 验证配置参数
validate_prelnc_config <- function(config) {
  errors <- c()
  warnings <- c()
  
  input_type <- tolower(config$input_type %||% "fastq")

  # 检查必要的参数
  required_fields <- c("input_type", "samples_dirs", "genome_file", "gtf_file", "output_path", "lncrna_name", "sequencing_platform")
  missing_fields <- setdiff(required_fields, names(config))
  if (length(missing_fields) > 0) {
    errors <- c(errors, paste("Missing required fields:", 
                             paste(missing_fields, collapse = ", ")))
  }
  
  if (!input_type %in% c("fastq")) {
    errors <- c(errors, "input_type must be fastq in current prelnc design")
  }
  if (!tolower(config$sequencing_platform %||% "") %in% c("10x", "smartseq2", "dropseq", "generic")) {
    errors <- c(errors, "sequencing_platform must be one of: 10x, smartseq2, dropseq, generic")
  }

  if (is.null(config$samples_dirs) || !nzchar(config$samples_dirs)) {
    errors <- c(errors, "samples_dirs is required")
  } else if (!dir.exists(config$samples_dirs)) {
    warnings <- c(warnings, "samples_dirs does not exist")
  }
  
  if (!file.exists(config$genome_file)) {
    errors <- c(errors, "Genome file does not exist")
  }
  
  if (!file.exists(config$gtf_file)) {
    warnings <- c(warnings, "GTF file does not exist")
  }
  
  # 检查lncRNA命名
  if (!grepl("^[A-Za-z0-9-]+$", config$lncrna_name)) {
    errors <- c(errors, "lncRNA name prefix contains invalid characters. Use only letters, digits, and hyphens.")
  }
  
  return(list(errors = errors, warnings = warnings))
}

validate_count_config <- function(config) {
  errors <- c()
  warnings <- c()
  
  required_fields <- c("sequencing_platform", "count_engine", "samples_dirs", "genome_file", "known_gtf", "lnc_gtf", "output_path")
  missing_fields <- setdiff(required_fields, names(config))
  if (length(missing_fields) > 0) {
    errors <- c(errors, paste("Missing required fields:", 
                             paste(missing_fields, collapse = ", ")))
  }
  
  if (!dir.exists(config$samples_dirs)) {
    warnings <- c(warnings, "Samples directory does not exist")
  }
  
  if (!file.exists(config$lnc_gtf)) {
    warnings <- c(warnings, "lncRNA GTF file does not exist (may not be generated yet)")
  }
  if (!is.null(config$combined_gtf) && nzchar(config$combined_gtf) && !file.exists(config$combined_gtf)) {
    warnings <- c(warnings, "combined_gtf does not exist (will be rebuilt from known_gtf + lnc_gtf)")
  }
  if (!tolower(config$sequencing_platform %||% "") %in% c("10x", "smartseq2", "dropseq", "generic")) {
    errors <- c(errors, "sequencing_platform must be one of: 10x, smartseq2, dropseq, generic")
  }
  if (!tolower(config$count_engine %||% "") %in% c("cellranger", "featurecounts", "starsolo")) {
    errors <- c(errors, "count_engine must be one of: cellranger, featurecounts, starsolo")
  }
  
  return(list(errors = errors, warnings = warnings))
}

validate_dataprocess_config <- function(config) {
  errors <- c()
  warnings <- c()
  
  # 检查必要的目录
  if (!dir.exists(config$counts_dir)) {
    warnings <- c(warnings, "Counts directory does not exist")
  }
  
  # 检查注释方法相关参数
  if (config$anno_method == "SingleR" && !file.exists(config$ref_file)) {
    errors <- c(errors, "SingleR reference file not found")
  }
  
  if (config$anno_method == "scMM") {
    if (config$tissue == "") {
      errors <- c(errors, "Tissue type must be specified for scMM")
    }
    if (!file.exists(config$marker_gene_file)) {
      warnings <- c(warnings, "Marker gene file not found for scMM")
    }
  }
  
  return(list(errors = errors, warnings = warnings))
}

# 创建目录并设置权限
create_output_directories <- function(base_dir, subdirs = NULL) {
  if (!dir.exists(base_dir)) {
    dir.create(base_dir, recursive = TRUE, mode = "0755")
  }
  
  if (!is.null(subdirs)) {
    for (subdir in subdirs) {
      full_path <- file.path(base_dir, subdir)
      if (!dir.exists(full_path)) {
        dir.create(full_path, recursive = TRUE, mode = "0755")
      }
    }
  }
  
  return(base_dir)
}

# 安全读取YAML配置文件
safe_read_yaml <- function(filepath) {
  if (!file.exists(filepath)) {
    stop(paste("Config file not found:", filepath))
  }
  
  tryCatch({
    yaml::read_yaml(filepath)
  }, error = function(e) {
    stop(paste("Error reading YAML file:", e$message))
  })
}

# 安全写入YAML配置文件
safe_write_yaml <- function(data, filepath) {
  tryCatch({
    yaml::write_yaml(data, filepath)
    return(TRUE)
  }, error = function(e) {
    warning(paste("Error writing YAML file:", e$message))
    return(FALSE)
  })
}

# 获取系统资源信息
get_system_resources <- function() {
  list(
    cpus = parallel::detectCores(),
    memory = as.numeric(system("grep MemTotal /proc/meminfo | awk '{print $2}'", intern = TRUE)) / 1024 / 1024,
    disk = as.numeric(system("df . | tail -1 | awk '{print $4}'", intern = TRUE)) / 1024 / 1024
  )
}

# 生成唯一的配置文件名
generate_config_filename <- function(module, prefix = "config") {
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  module_map <- c(
    "prelnc" = "LncPre",
    "count" = "Count",
    "dataprocess" = "dataProcess",
    "function" = "function"
  )
  
  if (module %in% names(module_map)) {
    return(paste0(prefix, "_", module_map[module], "_", timestamp, ".yaml"))
  } else {
    return(paste0(prefix, "_", module, "_", timestamp, ".yaml"))
  }
}

# 加载示例配置
load_example_config <- function(module) {
  examples <- list(
    prelnc = list(
      input_type = "fastq",
      sequencing_platform = "10x",
      samples_dirs = "/path/to/your/fastq_dir",
      samples_info = "/path/to/samples_info.txt",
      read_layout = "auto",
      project_name = "example_project",
      output_path = "/path/to/output",
      genome_file = "/path/to/genome.fa",
      gtf_file = "/path/to/annotation.gtf",
      lncrna_name = "exampleLnc",
      aligner_for_discovery = "hisat2",
      assembler = "stringtie",
      threads = 4,
      dry_run = TRUE
    ),
    count = list(
      sequencing_platform = "10x",
      count_engine = "cellranger",
      samples_dirs = "/path/to/your/samples",
      project_name = "example_project",
      output_path = "/path/to/output/count",
      genome_file = "/path/to/genome.fa",
      known_gtf = "/path/to/annotation.gtf",
      lnc_gtf = "/path/to/lnc_annotation.gtf",
      combined_gtf = "/path/to/combined_mRNA_lncRNA.gtf",
      threads = 4
    )
    # 其他模块的示例配置...
  )
  
  if (module %in% names(examples)) {
    return(examples[[module]])
  } else {
    return(NULL)
  }
}

# 检查文件路径是否可写
is_path_writable <- function(path) {
  if (file.exists(path)) {
    return(file.access(path, 2) == 0)
  } else {
    # 检查父目录是否可写
    parent_dir <- dirname(path)
    if (dir.exists(parent_dir)) {
      return(file.access(parent_dir, 2) == 0)
    } else {
      return(FALSE)
    }
  }
}

# 初始化函数
.initialize <- function() {
  # 检查scLncR
  scLncR_info <- get_scLncR_info()
  
  # 创建必要的目录
  app_dirs <- c("configs", "logs", "tmp", "outputs")
  for (dir in app_dirs) {
    if (!dir.exists(dir)) {
      dir.create(dir, recursive = TRUE, showWarnings = FALSE)
    }
  }
  
  # 返回初始化状态
  return(list(
    scLncR = scLncR_info,
    app_dirs = app_dirs,
    initialized = TRUE
  ))
}

# 执行初始化
init_status <- .initialize()
