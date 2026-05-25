# shiny_app/run_shiny.R
# 简化的启动脚本

#!/usr/bin/env Rscript

cat("========================================\n")
cat("   scLncR Shiny GUI\n")
cat("========================================\n\n")

get_script_dir <- function() {
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg) > 0) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[1]), mustWork = FALSE)))
  }
  script_path <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
  if (!is.null(script_path) && nzchar(script_path)) {
    return(dirname(normalizePath(script_path, mustWork = FALSE)))
  }
  getwd()
}

# 设置工作目录
script_dir <- get_script_dir()

# 尝试找到scLncR安装目录
find_scLncR_home <- function() {
  # 尝试几种方法
  possibilities <- c(
    dirname(dirname(Sys.which("scLncR"))),
    Sys.getenv("SC_LNCR_HOME"),
    normalizePath(file.path(script_dir, "..")),
    getwd()
  )
  
  for (path in possibilities) {
    if (path != "" && dir.exists(path)) {
      # 检查是否是scLncR目录
      expected_files <- c("scLncR", "R", "shiny_app")
      found <- sapply(file.path(path, expected_files), 
                     function(x) file.exists(x) || dir.exists(x))
      
      if (sum(found) >= 2) {
        return(normalizePath(path))
      }
    }
  }
  
  return(NULL)
}

load_shiny_dependencies <- function(scLncR_home) {
  loader <- file.path(scLncR_home, "R", "modules", "scLncR_load_package.R")
  if (!file.exists(loader)) {
    stop("Package loader not found: ", loader)
  }
  source(loader, local = TRUE)
  scLncR_load_shiny_packages()
}

# 主函数
main <- function() {
  # 查找scLncR
  scLncR_home <- find_scLncR_home()
  
  if (is.null(scLncR_home)) {
    cat("ERROR: Could not find scLncR installation.\n")
    cat("Please make sure scLncR is installed and in your PATH.\n")
    cat("You can also set the SC_LNCR_HOME environment variable.\n")
    return(invisible(FALSE))
  }
  
  cat(sprintf("Found scLncR at: %s\n", scLncR_home))
  
  # 设置环境变量
  Sys.setenv(SC_LNCR_HOME = scLncR_home)
  load_shiny_dependencies(scLncR_home)
  
  # 切换到shiny_app目录
  shiny_dir <- file.path(scLncR_home, "shiny_app")
  
  if (!dir.exists(shiny_dir)) {
    cat(sprintf("ERROR: Shiny app directory not found: %s\n", shiny_dir))
    return(invisible(FALSE))
  }
  
  setwd(shiny_dir)
  
  cat("Loading required packages...\n")
  
  cat("\n")
  cat("Server starting...\n")
  cat("Open your browser and go to: http://localhost:3838\n")
  cat("Press Ctrl+C to stop the server\n")
  cat("========================================\n\n")
  
  # 运行应用
  shiny::runApp(
    appDir = ".",
    host = "0.0.0.0",
    port = 3838,
    launch.browser = FALSE,
    display.mode = "normal"
  )
}

# 执行
if (interactive()) {
  main()
} else {
  tryCatch({
    main()
  }, interrupt = function(e) {
    cat("\n\nShiny server stopped.\n")
    quit(save = "no", status = 0)
  })
}
