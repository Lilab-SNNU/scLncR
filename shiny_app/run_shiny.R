# shiny_app/run_shiny.R
# 简化的启动脚本

#!/usr/bin/env Rscript

cat("========================================\n")
cat("   scLncR Shiny GUI\n")
cat("========================================\n\n")

# 设置工作目录
script_dir <- dirname(sys.frame(1)$ofile)
if (length(script_dir) == 0) {
  script_dir <- getwd()
}

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
  
  # 切换到shiny_app目录
  shiny_dir <- file.path(scLncR_home, "shiny_app")
  
  if (!dir.exists(shiny_dir)) {
    cat(sprintf("ERROR: Shiny app directory not found: %s\n", shiny_dir))
    return(invisible(FALSE))
  }
  
  setwd(shiny_dir)
  
  cat("Loading required packages...\n")
  required_packages <- c("shiny", "shinydashboard", "shinyjs", "yaml", 
                         "shinyWidgets", "shinyFiles", "DT")
  
  for (pkg in required_packages) {
    if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
      cat(sprintf("Installing %s...\n", pkg))
      install.packages(pkg, repos = "https://cloud.r-project.org")
      library(pkg, character.only = TRUE)
    }
  }
  
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