# shiny_app/launch.R
# 这个文件专门处理路径问题并启动Shiny应用

#!/usr/bin/env Rscript

cat("========================================\n")
cat("   scLncR Shiny GUI Launcher\n")
cat("========================================\n\n")

# 检查并安装必要的包
required_packages <- c("shiny", "shinydashboard", "shinyjs", "yaml", 
                       "shinyWidgets", "shinyFiles", "shinyalert", "DT")

for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    cat(sprintf("Installing %s...\n", pkg))
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
}

# 获取scLncR安装目录
get_scLncR_home <- function() {
  # 方法1: 从系统PATH中查找
  scLncR_path <- Sys.which("scLncR")
  if (scLncR_path != "") {
    install_dir <- dirname(scLncR_path)
    cat(sprintf("Found scLncR at: %s\n", scLncR_path))
    cat(sprintf("Installation directory: %s\n", install_dir))
    return(install_dir)
  }
  
  # 方法2: 从环境变量获取
  if (Sys.getenv("SC_LNCR_HOME") != "") {
    install_dir <- Sys.getenv("SC_LNCR_HOME")
    cat(sprintf("Using SC_LNCR_HOME: %s\n", install_dir))
    return(install_dir)
  }
  
  # 方法3: 猜测当前目录
  if (file.exists("scLncR")) {
    cat("Found scLncR in current directory\n")
    return(normalizePath("."))
  }
  
  # 方法4: 用户输入
  cat("scLncR not found in PATH. Please enter the installation directory:\n")
  install_dir <- readline("scLncR installation directory: ")
  
  if (!dir.exists(install_dir)) {
    stop("Invalid directory. Please make sure scLncR is properly installed.")
  }
  
  return(normalizePath(install_dir))
}

# 主启动函数
launch_shiny <- function() {
  cat("Setting up environment...\n")
  
  # 获取scLncR安装目录
  scLncR_home <- get_scLncR_home()
  
  # 设置环境变量
  Sys.setenv(SC_LNCR_HOME = scLncR_home)
  
  # 切换到shiny_app目录
  shiny_dir <- file.path(scLncR_home, "shiny_app")
  if (!dir.exists(shiny_dir)) {
    stop(sprintf("Shiny app directory not found: %s", shiny_dir))
  }
  
  setwd(shiny_dir)
  
  cat(sprintf("Launching Shiny app from: %s\n", shiny_dir))
  cat("Server will be available at: http://localhost:3838\n")
  cat("Press Ctrl+C to stop the server\n\n")
  
  # 运行Shiny应用
  shiny::runApp(
    appDir = ".",
    host = "0.0.0.0",
    port = 3838,
    launch.browser = FALSE
  )
}

# 执行启动
if (interactive()) {
  launch_shiny()
} else {
  # 如果通过Rscript运行，捕获中断信号
  tryCatch({
    launch_shiny()
  }, interrupt = function(e) {
    cat("\n\nShiny server stopped by user.\n")
    quit(save = "no", status = 0)
  })
}