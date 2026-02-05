# shiny_app/launch.R
#!/usr/bin/env Rscript

cat("
========================================
     scLncR Shiny Application
========================================
")

# 设置工作目录
setwd(dirname(sys.frame(1)$ofile))

# 检查并安装必要的包
required_packages <- c(
  "shiny", "shinydashboard", "shinyjs", "shinyWidgets",
  "shinyFiles", "yaml", "DT", "fs", "processx"
)

cat("Checking required packages...\n")
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    cat(paste("Installing", pkg, "...\n"))
    install.packages(pkg, repos = "https://cloud.r-project.org")
    library(pkg, character.only = TRUE)
  }
}

# 设置最大上传大小
options(shiny.maxRequestSize = 500*1024^2)  # 500MB

# 设置文件选择器根目录
roots <- c(
  Home = normalizePath("~"),
  Current = normalizePath("."),
  Root = "/"
)

# 保存roots到全局环境
assign("roots", roots, envir = .GlobalEnv)

cat("\nStarting Shiny application...\n")
cat("• Server: http://localhost:3838\n")
cat("• Press Ctrl+C to stop the server\n")
cat("• Logs will be saved to ./logs/\n\n")

# 运行应用
shiny::runApp(
  appDir = ".",
  host = "0.0.0.0",
  port = 3838,
  launch.browser = FALSE
)