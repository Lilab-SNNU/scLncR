# shiny_app/run_app.R
# 启动Shiny应用的脚本

#!/usr/bin/env Rscript

# 设置工作目录为脚本所在目录
script_dir <- dirname(sys.frame(1)$ofile)
if (length(script_dir) == 0) {
  script_dir <- getwd()
}
setwd(script_dir)

# 检查并安装必要的包
required_packages <- c("shiny", "shinydashboard", "shinyjs", "yaml", 
                       "shinyBS", "shinyWidgets")

cat("Checking required packages...\n")
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    cat(paste("Installing", pkg, "...\n"))
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
}

# 运行Shiny应用
cat("\nStarting scLncR Shiny application...\n")
cat("Server will be available at: http://localhost:3838\n")
cat("Press Ctrl+C to stop the server\n\n")

shiny::runApp(host = "0.0.0.0", port = 3838)