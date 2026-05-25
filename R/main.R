# R/main.R
# Main dispatcher for scLncR CLI

# Parse internal --script-dir argument
script_dir <- Sys.getenv("SC_LNCR_HOME", unset = NA)
if (is.na(script_dir) || !dir.exists(script_dir)) {
  stop("Internal error: SC_LNCR_HOME not set or invalid.")
}

# Parse user-facing arguments
user_args <- commandArgs(trailingOnly = TRUE)

# -----------------------------
# Global help: scLncR -h or scLncR --help or no args
# -----------------------------
if (length(user_args) == 0 || 
    user_args[1] %in% c("-h", "--help")) {
  
  cat("scLncR v0.1.0 - Single-cell lncRNA Discovery Pipeline\n")
  cat("Usage: scLncR <command> [options]\n\n")
  
  cat("Available commands:\n")
  cat("  qc            Raw FASTQ quality control (FastQC + MultiQC; no trimming)\n")
  cat("  prelnc        Candidate lncRNA discovery (raw FASTQ-first, technology-aware)\n")
  cat("  count         lncRNA-aware quantification (raw FASTQ-first interface)\n")
  cat("  dataProcess   ScRNA-seq expression count preprocess and annotation \n")
  cat("  function      DownStream analysis to explore lncRNA function \n")
  cat("  shiny         Launch the scLncR Shiny GUI \n\n")
  
  cat("Examples:\n")
  cat("  scLncR prelnc -c config.yaml\n")
  cat("  scLncR prelnc --help\n\n")
  
  cat("Options:\n")
  cat("  -h, --help    Show this global help message\n\n")
  
  cat("Note: Each command has its own --help. For example:\n")
  cat("      scLncR prelnc --help\n")
  
  q()
}

# -----------------------------
# Dispatch to subcommands
# -----------------------------
cmd <- user_args[1]
sub_args <- user_args[-1]
# Source modules only when needed (or source all if few)
source(file.path(script_dir, "R", "modules", "scLncR_load_package.R"))

if (cmd == "qc") {
  scLncR_load_qc_packages()
  source(file.path(script_dir, "R", "utils", "config_loader_qc.R"))
  source(file.path(script_dir, "R", "modules", "scLncR_QC.R"))
  run_qc(sub_args, script_dir)
} else if (cmd == "prelnc") {
  scLncR_load_prelnc_packages()
  source(file.path(script_dir, "R", "utils", "config_loader_prelnc.R"))
  source(file.path(script_dir, "R", "modules", "scLncR_LncPre.R"))
  run_prelnc(sub_args, script_dir)
} else if (cmd == "count") {
  scLncR_load_count_packages()
  source(file.path(script_dir, "R", "utils", "config_loader_count.R"))
  source(file.path(script_dir, "R", "modules", "scLncR_Count.R"))
  run_count(sub_args, script_dir)
} else if (cmd == "dataProcess") {
  source(file.path(script_dir, "R", "utils", "config_loader_dataProcess.R"))
  source(file.path(script_dir, "R", "modules", "scLncR_dataProcess.R"))
  run_dataProcess(sub_args, script_dir)
} else if (cmd == "function") {
  source(file.path(script_dir, "R", "utils", "config_loader_function.R"))
  source(file.path(script_dir, "R", "modules", "scLncR_dataProcess.R"))
  source(file.path(script_dir, "R", "modules", "scLncR_function.R"))
  run_function(sub_args, script_dir)
} else if (cmd == "shiny") {
  if (length(sub_args) > 0 && any(sub_args %in% c("-h", "--help"))) {
    cat("scLncR shiny - launch the scLncR Shiny GUI\n")
    cat("Usage: scLncR shiny\n\n")
    cat("The GUI starts on host 0.0.0.0 and port 3838.\n")
    q(status = 0)
  }
  scLncR_load_shiny_packages()
  shiny::runApp(
    appDir = file.path(script_dir, "shiny_app"),
    host = "0.0.0.0",
    port = 3838,
    launch.browser = FALSE
  )
} else {
  cat("Error: Unknown command '", cmd, "'.\n", sep = "")
  cat("Run 'scLncR -h' for available commands.\n")
  q(status = 1)
}
