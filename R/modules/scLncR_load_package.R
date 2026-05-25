# =============================================================================
# scLncR - centralized R package loading utilities
# =============================================================================

#' Check whether an R package is available
#'
#' @description
#' Checks package availability without attaching it. Missing required packages
#' stop with a clear message; missing optional packages produce a warning.
#'
#' @param pkg Package name.
#' @param required Logical. Whether a missing package should stop execution.
#' @return Invisibly returns TRUE when the package is available, otherwise FALSE.
#' @export
scLncR_check_package <- function(pkg, required = TRUE) {
  ok <- requireNamespace(pkg, quietly = TRUE)
  if (!ok) {
    msg <- sprintf(
      "Required R package '%s' is not installed or not available in this environment.",
      pkg
    )
    if (isTRUE(required)) {
      stop(msg, call. = FALSE)
    } else {
      warning(msg, call. = FALSE)
    }
  }
  invisible(ok)
}

#' Load R packages required by scLncR modules
#'
#' @description
#' Provides centralized package loading utilities for scLncR. Module-specific
#' functions load only the packages required by the corresponding workflow.
#'
#' @param pkgs A character vector of package names.
#' @param required Logical. Whether missing packages should stop execution.
#' @return Invisibly returns TRUE when loading succeeds.
#' @export
scLncR_load_packages <- function(pkgs, required = TRUE) {
  pkgs <- unique(as.character(pkgs))
  pkgs <- pkgs[nzchar(pkgs)]
  for (pkg in pkgs) {
    ok <- scLncR_check_package(pkg, required = required)
    if (isTRUE(ok)) {
      suppressPackageStartupMessages(
        suppressWarnings(
          library(pkg, character.only = TRUE)
        )
      )
    }
  }
  invisible(TRUE)
}

#' Load core packages for CLI parsing and YAML configs
#' @return Invisibly returns TRUE.
#' @export
scLncR_load_core_packages <- function() {
  scLncR_load_packages(c("optparse", "yaml"), required = TRUE)
}

#' Load packages for raw FASTQ QC
#' @return Invisibly returns TRUE.
#' @export
scLncR_load_qc_packages <- function() {
  scLncR_load_packages("optparse", required = TRUE)
}

#' Load packages for prelnc candidate discovery
#' @return Invisibly returns TRUE.
#' @export
scLncR_load_prelnc_packages <- function() {
  scLncR_load_packages(c("optparse", "yaml"), required = TRUE)
}

#' Load packages for count workflows
#' @return Invisibly returns TRUE.
#' @export
scLncR_load_count_packages <- function() {
  scLncR_load_packages(c("optparse", "yaml"), required = TRUE)
}

#' Load packages for dataProcess workflows
#' @return Invisibly returns TRUE.
#' @export
scLncR_load_dataProcess_packages <- function(annotation_method = NULL) {
  pkgs <- c(
    "optparse", "stringr", "gridExtra", "psych", "pheatmap",
    "Seurat", "dplyr", "ggsci", "tidyverse", "patchwork",
    "ggplot2", "reshape2", "ggrepel"
  )
  annotation_method <- annotation_method %||% ""
  if (identical(annotation_method, "SingleR")) pkgs <- c(pkgs, "SingleR")
  if (identical(annotation_method, "scMM")) pkgs <- c(pkgs, "scMayoMap")
  scLncR_load_packages(pkgs, required = TRUE)
}

#' Load packages for downstream function analysis
#' @return Invisibly returns TRUE.
#' @export
scLncR_load_function_packages <- function(modules = character(0)) {
  pkgs <- c(
    "optparse", "Seurat", "dplyr", "tibble", "ggplot2",
    "reshape2", "patchwork", "scales", "Biobase", "VGAM"
  )
  modules <- unique(tolower(as.character(modules)))
  if ("monocle2" %in% modules) pkgs <- c(pkgs, "monocle")
  if ("wgcna" %in% modules) pkgs <- c(pkgs, "hdWGCNA", "WGCNA", "UCell")
  scLncR_load_packages(pkgs, required = TRUE)
}

#' Load packages for normalization benchmark
#' @return Invisibly returns TRUE.
#' @export
scLncR_load_normalization_benchmark_packages <- function(strategies = character(0)) {
  pkgs <- c("Seurat", "ggplot2", "dplyr", "reshape2", "yaml")
  strategies <- unique(tolower(as.character(strategies)))
  if ("scran" %in% strategies) {
    pkgs <- c(pkgs, "SingleCellExperiment", "scran", "scuttle", "SummarizedExperiment")
  }
  scLncR_load_packages(pkgs, required = TRUE)
}

#' Load packages for the scLncR Shiny GUI
#'
#' @description
#' Loads packages required by the Shiny graphical interface. Missing packages
#' stop with a clear message; scLncR does not install packages at runtime.
#'
#' @return Invisibly returns TRUE.
#' @export
scLncR_load_shiny_packages <- function() {
  scLncR_load_packages(
    c(
      "shiny", "shinydashboard", "shinyjs", "yaml",
      "shinyWidgets", "shinyFiles", "DT", "fs"
    ),
    required = TRUE
  )
}

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0) return(y)
  x
}
