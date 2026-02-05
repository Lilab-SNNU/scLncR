# =============================================================================
# scLncR - Smart Config Loader for LncExplore (Dynamic Module Validation)
# Author: [Yin SW]
# Description: Validates ONLY parameters for modules specified in run_modules.
# Dependencies: yaml
# =============================================================================


load_LncExplore_config <- function(config_path) {
  if (!requireNamespace("yaml", quietly = TRUE)) {
    stop("Package 'yaml' is required. Install with: install.packages('yaml')")
  }
  if (!file.exists(config_path)) stop("Config file not found: ", config_path)
  
  cfg <- yaml::read_yaml(config_path)
  
  # -----------------------------
  # 1. Required global parameter: input_seurat (depends on all modules)
  # -----------------------------
  if (is.null(cfg$input_seurat) || !file.exists(cfg$input_seurat)) {
    stop("Global parameter 'input_seurat' is required and must be a valid path.")
  }
  if (!grepl("\\.rds$", cfg$input_seurat, ignore.case = TRUE)) {
    warning("Input file lacks .rds extension. Verify it is a Seurat object.")
  }

  # -----------------------------
  # 2. Determine the modules to run (flexible configuration supported)
  # -----------------------------
  ALL_MODULES <- c("location", "monocle2", "wgcna")
  
  if (is.null(cfg$run_modules) || length(cfg$run_modules) == 0) {
    message("⚠️  No 'run_modules' specified. Defaulting to ALL modules: ", 
            paste(ALL_MODULES, collapse = ", "))
    run_modules <- ALL_MODULES
  } else {
    # Standardization: Convert to lowercase + Deduplication
    run_modules <- unique(tolower(trimws(unlist(strsplit(
      paste(cfg$run_modules, collapse = ","), ","
    )))))
    
    # Validate module name
    invalid <- setdiff(run_modules, ALL_MODULES)
    if (length(invalid) > 0) {
      stop("Invalid module(s) in 'run_modules': ", paste(invalid, collapse = ", "), 
           ". Valid options: ", paste(ALL_MODULES, collapse = ", "))
    }
    message("✅ Running selected modules: ", paste(run_modules, collapse = ", "))
  }
  
  # -----------------------------
  # 3. Validate each module as needed (only validate modules in run_modules).
  # -----------------------------
  # The validation logic is encapsulated as an internal function (to keep the main flow clear).
  validate_location <- function(loc_cfg) {
    req <- c("LOG2FC_THRESH", "PADJ_THRESH", "output_path", "lncRNA_name")
    if (any(!req %in% names(loc_cfg))) stop("Missing required 'location' parameters")
    if (loc_cfg$LOG2FC_THRESH <= 0) stop("'location$LOG2FC_THRESH' must be > 0")
    if (loc_cfg$PADJ_THRESH <= 0 || loc_cfg$PADJ_THRESH > 1) stop("'location$PADJ_THRESH' must be in (0,1]")
    if (grepl("[^a-zA-Z0-9-]", loc_cfg$lncRNA_name)) stop("Invalid chars in lncRNA_name")
    if (!dir.exists(dirname(loc_cfg$output_path))) dir.create(dirname(loc_cfg$output_path), recursive = TRUE)
  }
  
  validate_monocle2 <- function(mono_cfg) {
    req <- c("qval", "reduceDimensionMethod", "output_path", "hub_genes")
    if (any(!req %in% names(mono_cfg))) stop("Missing required 'monocle2' parameters")
    if (mono_cfg$qval <= 0 || mono_cfg$qval > 1) stop("'monocle2$qval' must be in (0,1]")
    if (!mono_cfg$reduceDimensionMethod %in% c("DDRTree", "ICA", "tSNE")) 
      stop("Invalid 'monocle2$reduceDimensionMethod'")
    if (!dir.exists(dirname(mono_cfg$output_path))) dir.create(dirname(mono_cfg$output_path), recursive = TRUE)
  }
  
  validate_wgcna <- function(wgcna_cfg) {
    req <- c("output_path", "pro_name", "gene_select_method")
    if (any(!req %in% names(wgcna_cfg))) stop("Missing required 'wgcna' parameters")
    if (!wgcna_cfg$gene_select_method %in% c("fraction", "variance")) 
      stop("Invalid 'wgcna$gene_select_method'")
    if (!dir.exists(dirname(wgcna_cfg$output_path))) dir.create(dirname(wgcna_cfg$output_path), recursive = TRUE)
    # cell_types and lnc_name can be empty; only the type is checked.
    if (!is.null(wgcna_cfg$cell_types) && !is.character(unlist(wgcna_cfg$cell_types))) 
      stop("'wgcna$cell_types' must be character vector or list")
  }
  
  # Perform dynamic verification
  if ("location" %in% run_modules) {
    if (is.null(cfg$location)) stop("Module 'location' selected but 'location' section missing in config")
    validate_location(cfg$location)
    message("✓ Validated: Location Analysis parameters")
  }
  
  if ("monocle2" %in% run_modules) {
    if (is.null(cfg$monocle2)) stop("Module 'monocle2' selected but 'monocle2' section missing in config")
    validate_monocle2(cfg$monocle2)
    message("✓ Validated: Monocle2 Trajectory parameters")
  }
  
  if ("wgcna" %in% run_modules) {
    if (is.null(cfg$wgcna)) stop("Module 'wgcna' selected but 'wgcna' section missing in config")
    validate_wgcna(cfg$wgcna)
    message("✓ Validated: WGCNA parameters")
  }
  
  # -----------------------------
  # 4. Inject the runtime plan into the configuration (for use by run_LncExplore.R)
  # -----------------------------
  cfg$.run_modules <- run_modules  # Private field, marking the actual module to be run.
  
  message("\n✨ Configuration validated successfully for modules: ", 
          paste(run_modules, collapse = ", "))
  message("   Input Seurat: ", cfg$input_seurat)
  return(cfg)
}