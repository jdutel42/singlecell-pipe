#!/usr/bin/env Rscript

# ============================================================
# QC Filtering - scRNA-seq
# ============================================================



#####################
#--- DESCRIPTION ---#
#####################

############################################################
# QC Filtering Pipeline - scRNA-seq
#
# Author : Jordan DUTEL
# Contact : jordan.dutel@inserm.fr
#
# Description :
# This script performs filtering on
# single-cell RNA-seq data stored as Seurat objects.
#
# The pipeline automatically:
#   1. Loads a Seurat object (.qs)
#   2. Applies QC rules based on metadata columns
#   3. Creates a new QC status column (pass/fail)
#   4. Optionally subsets passing cells
#   5. Saves the updated object
#
# Execution is designed for pipeline automation using
# callr::rscript with separate logging streams for:
# - Standard output logs
# - Error logs
#
# Inputs:
# - Seurat object stored in .qs format
# - Dataset name
#
# Outputs:
# Filtered object  are stored in:
# 02_data/Processed_data/<object>.qs
#
# Execution example:
#
# log_time <- format(Sys.time(), "%Y%m%d_%H%M%S")
#
# callr::rscript(
#   script = "qc_filter.R",
#   cmdargs = c(
#     "-i", "object.qs",
#     "-d", "DatasetName",
#     "-o", "Results_dir",
#     "-q", "After_QC"
#     "-s", "TRUE"
#   ),
#   stdout = file.path("06_log", "qc", "filter", paste0("output_qc_filter_", log_time, ".log")),
#   stderr = file.path("06_log", "qc", "filter", paste0("error_qc_filter_", log_time, ".log"))
# )
#
############################################################

#############
#--- LIB ---#
#############

suppressPackageStartupMessages({
  library(optparse)
  library(qs)
  library(Seurat)
  library(ggplot2)
  library(qs)
})


##################
#--- FUNCTION ---#
##################

run_qc_filter <- function(seu_file_path, dataset_name, qc_name, filter_rule, subset_seu = TRUE, save_seu = TRUE) {
  
  # ---- Load subset Seurat ----
  message("Loading Seurat object...")
  
  seu <- qs::qread(seu_file_path)
  
  if (!inherits(seu, "Seurat")) {
    stop("ERROR: Input object is not a Seurat object.")
  }
  
  meta <- seu@meta.data
  
  # --------------------------------------------------
  # Evaluate rule safely
  # --------------------------------------------------
  
  if (is.null(filter_rule)) {
    stop("Filter rule cannot be NULL")
  }
  
  # Make metadata variables available in evaluation environment
  eval_env <- as.list(meta)
  
  # Evaluate logical rule
  condition <- tryCatch(
    eval(str2lang(filter_rule), envir = eval_env),
    error = function(e) stop("Error evaluating filter rule: ", e$message)
  )
  

  
  if (!is.logical(condition)) {
    stop("Filter rule must evaluate to logical values")
  }
  
  # --------------------------------------------------
  # QC annotation
  # In scRNA-seq pipelines:
  # Usually PASS = FALSE condition
  # Here we assume rule = cells to REMOVE
  # --------------------------------------------------
  
  seu[[qc_name]] <- ifelse(condition, "fail", "pass")
  
  pass_n <- sum(seu@meta.data[[qc_name]] == "pass")
  fail_n <- sum(seu@meta.data[[qc_name]] == "fail")
  
  message("QC filtering applied:")
  message("  PASS: ", pass_n)
  message("  FAIL: ", fail_n)
  message("  Total: ", pass_n + fail_n)

  # --------------------------------------------------
  # Optional subset (keep pass cells)
  # --------------------------------------------------
  if (subset_seu) {
    message("Subsetting Seurat object (only keeping 'pass' cells)")
    
    keep_cells <- rownames(seu@meta.data)[seu@meta.data[[qc_name]] == "pass"]
    
    seu <- subset(seu, cells = keep_cells)
  }
  
  # --------------------------------------------------
  # Save
  # --------------------------------------------------
  
  if (save_seu) {
    
    message("Save updated seurat object...")
    
    qs::qsave(seu, seu_file_path)
    
    cat("Filtered seurat object saved here :", seu_file_path, "\n")
    
    message("Filtering done successfully.")
  }

  return(seu)
}


##############
#--- MAIN ---#
##############

main <- function() {
  
  # ------------------------------------------------------------
  # 1. CLI ARGUMENTS
  # ------------------------------------------------------------
  
  option_list <- list(
    
    make_option(
      c("-i", "--input"),
      type = "character",
      help = "Path to Seurat object (.qs)"
    ),
    
    make_option(
      c("-d", "--dataset"),
      type = "character",
      help = "Dataset name"
    ),
    
    make_option(
      c("-n", "--qc_name"),
      type = "character",
      help = "QC column name to create"
    ),
    
    make_option(
      c("-r", "--rule"),
      type = "character",
      default = NULL,
      help = "QC filtering rule (e.g. '(mito_status == \"high\" & ribo_status == \"high\") | lib_size == \"low\"')"
    ),
    
    make_option(
      c("-b", "--subset"),
      type = "logical",
      default = TRUE,
      help = "Subset Seurat object to only keep 'PASS' cells (TRUE/FALSE) [default = TRUE]"
    ),
    
    make_option(
      c("-s", "--save"),
      type = "logical",
      default = TRUE,
      help = "Save Seurat object"
    )
  )
  
  parser <- OptionParser(
    usage = "%prog -i object.qs -d dataset_name [options]",
    option_list = option_list
  )
  
  opt <- parse_args(parser)
  
  # ------------------------------------------------------------
  # 2. VALIDATION
  # ------------------------------------------------------------
  
  # ---- Required arguments ----
  if (is.null(opt$input) || is.null(opt$dataset)) {
    print_help(parser)
    stop("ERROR: --input and --dataset are required.", call. = FALSE)
  }
  
  # ---- Input file checks ----
  if (!file.exists(opt$input)) {
    stop("ERROR: Input file does not exist: ", opt$input, call. = FALSE)
  }
  
  if (file.info(opt$input)$size == 0) {
    stop("ERROR: Input file is empty.", call. = FALSE)
  }
  
  
  # ---- Logging ----
  message("\n========== Pipeline Configuration ==========")
  message("Input file : ", opt$input)
  message("Dataset    : ", opt$dataset)
  message("QC name  : ", opt$qc_name)
  message("Filter rule : ", opt$rule)
  message("Subset to PASS cells : ", opt$subset)
  message("Save Seurat : ", opt$save)
  message("===========================================\n")
  
  # ---- Run QC Mitochondria ----
  run_qc_filter(opt$input, opt$dataset, opt$qc_name, opt$rule, opt$subset, opt$save)
  
}

################
#--- SOURCE ---#
################

if (sys.nframe() == 0) {
  main()
}