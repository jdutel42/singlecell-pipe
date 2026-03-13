#!/usr/bin/env Rscript

# ============================================================
# QC Library Diversity (nFeatureRNA) - scRNA-seq
# ============================================================

#####################
#--- DESCRIPTION ---#
#####################

############################################################
# QC Library Diversity Pipeline - scRNA-seq
#
# Author : Jordan DUTEL
# Contact : jordan.dutel@inserm.fr
#
# Description :
# This script performs quality control summary analysis on
# single-cell RNA-seq data stored as Seurat objects.
#
# The pipeline automatically generates:
# - Detect high and low library diversity content outliers (> 3 MADs)
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
# QC reports are stored in:
# 03_result/<dataset>/QC/Library_Diversity/<qc_state>/
#
# Execution example:
#
# log_time <- format(Sys.time(), "%Y%m%d_%H%M%S")
#
# callr::rscript(
#   script = "qc_library_diversity.R",
#   cmdargs = c(
#     "-i", "object.qs",
#     "-d", "DatasetName",
#     "-o", "Results_dir",
#     "-q", "After_QC"
#     "-s", "TRUE"
#   ),
#   stdout = file.path("06_log", "qc", "library_diversity", paste0("output_qc_library_diversity_", log_time, ".log")),
#   stderr = file.path("06_log", "qc", "library_diversity", paste0("error_qc_library_diversity_", log_time, ".log"))
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
  library(dplyr)
  library(ggplot2)
  library(ggbeeswarm)
  library(scater)   # for isOutlier
})


################
#--- COLOR ----#
################

pal <- c(
  "#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00",
  "#A65628","#F781BF","#999999","#66C2A5","#FC8D62",
  "#8DA0CB","#E78AC3","#A6D854","#FFD92F","#E5C494",
  "#1B9E77","#D95F02","#7570B3","#E7298A","#66A61E",
  "#E6AB02","#A6761D","#666666","#8DD3C7","#BEBADA",
  "#FB8072","#80B1D3","#FDB462","#B3DE69","#FCCDE5",
  "#D9D9D9","#BC80BD","#CCEBC5","#FFED6F","#1F78B4",
  "#33A02C","#FB9A99","#E31A1C","#FDBF6F","#CAB2D6",
  "#6A3D9A","#FFFF99","#B15928","#8E0152","#DE77AE",
  "#7FBC41","#4EB265","#4393C3","#2166AC","#B2182B",
  "#F4A582","#92C5DE","#FDDBC7","#D1E5F0","#542788"
)


##################
#--- FUNCTION ---#
##################

run_qc_libdiv <- function(seu_file_path, dataset_name, batch_col, outdir, qc_state, save_seu = TRUE) {
  
  # ---- Load subset Seurat ----
  message("Loading Seurat object...")
  
  seu <- qs::qread(seu_file_path)
  
  if (!inherits(seu, "Seurat")) {
    stop("ERROR: Input object is not a Seurat object.")
  }
  
  meta <- seu[[]]
  
  if (!batch_col %in% colnames(meta)) {
    stop("ERROR: batch_col not found in metadata.", call. = FALSE)
  }
  
  # ---- Perform Library Diversity QC using isOutlier ----
  qc_lib_diversity <- meta %>%
    group_by(.data[[batch_col]]) %>%
    mutate(
      low_div  = as.logical(scater::isOutlier(nFeature_RNA, nmads = 3, type = "lower")),
      high_div = as.logical(scater::isOutlier(nFeature_RNA, nmads = 3, type = "higher")),
      libdiv_status = case_when(
        low_div  ~ "low",
        high_div ~ "high",
        TRUE     ~ "normal"
      )
    ) %>%
    ungroup() %>%
    filter(nFeature_RNA > 0)
  
  # Add diversity status to Seurat metadata
  seu <- AddMetaData(seu, qc_lib_diversity$libdiv_status, col.name = "libdiv_status")
  message("Library diversity outliers calculated and annotated in metadata.")
  
  # ---- Create save dir ----
  save_dir <- file.path(outdir, dataset_name, "QC", "Library_Diversity", "Before_QC")
  dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)
  message("Output directory created at: ", save_dir)
  
  # ---- Density plot ----
  message("Generating density plot...")
  p1 <- ggplot(qc_lib_diversity, aes(x = nFeature_RNA)) +
    geom_density(fill = "grey80", alpha = 0.6) +
    geom_rug(aes(color = libdiv_status), alpha = 0.6, sides = "b") +
    scale_x_log10() +
    scale_color_manual(values = c(low = "blue", normal = "green", high = "red")) +
    facet_wrap(~.data[[batch_col]]) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    labs(color = "Library diversity QC",
         x = "nFeatureRNA (log10)",
         y = "Density")
  
  ggsave(file.path(save_dir, tolower(paste0("densityplot_", qc_state, "_library_diversity_by_batch_qc_", dataset_name, ".png"))),
         plot = p1, 
         width = 12, 
         height = 6,
         dpi = 300)
  
  ggsave(file.path(save_dir, tolower(paste0("densityplot_", qc_state, "_library_diversity_by_batch_qc_", dataset_name, ".pdf"))),
         plot = p1, 
         width = 12, 
         height = 6)
  message("Density plot saved.")
  
  # ---- Violin + quasirandom plot ----
  message("Generating violin + quasirandom plot...")
  p2 <- ggplot(qc_lib_diversity,
               aes(x = .data[[batch_col]], y = nFeature_RNA, color = libdiv_status)) +
    geom_violin(fill = "grey90", color = NA) +
    geom_quasirandom(alpha = 0.5, size = 0.8) +
    scale_y_log10() +
    scale_color_manual(values = c(low = "blue", normal = "green", high = "red")) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
  
  ggsave(file.path(save_dir, tolower(paste0("violin_quasirandom_", qc_state, "_library_diversity_by_batch_qc_", dataset_name, ".png"))),
         plot = p2, 
         width = 12, 
         height = 6,
         dpi = 300)
  
  ggsave(file.path(save_dir, tolower(paste0("violin_quasirandom_", qc_state, "_library_diversity_by_batch_qc_", dataset_name, ".pdf"))),
         plot = p2, 
         width = 12, 
         height = 6)
  message("Violin + quasirandom plot saved.")
  
  # ---- Save annotated Seurat object ----
  message("Save updated seurat object...")
  if (save_seu) {
    qs::qsave(seu, file = seu_file_path)
    cat("Annotated seurat object save here :", seu_file_path, "\n")
  }
  message("Seurat object with library diversity QC annotations saved.")
  
  message("Library diversity QC completed successfully")
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
      c("-b", "--batch_col"),
      type = "character",
      help = "Batch column name"
    ),
    
    make_option(
      c("-o", "--outdir"),
      type = "character",
      default = "03_result",
      help = "Base output directory [default = 03_result]"
    ),
    
    make_option(
      c("-q", "--qc_state"),
      type = "character",
      default = NULL,
      help = "QC state: Before_QC or After_QC"
    ),
    
    make_option(
      c("-s", "--save"),
      type = "logical",
      default = TRUE,
      help = "TRUE or FALSE to save annotated Seurat object [default = TRUE]"
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
  
  # # ---- Write permissions check ----
  # if (!file.access(final_outdir, 2) == 0) {
  #   stop("ERROR: No write permission on output directory.", call. = FALSE)
  # }
  
  # ---- QC state validation ----
  valid_qc_states <- c("Before_QC", "After_QC")
  
  if (!is.null(opt$qc_state)) {
    if (!opt$qc_state %in% valid_qc_states) {
      stop(
        "ERROR: qc_state must be one of: ",
        paste(valid_qc_states, collapse = " | "),
        call. = FALSE
      )
    }
  }
  
  
  
  
  # # ---- colors ----
  # if (is.null(opt$palette)) {
  #   pal <- colors_55
  # } else {
  #   pal <- strsplit(opt$palette, ",")[[1]]
  # }
  #
  # # ---- Auto detect QC state ----
  # base_path <- file.path(
  #   opt$outdir,
  #   opt$dataset,
  #   "QC",
  #   "Cell_Counts"
  # )
  # 
  # if (is.null(opt$qc_state)) {
  #   
  #   if (dir.exists(file.path(base_path, "Before_QC"))) {
  #     qc_state <- "After_QC"
  #   } else {
  #     qc_state <- "Before_QC"
  #   }
  #   
  # } else {
  #   qc_state <- opt$qc_state
  # }
  # 
  # # ---- Output directory validation ----
  # final_outdir <- file.path(
  #   base_path,
  #   qc_state
  # )
  # 
  # dir.create(final_outdir, recursive = TRUE, showWarnings = FALSE)
  
  
  
  # ---- Logging ----
  message("\n========== Pipeline Configuration ==========")
  message("Input file : ", opt$input)
  message("Dataset    : ", opt$dataset)
  message("Batch col  : ", opt$batch_col)
  message("Output dir : ", opt$outdir)
  message("QC state   : ", opt$qc_state)
  message("Save Seurat : ", opt$save)
  message("===========================================\n")
  
  # ---- Run QC Mitochondria ----
  run_qc_libdiv(opt$input, opt$dataset, opt$batch_col, opt$outdir, opt$qc_state, opt$save)
  
}

################
#--- SOURCE ---#
################

if (sys.nframe() == 0) {
  main()
}