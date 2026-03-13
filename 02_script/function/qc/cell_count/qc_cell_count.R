#!/usr/bin/env Rscript

# ============================================================
# QC Cell Count - scRNA-seq
# ============================================================

#####################
#--- DESCRIPTION ---#
#####################

############################################################
# QC Cell Count Pipeline - scRNA-seq
#
# Author : Jordan DUTEL
# Contact : jordan.dutel@inserm.fr
#
# Description :
# This script performs quality control summary analysis on
# single-cell RNA-seq data stored as Seurat objects.
#
# The pipeline automatically generates:
# - Cell count statistics tables
# - Cell distribution barplots
# - QC UMAP visualizations
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
# 03_result/<dataset>/QC/Cell_Counts/<qc_state>/
#
# Execution example:
#
# log_time <- format(Sys.time(), "%Y%m%d_%H%M%S")
#
# callr::rscript(
#   script = "qc_cell_count.R",
#   cmdargs = c(
#     "-i", "object.qs",
#     "-d", "DatasetName",
#     "-o", "Results",
#     "-q", "After_QC"
#   ),
#   stdout = file.path("Logs", "qc", "cell_count", paste0("output_qc_cell_count_", log_time, ".log")),
#   stderr = file.path("Logs", "qc", "cell_count", paste0("error_qc_cell_count_", log_time, ".log"))
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
  library(tidyr)
  library(ggplot2)
  library(scales)
  library(htmlwidgets)
  library(DT)
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

# ------------------------------------------------------------
# CELL COUNT TABLE
# ------------------------------------------------------------

qc_cell_count_table <- function(seu, batch_col, condition_col) {
  
  message("Generating interactive table...")
  
  df <- seu[[]] %>%
    count(
      .data[[batch_col]],
      .data[[condition_col]],
      name = "Cells"
    ) %>%
    mutate(
      Percent_global = 
        percent(Cells / sum(Cells), accuracy = 0.1) 
      ) %>%
    group_by(.data[[condition_col]]) %>%
    mutate(
      Percent_within_condition =
        percent(Cells / sum(Cells), accuracy = 0.1)
    ) %>%
    ungroup() %>%
    rename(Batch = all_of(batch_col)) %>%
    pivot_wider(
      names_from = all_of(condition_col),
      values_from = Percent_within_condition,
      values_fill = "0%"
    ) %>% 
    mutate(
      # MM = ifelse(is.na(MM), "0%", MM),
      # HD = ifelse(is.na(HD), "0%", HD),
      Total = rowSums(select(., contains("Cells")))
    ) %>%
    bind_rows(
      summarise(
        .,
        across(where(is.numeric), sum),
        Batch = "Total"
      )
    )
  
  datatable(
    df,
    options = list(pageLength = 20, autoWidth = TRUE),
    rownames = FALSE
  )
}

# ------------------------------------------------------------
# BARPLOTS
# ------------------------------------------------------------

qc_barplot <- function(seu, colname, save_dir, dataset, qc_state, pal) {
  
  message("Generating barplots...")
  
  p <- ggplot(
    seu[[]],
    aes(x = .data[[colname]], fill = .data[[colname]])
  ) +
    geom_bar() +
    scale_fill_manual(values = pal, drop = FALSE) +
    theme_classic() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(hjust = 0.5, face = "bold")
    ) +
    ggtitle(paste("Number of cells -", colname)) +
    NoLegend()
  
  ggsave(
    filename = file.path(
      save_dir,
      tolower(paste0("barplot_", colname, "_", dataset, "_", qc_state, ".png"))
    ),
    plot = p,
    width = 6,
    height = 4,
    dpi=300
  )
  
  ggsave(
    filename = file.path(
      save_dir,
      tolower(paste0("barplot_", colname, "_", dataset, "_", qc_state, ".pdf"))
    ),
    plot = p,
    width = 6,
    height = 4
  )
  
  message("Barplots completed successfully.")
  
}

# ------------------------------------------------------------
# UMAP FACET
# ------------------------------------------------------------

qc_umap <- function(seu, group_col, save_dir, dataset, qc_state, pal) {
  
  
  if (!"umap" %in% names(seu@reductions)) {
    message("No UMAP reduction found. Skipping UMAP plot.")
    return(NULL)
  }
  
  message("Generating UMAP plot...")
  
  p <- DimPlot(seu, group.by = group_col) +
    facet_wrap(stats::as.formula(paste("~", group_col)), ncol = 6) +
    scale_color_manual(values = pal, drop = FALSE) +
    ggtitle(paste("UMAP by", group_col)) +
    NoLegend()
  
  ggsave(
    filename = file.path(
      save_dir,
      tolower(paste0("UMAP_", group_col, "_", dataset, "_", qc_state, ".png"))
    ),
    plot = p,
    width = 20,
    height = 10,
    dpi=300
  )
  
  ggsave(
    filename = file.path(
      save_dir,
      tolower(paste0("UMAP_", group_col, "_", dataset, "_", qc_state, ".pdf"))
    ),
    plot = p,
    width = 20,
    height = 10,
  )
  
  message("QC Cell Counts completed successfully.")
  
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
      c("-b", "--batch_col"),
      type = "character",
      default = "emulsion",
      help = "Batch column name [default = emulsion]"
    ),
    
    make_option(
      c("-c", "--condition_col"),
      type = "character",
      default = "condition",
      help = "Condition column name [default = condition]"
    ),
  
    make_option(
      c("-p", "--palette"),
      type = "character",
      default = pal,
      help = "Color palette to use"
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
  
  
  
  # ---- QC state validation ----
  valid_qc_states <- c("Before_QC", "before_qc", "After_QC", "after_qc")
  
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
  # 
  # # ---- Write permissions check ----
  # if (!file.access(final_outdir, 2) == 0) {
  #   stop("ERROR: No write permission on output directory.", call. = FALSE)
  # }
  # 
  # ---- Logging ----
  message("\n========== Pipeline Configuration ==========")
  message("Input file : ", opt$input)
  message("Dataset    : ", opt$dataset)
  message("QC state   : ", opt$qc_state)
  message("Output dir : ", opt$outdir)
  message("===========================================\n")
  
  # ------------------------------------------------------------
  # 3. LOAD SEURAT OBJECT
  # ------------------------------------------------------------
  
  message("Loading Seurat object...")
  seu <- qs::qread(opt$input)
  
  meta <- seu[[]]
  
  if (!opt$batch_col %in% colnames(meta)) {
    stop("ERROR: batch_col not found in metadata.", call. = FALSE)
  }
  
  if (!opt$condition_col %in% colnames(meta)) {
    stop("ERROR: condition_col not found in metadata.", call. = FALSE)
  }
  
  if (!inherits(seu, "Seurat")) {
    stop("ERROR: Input object is not a Seurat object.")
  }
  
  # ---- Create save dir ----
  save_dir <- file.path(opt$outdir, opt$dataset, "QC", "Cell_Counts", opt$qc_state)
  dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)
  message("Output directory created at: ", save_dir)
  
  # ------------------------------------------------------------
  # 4. CELL COUNT TABLE
  # ------------------------------------------------------------

  table_dt <- qc_cell_count_table(
    seu,
    opt$batch_col,
    opt$condition_col
  )
  
  saveWidget(
    table_dt,
    file = file.path(
      save_dir,
      tolower(paste0("Cell_Count_Table_", opt$dataset, "_", opt$qc_state, ".html"))
    ),
    selfcontained = TRUE
  )
  
  # ------------------------------------------------------------
  # 5. BARPLOTS
  # ------------------------------------------------------------
  
  qc_barplot(seu, opt$batch_col, save_dir, opt$dataset, opt$qc_state, pal)
  qc_barplot(seu, opt$condition_col, save_dir, opt$dataset, opt$qc_state, pal)
  
  # ------------------------------------------------------------
  # 6. UMAP FACET
  # ------------------------------------------------------------

  qc_umap(seu, opt$batch_col, save_dir, opt$dataset, opt$qc_state, pal)

  
  
  message("Cellcount QC completed successfully")
  
}

################
#--- SOURCE ---#
################

if (sys.nframe() == 0) {
  main()
}
