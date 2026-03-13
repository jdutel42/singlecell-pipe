# ============================================================
#   QC Cell counts - scRNA-seq
# Seurat
# ============================================================

##############
#--- HELP ---#
##############

print_help <- function() {
  cat("
Usage:
Rscript qc_cell_count.R <seu_file> <dataset_name> [output_path] [qc_state]

Arguments:
seu_file : path to Seurat object (required)
dataset_name : dataset name (required)
output_path : base output directory (optional)
qc_state : 'Before_QC' or 'After_QC' (optional)

Behavior:

If output_path not provided:
./03_result/<dataset_name>/QC/Cell_Counts/<qc_state>/

If qc_state not provided:
If Before_QC exists → qc_state = After_QC
Else → qc_state = Before_QC

Examples:
Rscript qc_cell_count.R ./seu_CD4.rds CD4_Tcells
Rscript qc_cell_count.R ./seu_CD4.rds CD4_Tcells ./Results
Rscript qc_cell_count.R ./seu_CD4.rds CD4_Tcells ./Results After_QC
")
}

##############
#--- ARGS ---#
##############

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2 || "--help" %in% args) {
  print_help()
  quit(status = 0)
}

seu_file <- args[1] # required: path to Seurat object

dataset_name <- args[2] # required: dataset name (e.g., CD4_Tcells)

# optional arguments

outdir_arg <- if (length(args) >= 3) args[3] else NULL
qc_state_arg <- if (length(args) >= 4) args[4] else NULL

save_seu <- TRUE

# ---- Validate input ----
  
  if (!file.exists(seu_file)) {
    stop("Error: Seurat file does not exist.")
  }

# ---- Detect QC state if not provided ----
  
  if (is.null(qc_state_arg)) {
    
    base_path <- file.path(
      "03_result",
      dataset_name,
      "QC",
      "Cell_Counts"
    )
    
    before_path <- file.path(base_path, "Before_QC")
    
    if (dir.exists(before_path)) {
      qc_state <- "After_QC"
    } else {
      qc_state <- "Before_QC"
    }
    
  } else {
    
    if (!qc_state_arg %in% c("Before_QC", "After_QC")) {
      stop("Error: qc_state must be 'Before_QC' or 'After_QC'.")
    }
    
    qc_state <- qc_state_arg
  }

# ---- Define output directory ----
  
  if (is.null(outdir_arg)) {
    
    outdir <- file.path(
      "03_result",
      dataset_name,
      "QC",
      "Cell_Counts",
      qc_state
    )
    
  } else {
    
    outdir <- file.path(
      outdir_arg,
      dataset_name,
      "QC",
      "Cell_Counts",
      qc_state
    )
  }

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

cat("Seurat file path :", seu_file, "\n")
cat("Dataset name :", dataset_name, "\n")
cat("Output dir :", outdir, "\n")
cat("QC state :", qc_state, "\n")



#############
#--- LIB ---#
#############

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
})



# ------------------------------------------------------------
# 1. Table of cell counts
# ------------------------------------------------------------
qc_cell_count_table <- function(
    seu,
    group_col
) {
  stopifnot(group_col %in% colnames(seu@meta.data))
  seu@meta.data %>%
    dplyr::group_by(.data[[group_col]]) %>%
    dplyr::summarise(cell_count = dplyr::n(), .groups = "drop")
}

# ------------------------------------------------------------
# 2. Barplots of cell counts
# ------------------------------------------------------------
qc_cell_count_barplot <- function(
    seu,
    group_cols,
    qc_state,
    outdir = NULL,
    dataset_name = "dataset",
    palette = NULL
) {
  plots <- list()
  i <- 1
  
  if (!is.null(outdir)) {
    dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  }
  
  for (group_col in group_cols) {
    stopifnot(group_col %in% colnames(seu@meta.data))
    
    p <- ggplot(
      seu@meta.data,
      aes(x = .data[[group_col]], fill = .data[[group_col]])
    ) +
      geom_bar() +
      theme_classic() +
      theme(
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        plot.title = element_text(hjust = 0.5, face = "bold")
      ) +
      ggtitle(paste("NCells -", group_col)) +
      NoLegend()
    
    if (!is.null(palette)) {
      p <- p + scale_fill_manual(values = palette, drop = FALSE)
    }
    
    plots[[i]] <- p
    i <- i + 1
    
    if (!is.null(outdir)) {
      ggsave(
        plot = p,
        filename = file.path(
          outdir,
          paste0("Cell_Counts_By_", group_col, "_", dataset_name, "_", qc_state, "_bis.png")
        ),
        width = 6,
        height = 4
      )
    }
  }
  
  return(plots)
}

# ------------------------------------------------------------
# 3. UMAP facetted by group
# ------------------------------------------------------------
qc_cell_count_umap <- function(
    seu,
    group_col,
    outdir = NULL,
    qc_state,
    dataset_name = "dataset",
    reduction = "umap",
    ncol = 6
) {
  stopifnot(group_col %in% colnames(seu@meta.data))
  
  p <- DimPlot(seu, group.by = group_col, reduction = reduction) +
    facet_wrap(stats::as.formula(paste("~", group_col)), ncol = ncol) +
    ggtitle(paste("UMAP by", group_col, qc_state)) +
    NoLegend()
  
  if (!is.null(outdir)) {
    dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
    ggsave(
      plot = p,
      filename = file.path(
        outdir,
        paste0("UMAP_Cell_Counts_by_", group_col, "_", dataset_name, "_", qc_state, "_bis.png")
      ),
      width = 23,
      height = 12
    )
  }
  
  return(p)
}

# ------------------------------------------------------------
# 4️⃣ Fonction principale généralisée pour n'importe quel subset
# ------------------------------------------------------------
qc_cell_counts <- function(
    seu_subset,
    dataset_name = "dataset",
    qc_state,
    features_barplot = c("condition", "id", "emulsion"),
    umap_group_col = "emulsion",
    outdir = NULL,
    palette = NULL,
    reduction = "umap",
    umap_ncol = 6
) {
  all_plots <- list()
  i <- 1
  
  # --- Table résumé par emulsion ---
  print(qc_cell_count_table(seu_subset, umap_group_col))
  
  # --- Barplots ---
  barplots <- qc_cell_count_barplot(
    seu = seu_subset,
    group_cols = features_barplot,
    outdir = if (!is.null(outdir)) file.path(outdir, dataset_name, "QC", "Cell_Counts", qc_state) else NULL,
    dataset_name = dataset_name,
    palette = palette,
    qc_state = qc_state
  )
  
  for (p in barplots) {
    all_plots[[i]] <- p
    i <- i + 1
  }
  
  # --- UMAP facetted ---
  p_umap <- qc_cell_count_umap(
    seu = seu_subset,
    group_col = umap_group_col,
    outdir = if (!is.null(outdir)) file.path(outdir, dataset_name, "QC", "Cell_Counts", qc_state) else NULL,
    dataset_name = dataset_name,
    reduction = reduction,
    ncol = umap_ncol,
    qc_state = qc_state
  )
  
  all_plots[[i]] <- p_umap
  
  return(all_plots)
}

# Si le script est appelé directement, on exécute la fonction principale
if (!interactive()) {
  seu <- qs::qread(seu_file)
  qc_cell_counts(
    seu_subset = seu,
    dataset_name = dataset_name,
    outdir = outdir,
    qc_state = qc_state_arg,
    palette = NULL
  )
} else {
  cat("Script sourced. Define your Seurat object and call qc_cell_counts() with appropriate parameters.\n")
}
