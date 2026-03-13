# ============================================================
# QC Library size (nCount_RNA) - scRNA-seq
# Seurat + scater
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(scater)
  library(dplyr)
  library(ggplot2)
  library(ggbeeswarm)
})

# ------------------------------------------------------------
# 1. Compute library size QC using isOutlier (low + high)
# ------------------------------------------------------------
qc_libsize_compute <- function(
    seu,
    batch_col,
    feature = "nCount_RNA",
    nmads = 3,
    meta_colname = "libsize_status"
) {
  
  stopifnot(batch_col %in% colnames(seu@meta.data))
  stopifnot(feature %in% colnames(seu@meta.data))
  
  meta <- seu@meta.data %>%
    dplyr::filter(.data[[feature]] > 0)
  
  libsize <- meta[[feature]]
  
  meta <- meta %>%
    dplyr::mutate(
      high_outlier = isOutlier(
        libsize,
        nmads = nmads,
        type = "higher"
      ),
      low_outlier = isOutlier(
        libsize,
        nmads = nmads,
        type = "lower"
      ),
      !!meta_colname := dplyr::case_when(
        high_outlier ~ "high",
        low_outlier  ~ "low",
        TRUE         ~ "normal"
      )
    )
  
  seu <- AddMetaData(
    seu,
    metadata = meta[[meta_colname]],
    col.name = meta_colname
  )
  
  return(seu)
}

# ------------------------------------------------------------
# 2. Visualization
# ------------------------------------------------------------
qc_libsize_plots <- function(
    seu,
    batch_col,
    outdir,
    dataset_name = "dataset",
    feature = "nCount_RNA",
    meta_status = "libsize_status"
) {
  
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  
  meta <- seu@meta.data
  
  # ---- Density plot
  p1 <- ggplot(meta, aes(x = .data[[feature]])) +
    geom_density(fill = "grey80", alpha = 0.6) +
    geom_rug(
      aes(color = .data[[meta_status]]),
      alpha = 0.6,
      sides = "b"
    ) +
    scale_x_log10() +
    scale_color_manual(
      values = c(
        low = "blue",
        normal = "green",
        high = "red"
      )
    ) +
    facet_wrap(
      stats::as.formula(paste("~", batch_col))
    ) +
    theme_classic() +
    theme(
      axis.text.x = element_text(
        angle = 45,
        vjust = 1,
        hjust = 1
      )
    ) +
    labs(
      color = "Library size QC",
      x = "nCount_RNA (log10)",
      y = "Density"
    )
  
  ggsave(
    plot = p1,
    filename = file.path(
      outdir,
      paste0(
        "DensityPlot_library_size_by_",
        batch_col,
        "_",
        dataset_name,
        ".png"
      )
    ),
    width = 12,
    height = 6
  )
  
  # ---- Violin plot
  p2 <- ggplot(
    meta,
    aes(
      x = .data[[batch_col]],
      y = .data[[feature]],
      color = .data[[meta_status]]
    )
  ) +
    geom_violin(fill = "grey90", color = NA) +
    geom_quasirandom(alpha = 0.5, size = 0.8) +
    scale_y_log10() +
    scale_color_manual(
      values = c(
        low = "blue",
        normal = "green",
        high = "red"
      )
    ) +
    theme_classic() +
    theme(
      axis.text.x = element_text(
        angle = 45,
        vjust = 1,
        hjust = 1
      )
    )
  
  ggsave(
    plot = p2,
    filename = file.path(
      outdir,
      paste0(
        "ViolinPlot_library_size_by_",
        batch_col,
        "_",
        dataset_name,
        ".png"
      )
    ),
    width = 12,
    height = 6
  )
  
  invisible(list(density = p1, violin = p2))
}
