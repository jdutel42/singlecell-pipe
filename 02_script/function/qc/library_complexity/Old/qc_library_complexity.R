# ==============================================================
# QC – Library Complexity (CD4 T cells)
# Encapsulated functions (same pattern as other QC blocks)
# ==============================================================\n
library(dplyr)
library(ggplot2)
library(scater)

# --------------------------------------------------------------
# 1. Compute library complexity + outlier status
# --------------------------------------------------------------
compute_libcomplex_qc <- function(seu,
                                  group_var = "emulsion",
                                  nmads = 3) {
  
  meta <- seu@meta.data
  
  meta$log10_Complexity <- log10(meta$nFeature_RNA) / log10(meta$nCount_RNA)
  
  meta <- meta %>%
    dplyr::group_by(.data[[group_var]]) %>%
    dplyr::mutate(
      libcomplex_status = dplyr::case_when(
        isOutlier(log10_Complexity, nmads = nmads, type = "lower") ~ "low",
        isOutlier(log10_Complexity, nmads = nmads, type = "higher") ~ "high",
        TRUE ~ "normal"
      )
    ) %>%
    dplyr::ungroup()
  
  seu <- AddMetaData(seu, meta$log10_Complexity, col.name = "log10_Complexity")
  seu <- AddMetaData(seu, meta$libcomplex_status, col.name = "libcomplex_status")
  
  return(seu)
}

# --------------------------------------------------------------
# 2. Density plot (global)
# --------------------------------------------------------------
plot_libcomplex_density <- function(seu, outdir, prefix) {
  
  p <- ggplot(seu@meta.data, aes(x = log10_Complexity)) +
    geom_density(fill = "grey80", alpha = 0.6) +
    theme_classic() +
    labs(
      title = "Library complexity",
      x = "log10(Genes / UMI)",
      y = "Density"
    )
  
  ggsave(
    filename = file.path(outdir, paste0(prefix, "_Density_Library_Complexity.png")),
    plot = p,
    width = 10,
    height = 5
  )
  
  return(p)
}

# --------------------------------------------------------------
# 3. Scatter plot nCount vs nFeature (colored by complexity)
# --------------------------------------------------------------
plot_libcomplex_scatter <- function(seu, outdir, prefix) {
  
  p <- ggplot(seu@meta.data,
              aes(x = nCount_RNA,
                  y = nFeature_RNA,
                  color = libcomplex_status)) +
    geom_point(alpha = 0.6, size = 0.5) +
    geom_smooth(method = "lm", color = "black", linewidth = 0.6) +
    scale_x_log10() +
    scale_y_log10() +
    scale_color_manual(values = c(
      low = "blue3",
      normal = "grey70",
      high = "red3"
    )) +
    theme_classic() +
    labs(
      title = "Library complexity",
      x = "UMIs per cell",
      y = "Genes per cell",
      color = "Complexity"
    )
  
  ggsave(
    filename = file.path(outdir, paste0(prefix, "_Scatter_Library_Complexity.png")),
    plot = p,
    width = 10,
    height = 5
  )
  
  return(p)
}

# --------------------------------------------------------------
# 4. Density by group (emulsion)
# --------------------------------------------------------------
plot_libcomplex_density_by_group <- function(seu,
                                             group_var = "emulsion",
                                             outdir,
                                             prefix) {
  
  p <- ggplot(seu@meta.data,
              aes(x = log10_Complexity, fill = .data[[group_var]])) +
    geom_density(alpha = 0.4) +
    facet_wrap(as.formula(paste("~", group_var))) +
    theme_classic() +
    theme(legend.position = "none") +
    labs(
      title = "Library complexity by group",
      x = "log10(Genes / UMI)"
    )
  
  ggsave(
    filename = file.path(outdir, paste0(prefix, "_Density_Library_Complexity_By_", group_var, ".png")),
    plot = p,
    width = 12,
    height = 6
  )
  
  return(p)
}

# --------------------------------------------------------------
# 5. Violin plot by group
# --------------------------------------------------------------
plot_libcomplex_violin <- function(seu,
                                   group_var = "emulsion",
                                   outdir,
                                   prefix) {
  
  p <- ggplot(seu@meta.data,
              aes(x = .data[[group_var]],
                  y = log10_Complexity,
                  color = libcomplex_status)) +
    geom_violin(fill = "grey90", color = NA) +
    geom_quasirandom(alpha = 0.4, size = 0.6) +
    scale_color_manual(values = c(
      low = "blue",
      normal = "green",
      high = "red"
    )) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(
      title = "Library complexity by group",
      y = "log10(Genes / UMI)",
      x = NULL
    )
  
  ggsave(
    filename = file.path(outdir, paste0(prefix, "_Violin_Library_Complexity_By_", group_var, ".png")),
    plot = p,
    width = 12,
    height = 6
  )
  
  return(p)
}

# --------------------------------------------------------------
# Example run
# --------------------------------------------------------------
# seu_CD4 <- compute_libcomplex_qc(seu_CD4)
# plot_libcomplex_density(seu_CD4, OUTDIR, "CD4")
# plot_libcomplex_scatter(seu_CD4, OUTDIR, "CD4")
# plot_libcomplex_density_by_group(seu_CD4, "emulsion", OUTDIR, "CD4")
# plot_libcomplex_violin(seu_CD4, "emulsion", OUTDIR, "CD4")
