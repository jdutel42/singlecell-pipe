qc_cell_counts <- function(seu,
                           group_vars,
                           palette,
                           outdir,
                           prefix) {
  
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  
  for (var in group_vars) {
    p1 <- seu@meta.data %>%
      ggplot(aes(x = .data[[var]], fill = .data[[var]])) +
      geom_bar() +
      scale_fill_manual(values = palette, drop = FALSE) +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      ggtitle(paste("NCells -", var))
    
    ggsave(
      file.path(outdir, paste0(prefix, "_CellCounts_by_", var, ".png")),
      p1,
      width = 6,
      height = 4
    )
  }
  
  p2 <- DimPlot(seu, group.by = param$batch_var) +
    facet_wrap(~param$batch_var, ncol = 4) +
    ggtitle(paste0("UMAP by ", param$batch_var, " before QC CD4 Tcells")) + 
    NoLegend()
  
  ggsave(
    file.path(outdir, paste0("UMAP by ", param$batch_var, " before QC CD4 Tcells.png"), 
                             p2,
                             width = 10, 
                             height = 8)
    )
}
