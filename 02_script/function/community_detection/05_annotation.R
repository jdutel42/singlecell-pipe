# ============================================================
# 05_annotation.R
# Visualisations d'annotation : FeaturePlot ADT/RNA, heatmap markers,
# volcano plots — appelées sur le clustering FINAL choisi
# ============================================================


# ── 1. FeaturePlot RNA + ADT par paire ───────────────────────────────────────
#
# CORRECTION vs code original :
#   - opérait sur `seu_CD4_QC_Moderate` (hardcodé) au lieu de l'objet courant
#   - on passe maintenant l'objet explicitement
#
# @param seu       Objet Seurat final (avec UMAP et clustering)
# @param adt_dict  ADT_RNA_DICT
# @param out_dir   Dossier de sortie
# ─────────────────────────────────────────────────────────────────────────────
plot_adt_rna_feature <- function(seu, adt_dict = ADT_RNA_DICT, out_dir) {
  make_dir(out_dir)
  n_clust    <- length(unique(Idents(seu)))
  colors_use <- pick_colors(n_clust)

  for (prot in names(adt_dict)) {
    gene <- adt_dict[[prot]]

    DefaultAssay(seu) <- "ADT"
    p_adt <- FeaturePlot(seu, features = prot, reduction = "umap") +
      ggtitle(paste0(prot, " (ADT)"))

    DefaultAssay(seu) <- "RNA"
    if (gene %in% rownames(seu[["RNA"]])) {
      p_rna <- FeaturePlot(seu, features = gene, reduction = "umap") +
        ggtitle(paste0(gene, " (RNA)"))
    } else {
      p_rna <- ggplot() + ggtitle(paste0(gene, " (RNA) – non trouvé")) + theme_void()
    }

    p_clust <- DimPlot(seu, reduction = "umap", group.by = "seurat_clusters") +
      ggtitle("Clusters") +
      scale_color_manual(values = colors_use)

    p <- (p_adt | p_rna) / p_clust

    save_plot(p,
              file.path(out_dir, paste0("ADT_RNA_", prot, "_", gene)),
              width = 12, height = 14)
  }
  invisible(NULL)
}


# ── 2. Heatmap top markers ────────────────────────────────────────────────────
#
# @param seu       Objet Seurat final
# @param markers   data.frame issu de FindAllMarkers (filtré)
# @param assay     Assay pour scale.data
# @param n_top     N top markers par cluster
# @param out_dir   Dossier de sortie
# ─────────────────────────────────────────────────────────────────────────────
plot_marker_heatmap <- function(seu, markers, assay = "SCT", n_top = 10, out_dir) {
  make_dir(out_dir)

  genes_scaled <- rownames(seu[[assay]]@scale.data)
  top_n <- markers |>
    dplyr::filter(gene %in% genes_scaled) |>
    dplyr::group_by(cluster) |>
    dplyr::slice_max(avg_log2FC, n = n_top) |>
    dplyr::arrange(cluster, dplyr::desc(avg_log2FC))

  genes    <- unique(top_n$gene)
  clusters <- seu$seurat_clusters
  ord      <- order(clusters)

  mat <- GetAssayData(seu, assay = assay, layer = "scale.data")[genes, ord]
  mat <- pmin(pmax(mat, -2.5), 2.5)   # cap pour la lisibilité

  cluster_ha <- ComplexHeatmap::HeatmapAnnotation(
    Cluster = clusters[ord],
    col = list(
      Cluster = structure(
        pick_colors(length(unique(clusters))),
        names = levels(factor(clusters))
      )
    )
  )

  ht <- ComplexHeatmap::Heatmap(
    mat,
    name            = "Scaled expr.",
    top_annotation  = cluster_ha,
    show_column_names = FALSE,
    show_row_names    = TRUE,
    cluster_columns   = FALSE,
    cluster_rows      = TRUE,
    row_names_gp      = grid::gpar(fontsize = 8),
    column_title      = paste0("Top ", n_top, " markers per cluster")
  )

  for (ext in c("pdf", "png")) {
    if (ext == "pdf") {
      pdf(file.path(out_dir, "top_markers_heatmap.pdf"), width = 12, height = 10)
    } else {
      png(file.path(out_dir, "top_markers_heatmap.png"),
          width = 2000, height = 1600, res = 300)
    }
    ComplexHeatmap::draw(ht)
    dev.off()
  }

  invisible(ht)
}


# ── 3. Volcano plots par cluster ─────────────────────────────────────────────
#
# @param markers   data.frame complet de FindAllMarkers
# @param out_dir   Dossier de sortie
# @param n_label   N gènes à annoter (up + down)
# ─────────────────────────────────────────────────────────────────────────────
plot_volcano_per_cluster <- function(markers, out_dir, n_label = 10) {
  make_dir(out_dir)

  for (i in unique(markers$cluster)) {
    df <- markers |>
      dplyr::filter(cluster == i) |>
      dplyr::mutate(
        Status = dplyr::case_when(
          p_val_adj < 0.05 & avg_log2FC >  0.5 ~ "Up",
          p_val_adj < 0.05 & avg_log2FC < -0.5 ~ "Down",
          TRUE                                  ~ "NS"
        )
      )

    top_genes <- dplyr::bind_rows(
      df |> dplyr::filter(Status == "Up")   |> dplyr::arrange(dplyr::desc(avg_log2FC)) |> dplyr::slice_head(n = n_label),
      df |> dplyr::filter(Status == "Down") |> dplyr::arrange(avg_log2FC)              |> dplyr::slice_head(n = n_label)
    )

    p <- ggplot(df, aes(avg_log2FC, -log10(p_val_adj + 1e-300))) +
      geom_point(aes(color = Status), size = 1, alpha = 0.7) +
      geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed") +
      geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
      ggrepel::geom_text_repel(data = top_genes, aes(label = gene),
                               size = 3, max.overlaps = Inf) +
      scale_color_manual(values = c(Up = "#D62728", Down = "#1F77B4", NS = "grey70")) +
      labs(title  = paste("Cluster", i, "— DE genes"),
           x      = "avg log2FC",
           y      = "-log10(adj p-value)",
           color  = "Regulation") +
      theme_classic(base_size = 12)

    save_plot(p, file.path(out_dir, paste0("volcano_cluster_", i)), 8, 6)
  }
  invisible(NULL)
}


# ── 4. Module score UMAP overlay ─────────────────────────────────────────────
#
# @param seu         Objet Seurat avec colonnes ModuleScore_X dans @meta.data
# @param score_names Noms des colonnes de scores (retourné par compute_module_scores)
# @param out_dir     Dossier de sortie
# ─────────────────────────────────────────────────────────────────────────────
plot_module_score_umaps <- function(seu, score_names, out_dir) {
  make_dir(out_dir)

  for (sc in score_names) {
    if (!sc %in% colnames(seu@meta.data)) next

    p <- FeaturePlot(seu, features = sc, reduction = "umap",
                     cols = c("lightgrey", "#D62728"), order = TRUE) +
      ggtitle(gsub("ModuleScore_", "", sc))

    save_plot(p, file.path(out_dir, paste0("ModuleScore_UMAP_", sc)), 7, 6)
  }
  invisible(NULL)
}
