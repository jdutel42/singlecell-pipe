# ============================================================
# 03_metrics_bio.R
# Métriques biologiques : qualité des markers, module scores, cohérence ADT/RNA
# ============================================================

# ── 1. Qualité des markers différentiels ─────────────────────────────────────
#
# Score biologique basé sur FindAllMarkers :
#   - Fraction de clusters avec ≥ N markers significatifs (couverture)
#   - log2FC médian des top markers
#   - Fraction de clusters sans markers (= clusters biologiquement vides)
#
# IMPORTANT : ne PAS appeler cette fonction dans le grid search complet —
# elle est coûteuse (Wilcoxon sur tout le dataset × n_clusters).
# L'utiliser sur les top-K candidats sélectionnés par les métriques mathématiques.
#
# @param seu          Objet Seurat avec clustering (Idents définis)
# @param assay        Assay à utiliser ("SCT" recommandé, ou "RNA")
# @param min_pct      Seuil min.pct pour FindAllMarkers
# @param lfc_thresh   Seuil logfc.threshold
# @param n_top        Nombre de top markers à considérer par cluster
# @return             Liste : $markers (data.frame complet),
#                             $score (score numérique ∈ [0,1]),
#                             $summary (tibble par cluster)
# ─────────────────────────────────────────────────────────────────────────────
compute_marker_quality <- function(seu, assay = "SCT",
                                   min_pct = 0.25, lfc_thresh = 0.5,
                                   n_top = 10) {
  # Exclure gènes mitochondriaux, ribosomaux (bruit technique)
  markers <- tryCatch(
    FindAllMarkers(
      seu,
      assay           = assay,
      slot            = "data",
      test.use        = "wilcox",
      only.pos        = TRUE,
      min.pct         = min_pct,
      logfc.threshold = lfc_thresh,
      verbose         = FALSE
    ),
    error = function(e) {
      warning("FindAllMarkers échoué : ", conditionMessage(e))
      return(NULL)
    }
  )

  if (is.null(markers) || nrow(markers) == 0) {
    return(list(markers = NULL, score = 0, summary = NULL))
  }

  markers_filt <- markers |>
    dplyr::filter(
      p_val_adj < 0.05,
      avg_log2FC > lfc_thresh,
      !grepl("^MT-|^RPL|^RPS|^MALAT1", gene)
    )

  n_clusters_total <- length(unique(Idents(seu)))

  summary_df <- markers_filt |>
    dplyr::group_by(cluster) |>
    dplyr::summarise(
      n_markers      = dplyr::n(),
      median_lfc     = median(avg_log2FC),
      top_genes      = paste(head(gene[order(-avg_log2FC)], n_top), collapse = ", "),
      .groups        = "drop"
    )

  # Fraction de clusters avec au moins 5 markers significatifs
  n_covered   <- sum(summary_df$n_markers >= 5)
  coverage    <- n_covered / n_clusters_total

  # Score : coverage × log2FC médian global normalisé (plafonné à 3)
  median_lfc_global <- median(markers_filt$avg_log2FC, na.rm = TRUE)
  lfc_score         <- min(median_lfc_global / 3, 1)

  score <- 0.6 * coverage + 0.4 * lfc_score

  list(
    markers = markers_filt,
    score   = round(score, 4),
    summary = summary_df
  )
}


# ── 2. Module scores (signatures biologiques) ────────────────────────────────
#
# AddModuleScore pour chaque signature définie dans BIO_SIGNATURES.
# Retourne un score de cohérence : dans quelle mesure les signatures
# sont-elles enrichies dans des clusters distincts (pas diffus) ?
#
# Logique du score de cohérence :
#   Pour chaque signature, on calcule le max du score moyen par cluster.
#   Si une signature est clairement enrichie dans UN cluster (max élevé,
#   variance inter-cluster élevée), c'est biologiquement cohérent.
#
# @param seu        Objet Seurat avec UMAP et clustering
# @param signatures Liste nommée de vecteurs de gènes (BIO_SIGNATURES)
# @return           Liste : $seu (avec colonnes ModuleScore_X),
#                           $score (cohérence globale ∈ [0,1]),
#                           $per_sig (data.frame : max_mean et spread par signature)
# ─────────────────────────────────────────────────────────────────────────────
compute_module_scores <- function(seu, signatures = BIO_SIGNATURES) {
  score_names <- character(length(signatures))

  for (i in seq_along(signatures)) {
    sig_name   <- names(signatures)[i]
    genes_sig  <- intersect(signatures[[i]], rownames(seu))

    if (length(genes_sig) < 2) {
      message(sprintf("  Signature '%s' : < 2 gènes trouvés, ignorée.", sig_name))
      next
    }

    col_name      <- paste0("ModuleScore_", sig_name)
    score_names[i] <- col_name

    seu <- AddModuleScore(
      seu,
      features = list(genes_sig),
      name     = col_name,
      seed     = 42
    )
    # AddModuleScore ajoute un suffixe "1" — on renomme proprement
    colnames(seu@meta.data)[colnames(seu@meta.data) == paste0(col_name, "1")] <- col_name
  }

  score_names <- score_names[score_names != ""]

  # Score de cohérence : spread inter-cluster normalisé
  meta <- seu@meta.data
  cluster_col <- "seurat_clusters"

  per_sig <- lapply(score_names, function(sc) {
    if (!sc %in% colnames(meta)) return(NULL)
    cluster_means <- meta |>
      dplyr::group_by(.data[[cluster_col]]) |>
      dplyr::summarise(mean_score = mean(.data[[sc]], na.rm = TRUE), .groups = "drop")

    spread   <- max(cluster_means$mean_score) - min(cluster_means$mean_score)
    max_mean <- max(cluster_means$mean_score)

    data.frame(
      signature = sc,
      spread    = round(spread, 4),
      max_mean  = round(max_mean, 4)
    )
  }) |> dplyr::bind_rows()

  # Score global : moyenne des spreads normalisés (plafonné à 1)
  coherence_score <- if (nrow(per_sig) > 0) {
    mean(pmin(per_sig$spread / 0.5, 1))  # 0.5 = spread "bon"
  } else 0

  list(
    seu             = seu,
    score           = round(coherence_score, 4),
    per_sig         = per_sig,
    score_col_names = score_names
  )
}


# ── 3. Cohérence RNA / ADT ───────────────────────────────────────────────────
#
# Pour chaque paire (ADT, RNA) du dictionnaire, on calcule la corrélation
# de Spearman entre l'expression ADT et l'expression RNA au niveau CLUSTER
# (moyennes par cluster). Une forte corrélation indique que les clusters
# sont cohérents entre les deux modalités.
#
# @param seu      Objet Seurat avec assays "ADT" et "RNA"/"SCT"
# @param adt_dict Dictionnaire nommé : nom_ADT → nom_RNA (ADT_RNA_DICT)
# @param rna_assay  Assay RNA à utiliser
# @return         Liste : $cor_df (corrélations par paire),
#                         $score (corrélation médiane globale ∈ [0,1])
# ─────────────────────────────────────────────────────────────────────────────
compute_adt_rna_coherence <- function(seu, adt_dict = ADT_RNA_DICT,
                                      rna_assay = "SCT") {
  cluster_col <- "seurat_clusters"
  clusters    <- seu@meta.data[[cluster_col]]

  # Moyennes par cluster pour ADT et RNA
  adt_means <- tryCatch({
    DefaultAssay(seu) <- "ADT"
    adt_features <- intersect(names(adt_dict), rownames(seu[["ADT"]]))
    AverageExpression(seu, assays = "ADT",
                      features = adt_features,
                      group.by = cluster_col,
                      verbose = FALSE)$ADT
  }, error = function(e) { warning("ADT AverageExpression échoué"); NULL })

  rna_means <- tryCatch({
    DefaultAssay(seu) <- rna_assay
    rna_features <- unique(unlist(adt_dict))
    rna_features <- intersect(rna_features, rownames(seu[[rna_assay]]))
    AverageExpression(seu, assays = rna_assay,
                      features = rna_features,
                      group.by = cluster_col,
                      verbose = FALSE)[[rna_assay]]
  }, error = function(e) { warning("RNA AverageExpression échoué"); NULL })

  if (is.null(adt_means) || is.null(rna_means)) {
    return(list(cor_df = NULL, score = NA_real_))
  }

  cor_df <- lapply(names(adt_dict), function(adt_name) {
    rna_name <- adt_dict[[adt_name]]
    if (!adt_name %in% rownames(adt_means)) return(NULL)
    if (!rna_name %in% rownames(rna_means)) return(NULL)

    adt_vec <- as.numeric(adt_means[adt_name, ])
    rna_vec <- as.numeric(rna_means[rna_name, ])

    cor_val <- tryCatch(
      cor(adt_vec, rna_vec, method = "spearman"),
      error = function(e) NA_real_
    )

    data.frame(
      adt_feature = adt_name,
      rna_feature = rna_name,
      spearman_r  = round(cor_val, 4)
    )
  }) |> dplyr::bind_rows()

  median_cor <- median(cor_df$spearman_r, na.rm = TRUE)

  # Normalise [-1,1] → [0,1]
  score <- (median_cor + 1) / 2

  list(
    cor_df = cor_df,
    score  = round(score, 4)
  )
}


# ── 4. Score biologique composite ────────────────────────────────────────────
#
# Combine les trois métriques biologiques en un score unique.
# Appelé uniquement sur les candidats présélectionnés.
#
# @param marker_score     Score qualité markers ∈ [0,1]
# @param module_score     Score cohérence module scores ∈ [0,1]
# @param adt_rna_score    Score cohérence ADT/RNA ∈ [0,1]
# @param w_markers        Poids markers (défaut 0.5)
# @param w_modules        Poids module scores (défaut 0.3)
# @param w_adt            Poids cohérence ADT/RNA (défaut 0.2)
# @return                 Score biologique ∈ [0,1]
# ─────────────────────────────────────────────────────────────────────────────
compute_bio_score <- function(marker_score, module_score, adt_rna_score,
                              w_markers = 0.5, w_modules = 0.3, w_adt = 0.2) {
  # Traiter les NA comme 0 (métrique non disponible)
  scores <- c(
    w_markers * ifelse(is.na(marker_score), 0, marker_score),
    w_modules * ifelse(is.na(module_score), 0, module_score),
    w_adt     * ifelse(is.na(adt_rna_score), 0, adt_rna_score)
  )
  round(sum(scores), 4)
}
