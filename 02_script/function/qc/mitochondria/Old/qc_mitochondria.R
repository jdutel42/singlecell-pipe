# ============================================================
# QC Mitochondrial content - scRNA-seq
# Seurat - modulaire
# ============================================================


##############
#--- ARGS ---#
##############

args <- commandArgs(trailingOnly = TRUE)

seu_file <- args[1]         # ex: "path/to/seu_CD4.rds"
dataset_name <- args[2]     # ex: "CD4_Tcells"
outdir <- args[3]           # ex: "path/to/results"
datadir <- args[4]          # ex: "path/to/data"
save_seu <- TRUE            # par défaut TRUE

cat("Seurat file path   :", seu_file, "\n")
cat("Dataset name  :", dataset_name, "\n")
cat("Output dir    :", outdir, "\n")
cat("Data dir      :", datadir, "\n")
cat("Save Seurat   :", save_seu, "\n")

##############
#--- HELP ---#
##############

# Fonction pour afficher l'aide
print_help <- function() {
  cat("
Usage: Rscript qc_mitochondria.R <seu_file> <dataset_name> <outdir> <save_seu>

Arguments:
  seu_file      : chemin vers le fichier RDS de ton Seurat object
  dataset_name  : nom du dataset (ex: CD4_Tcells, CD8_Tcells)
  outdir        : dossier où sauvegarder les plots et l'objet annoté
  datadir       : chemin vers le dossier de données (ex: path/to/data)
  save_seu      : TRUE ou FALSE, pour sauvegarder l'objet Seurat annoté

Exemple:
  Rscript qc_mitochondria.R ./seu_CD4.rds CD4_Tcells ./Results ./Data TRUE
  ")
}

# Si aucun argument ou --help, afficher l'aide
if (length(args) == 0 || length(args) < 4 || "--help" %in% args) {
  print_help()
  quit(status = 0)
}

#############
#--- LIB ---#
#############

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(ggridges)
  library(scater)
  library(ggbeeswarm)  # pour geom_quasirandom
  library(viridis)
})

run_qc_mito <- function(seu_file, dataset_name, outdir, datadir, save_seu = TRUE) {

  # ---- Charger le subset Seurat ----
  seu <- qs::qread(seu_file)
  

  # ---- Calcul mitochondrie ----
  counts <- seu@assays$RNA@counts
  qc_mito <- scater::perCellQCMetrics(
    counts,
    subsets = list(mito = grep("^MT-", rownames(counts), value = TRUE))
  )
  seu <- AddMetaData(seu, as.data.frame(qc_mito))

  qc_mito2 <- data.frame(
    emulsion = seu@meta.data$emulsion,
    subsets_mito_percent = seu@meta.data$subsets_mito_percent
  )

  qc_mito2$high_mito_outliers <- scater::isOutlier(
    qc_mito2$subsets_mito_percent,
    nmads = 3,
    type = "higher"
  )

  qc_mito2 <- qc_mito2 %>%
    dplyr::group_by(emulsion) %>%
    dplyr::mutate(
      mito_status = dplyr::case_when(
        as.logical(high_mito_outliers) ~ "high",
        TRUE ~ "normal"
      )
    ) %>%
    dplyr::ungroup() %>%
    filter(!is.na(mito_status))

  seu <- AddMetaData(seu, qc_mito2$mito_status, col.name = "mito_status")

  # ---- Créer dossier de sauvegarde ----
  save_dir <- file.path(outdir, dataset_name, "QC", "Mitochondria", "Before_QC")
  dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)

  # ---- Ridgeplot ----
  p1 <- ggplot(seu@meta.data, aes(x = subsets_mito_percent, y = emulsion, fill = after_stat(x))) +
    ggridges::geom_density_ridges_gradient(scale = 3, rel_min_height = 0.01) +
    viridis::scale_fill_viridis(alpha = 0.9, option = "C", guide = "none") +
    labs(title = "Mitochondrial content by emulsion", x = "Percent mitochondrial reads", y = NULL) +
    theme_minimal(base_size = 12) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.spacing = unit(0.1, "lines"),
          strip.text.x = element_text(size = 8),
          plot.title = element_text(face = "bold"))

  ggsave(file.path(save_dir, paste0("Ridgeplot_mitochondrial_content_by_emulsion_before_QC_", dataset_name, "_bis.png")),
         plot = p1, width = 8, height = 5)

  # ---- FeaturePlot UMAP ----
  p2 <- FeaturePlot(seu, features = "subsets_mito_percent", reduction = "umap", order = TRUE, cols = c("lightgrey", "red"))
  ggsave(file.path(save_dir, paste0("UMAP_mitochondrial_content_before_QC_", dataset_name, "_bis.png")),
         plot = p2, width = 6, height = 5)

  # ---- Violin + jitter ----
  p3 <- ggplot(seu@meta.data, aes(x = emulsion, y = subsets_mito_percent, color = mito_status)) +
    geom_violin(fill = "grey90", color = NA) +
    ggbeeswarm::geom_quasirandom(alpha = 0.5, size = 0.8) +
    scale_y_log10() +
    scale_color_manual(values = c("normal" = "green3", "high" = "red3")) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

  ggsave(file.path(save_dir, paste0("ViolinPlot_mito_percent_by_emulsion_", dataset_name, "_bis.png")),
         plot = p3, width = 12, height = 6)

  # ---- Sauvegarder Seurat annoté ----
  if (save_seu) {
    save_file <- file.path(datadir, "Processed_data", paste0("seu_", dataset_name,"_QC_mitochondria.qs"))
    qs::qsave(seu, file = save_file)
    # saveRDS(seu, file = file.path(save_dir, paste0("seu_", dataset_name,"_QC_mitochondria.rds")))
    cat("Seurat annoté sauvegardé ici :", save_file, "\n")
  }

  return(TRUE)
}

run_qc_mito(seu_file, dataset_name, outdir, datadir, save_seu)

# ---- Exécuter si script lancé via callr ----
if (interactive()) {
  run_qc_mito(seu_file, dataset_name, outdir, datadir, save_seu)
  cat("Script exécuté en mode interactif. Vérifie les arguments et la fonction run_qc_mito() pour tester manuellement.\n")
}