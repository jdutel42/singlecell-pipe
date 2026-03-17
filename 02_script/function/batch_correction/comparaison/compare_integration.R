#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(readr)
  library(stringr)
  library(ggrastr)
  library(ggrepel)
})

##############
#--- FUNCTION ---#
##############

compare_metrics <- function(input_dirs, dataset_name, outdir, weights = c(iMetric=0.5, cMetric=0.5)) {
  

  
  # Fonction pour lire et nettoyer un CSV
  read_clean_csv <- function(file){
    df <- read_csv(file, col_names = TRUE, quote = "\"", show_col_types = FALSE)  # read_csv gère les guillemets
    names(df) <- str_replace_all(names(df), '"', '')      # supprime les " dans les noms
    df
  }
  
  # Dossiers contenant les CSV
  lisi_dir <- "LISI"
  purity_dir <- "KNN_purity"
  
  lisi_files <- file.path(input_dirs, lisi_dir, "lisi_scores.csv")
  purity_files <- file.path(input_dirs, purity_dir, "knn_purity_scores.csv")

  # Lire tous les CSV LISI
  lisi <- lapply(lisi_files, read_clean_csv) %>% bind_rows() %>%
    mutate(
      iLISI_scaled = (iLISI - min(iLISI)) / (max(iLISI) - min(iLISI)),
      cLISI_scaled = (cLISI - min(cLISI)) / (max(cLISI) - min(cLISI)),
      bacth_mixing_lisi = iLISI_scaled,
      condition_mixing_lisi = cLISI_scaled,
      Score_integration_LISI_scaled = bacth_mixing_lisi + (1 - condition_mixing_lisi)
    )
  
  # Lire tous les CSV Purity
  purity <- lapply(purity_files, read_clean_csv) %>% bind_rows() %>%
    mutate(
      iPurity_scaled = (iPurity - min(iPurity)) / (max(iPurity) - min(iPurity)),
      cPurity_scaled = (cPurity - min(cPurity)) / (max(cPurity) - min(cPurity)),
      batch_mixing_purity = 1 - iPurity_scaled,
      condition_mixing_purity = cPurity_scaled,
      Score_integration_Purity_scaled = batch_mixing_purity + (1 - condition_mixing_purity)
    )
  
  
  
  # Fusionner par "method"
  combined <- merge(lisi, purity, by="method") %>% 
    mutate(
      Batch_mixing_global = rowMeans(cbind(bacth_mixing_lisi, batch_mixing_purity)),
      Condition_mixing_global = rowMeans(cbind(condition_mixing_lisi, condition_mixing_purity)),
      Score_integration_global = (Score_integration_LISI_scaled + Score_integration_Purity_scaled) / 2
    )
  
  # Meilleur modèle
  best_model <- combined %>% filter(Score_integration_global == max(Score_integration_global))
  
  message("Best model: ", best_model)
  
  save_dir <- file.path(outdir, dataset_name, "Integration", "Comparaison")
  dir.create(file.path(save_dir), recursive = TRUE, showWarnings = FALSE)
  
  # --- Sauvegarde ---
  write.csv(combined, file.path(save_dir,
                                "combined_metrics.csv"), row.names=FALSE)
  
  # --- Plot ---
  p1 <- ggplot(combined, aes(x=Score_integration_LISI_scaled, y=Score_integration_Purity_scaled, label=method)) +
    geom_point_rast(aes(color = Score_integration_global), size = 4) +
    geom_text_repel(size = 2) +
    scale_color_viridis_c(option="plasma") +
    coord_cartesian(clip = "off") +
    theme_classic() +
    theme(
      plot.margin = margin(10, 30, 10, 10)
    ) +
    labs(
      x="Global integration score (LISI)",
      y="Global integration score (Purity)",
      color="Score_integration_global",
      title="Comparaison synthétique des méthodes"
    )
  
  ggsave(file.path(save_dir,
                   "metric_comparison.png"), p1, width=8, height=6, dpi=300)
  
  ggsave(file.path(save_dir,
                   "metric_comparison.pdf"), p1, width=8, height=6)
  
  p2 <- ggplot(combined, aes(x=Batch_mixing_global, y=Condition_mixing_global, label=method)) +
    geom_point_rast(aes(color = Score_integration_global), size = 4) +
    geom_text_repel(size = 2) +
    scale_color_viridis_c(option="plasma") +
    coord_cartesian(clip = "off") +
    theme_classic() +
    theme(
      plot.margin = margin(10, 30, 10, 10)
    ) +
    labs(
      x="Batch_mixing",
      y="Condition_mixing",
      color="Score_integration_global",
      title="Comparaison synthétique des méthodes"
    )
  
  ggsave(file.path(save_dir,
                   "metric_comparison_iMetric.png"), p2, width=8, height=6, dpi=300)
  
  ggsave(file.path(save_dir,
                   "metric_comparison_iMetric.pdf"), p2, width=8, height=6)
  
  p3 <- ggplot(combined, aes(x=reorder(method, Score_integration_global), y=Score_integration_global, fill=Score_integration_global)) +
    geom_col() +
    coord_flip() +
    scale_fill_viridis_c(option="plasma") +
    theme_minimal() +
    labs(x="Méthode", y="Score Global", title="Synthèse des méthodes")
  
  ggsave(file.path(save_dir,
                   "metric_comparison_bar.png"), p3, width=8, height=6, dpi=300)
  
  ggsave(file.path(save_dir,
                   "metric_comparison_bar.pdf"), p3, width=8, height=6)
  
  return(list(combined_table=combined, best_method=best_model))
}


##############
#--- MAIN ---#
##############

main <- function(){
  
  option_list <- list(
    make_option(c("-i","--input_dirs"), type="character", help="Comma-separated list of metric folders"),
    make_option(c("-d","--dataset_name"), type="character", help="Dataset name"),
    make_option(c("-o","--outdir"), type="character", default="03_result/metric_comparison")
  )
  
  opt <- parse_args(OptionParser(option_list=option_list))
  
  dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)
  
  input_dirs <- str_split(opt$input_dirs, ",")[[1]] |> str_trim()
  
  res <- compare_metrics(input_dirs, opt$dataset_name, opt$outdir)
  
  message("Metric comparison finished. Best method: ", res$best_method$method)
  
}

if (sys.nframe() == 0) main()