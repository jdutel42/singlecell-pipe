#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(qs)
  library(Seurat)
  library(lisi)
  library(dplyr)
  library(ggplot2)
  library(stringr)
})

##################
#--- FUNCTION ---#
##################

run_lisi <- function(input_dir, dataset_name, batch_col, outdir, condition_col) {
  
  ########################
  # OUTPUT DIR
  ########################
  
  save_dir <- file.path(
    outdir,
    dataset_name,
    "Integration",
    "LISI"
  )
  
  dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)
  
  # ============================================================
  # FIND FILES
  # ============================================================
  
  files <- list.files(
    input_dir,
    pattern="\\.qs$",
    full.names=TRUE
  )
  
  results <- list()
  
  # ============================================================
  # LOOP FILES
  # ============================================================
  
  for(f in files){
    
    message("Processing: ", f)
    
    seu <- qs::qread(f)
    
    emb <- Embeddings(seu, "harmony")
    
    meta <- seu@meta.data
    
    df <- data.frame(
      batch = meta[[batch_col]],
      condition = meta[[condition_col]]
    )
    
    lisi_scores <- compute_lisi(
      emb,
      df,
      c("batch","condition")
    )
    
    ilisi <- median(lisi_scores$batch)
    clisi <- median(lisi_scores$condition)
    score <- ilisi / clisi
    
    method <- basename(f) |> str_remove(".qs")
    
    results[[method]] <- data.frame(
      method = method,
      iLISI = ilisi,
      cLISI = clisi,
      Score = score
    )
    
  }
  
  results_df <- bind_rows(results)
  
  # ============================================================
  # SAVE SCORES
  # ============================================================
  
  write.csv(
    results_df,
    file.path(save_dir,"lisi_scores.csv"),
    row.names=FALSE
  )
  
  # ============================================================
  # PLOT
  # ============================================================
  
  p <- ggplot(results_df,
              aes(x=iLISI, y=cLISI, label=method)) +
    geom_point(size=4) +
    geom_text(vjust=-1) +
    theme_classic() +
    labs(
      title="Batch integration benchmark",
      x="Batch mixing (iLISI)",
      y="Condition conservation (cLISI)"
    )
  
  ggsave(
    file.path(save_dir,"lisi_comparison.png"),
    p,
    width=7,
    height=5,
    dpi=300
  )
  
  ggsave(
    file.path(save_dir,"lisi_comparison.pdf"),
    p,
    width=7,
    height=5
  )
  
  message("LISI computation finished")
  
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
      c("-i","--input_dir"), 
      type="character",
      help="Directory containing integration results"),
    
    make_option(
      c("-d","--dataset"),
      type="character",
      help="Dataset name"),
    
    make_option(
      c("-o","--outdir"), 
      type="character",
      default="03_result/lisi_scores"),
    
    make_option(
      c("-b","--batch_col"), 
      type="character",
      default="emulsion"),
    
    make_option(
      c("-c","--condition_col"), 
      type="character",
      default="condition")
    
  )
  
  opt <- parse_args(OptionParser(option_list=option_list))
  
  # ---- Logging ----
  message("\n========== Pipeline Configuration ==========")
  message("Input dir : ", opt$input_dir)
  message("Dataset    : ", opt$dataset)
  message("Batch col  : ", opt$batch_col)
  message("Output dir : ", opt$outdir)
  message("Condition col : ", opt$condition_col)
  message("===========================================\n")
  
  # ---- Run QC Mitochondria ----
  run_lisi(opt$input_dir, opt$dataset, opt$batch_col, opt$outdir, opt$condition_col)
  
}

################
#--- SOURCE ---#
################

if (sys.nframe() == 0) {
  main()
}