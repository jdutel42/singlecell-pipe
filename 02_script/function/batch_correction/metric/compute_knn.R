#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(qs)
  library(Seurat)
  library(FNN)      # pour get.knn
  library(dplyr)
  library(stringr)
})

##################
#--- FUNCTION ---#
##################

run_knn_purity <- function(input_dir, dataset_name, batch_col, condition_col, outdir, k=30){
  
  save_dir <- file.path(outdir, dataset_name, "Integration", "KNN_purity")
  dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)
  
  files <- list.files(input_dir, pattern="\\.qs$", full.names=TRUE)
  results <- list()
  
  for(f in files){
    message("Processing: ", f)
    seu <- qs::qread(f)
    emb <- Embeddings(seu, "harmony")
    meta <- seu@meta.data
    
    knn <- get.knn(emb, k=k)
    
    purities <- sapply(1:nrow(emb), function(i){
      neighbors <- knn$nn.index[i,]
      batch_same <- mean(meta[[batch_col]][neighbors] == meta[[batch_col]][i])
      cond_same  <- mean(meta[[condition_col]][neighbors] == meta[[condition_col]][i])
      c(batch=batch_same, condition=cond_same)
    })
    purities <- t(purities)
    
    ilocal <- median(purities[,1])
    clocal <- median(purities[,2])
    score  <- clocal / ilocal
    
    method <- basename(f) |> str_remove(".qs")
    results[[method]] <- data.frame(
      method=method,
      iPurity=ilocal,
      cPurity=clocal,
      Score=score
    )
  }
  
  results_df <- bind_rows(results)
  
  # ============================================================
  # SAVE SCORES
  # ============================================================
  write.csv(results_df,
            file.path(save_dir,"knn_purity_scores.csv"),
            row.names=FALSE)
  
  # ============================================================
  # PLOT
  # ============================================================
  
  p <- ggplot(results_df,
              aes(x=iPurity, y=cPurity, label=method)) +
    geom_point(size=4) +
    geom_text(vjust=-1) +
    theme_classic() +
    labs(
      title="Batch integration benchmark",
      x="Batch mixing (iPurity)",
      y="Condition conservation (cPurity)"
    )
  
  ggsave(
    file.path(save_dir,"knn_purity_comparison.png"),
    p,
    width=7,
    height=5
  )
  
  message("KNN-purity computation finished")
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
      c("-l","--label_col"),
      type="character",
      default="condition",
      help="Column in metadata for cell type labels (default: condition)"
      )
    
  )
  
  opt <- parse_args(OptionParser(option_list=option_list))
  
  # ---- Logging ----
  message("\n========== Pipeline Configuration ==========")
  message("Input dir : ", opt$input_dir)
  message("Dataset    : ", opt$dataset)
  message("Batch col  : ", opt$batch_col)
  message("Output dir : ", opt$outdir)
  message("Label col  : ", opt$label_col)
  message("===========================================\n")
  
  # ---- Run QC Mitochondria ----
  run_knn_purity(opt$input_dir, opt$dataset, opt$batch_col, opt$label_col, opt$outdir)
  
}

################
#--- SOURCE ---#
################

if (sys.nframe() == 0) {
  main()
}