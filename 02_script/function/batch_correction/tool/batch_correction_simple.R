#!/usr/bin/env Rscript

# ============================================================
# Harmony Integration + LISI evaluation - scRNA-seq
# ============================================================

#####################
#--- DESCRIPTION ---#
#####################

############################################################
# Harmony Batch Correction Pipeline
#
# Author : Jordan DUTEL
#
# Description :
# Performs batch correction using Harmony and evaluates
# integration quality using LISI scores before and after
# correction.
#
# Outputs :
# - Harmony embedding
# - Convergence plot
# - PCA / Harmony DimPlots
# - LISI score visualizations
#
############################################################

#############
#--- LIB ---#
#############

suppressPackageStartupMessages({
  library(optparse)
  library(qs)
  library(Seurat)
  library(harmony)
  library(lisi)
  library(ggplot2)
  library(dplyr)
})

##################
#--- FUNCTION ---#
##################

run_harmony <- function(seu_file,
                        dataset,
                        batch_col,
                        outdir,
                        dims_use = 11,
                        theta = 3,
                        sigma = 0.2,
                        save_seu = TRUE){
  
  message("Loading Seurat object...")
  seu <- qs::qread(seu_file)
  
  if(!inherits(seu,"Seurat")){
    stop("Input object must be Seurat")
  }
  
  dir.create(
    file.path(outdir,dataset,"Harmony"),
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  save_dir <- file.path(outdir,dataset,"Harmony")
  
  ###########################################################
  # LISI BEFORE CORRECTION
  ###########################################################
  
  message("Computing LISI before correction...")
  
  emb_raw <- Embeddings(seu,"pca")
  meta <- seu@meta.data
  
  lisi_raw <- compute_lisi(
    emb_raw,
    meta,
    c(batch_col,"condition","source")
  )
  
  seu$emulsion_LISI_before <- lisi_raw[[batch_col]]
  seu$condition_LISI_before <- lisi_raw$condition
  seu$source_LISI_before <- lisi_raw$source
  
  
  p_before <- FeaturePlot(
    seu,
    features=c("emulsion_LISI_before",
               "condition_LISI_before",
               "source_LISI_before"),
    cols=c("blue","yellow","red"),
    min.cutoff="q05",
    max.cutoff="q95"
  )
  
  ggsave(
    file.path(save_dir,
              "LISI_scores_before_correction.png"),
    plot=p_before,
    width=12,
    height=6
  )
  
  ###########################################################
  # RUN HARMONY
  ###########################################################
  
  message("Running Harmony...")
  
  DefaultAssay(seu) <- "SCT"
  
  seu <- RunHarmony(
    seu,
    group.by.vars=batch_col,
    reduction="pca",
    dims.use=1:dims_use,
    theta=theta,
    sigma=sigma,
    plot_convergence=TRUE,
    reduction.save="harmony",
    verbose=TRUE
  )
  
  p_converge <- last_plot()
  
  ggsave(
    file.path(save_dir,"Harmony_convergence.png"),
    plot=p_converge,
    width=12,
    height=6
  )
  
  ###########################################################
  # ELBOW
  ###########################################################
  
  p_elbow <- ElbowPlot(
    seu,
    ndims=50,
    reduction="harmony"
  ) + ggtitle("ElbowPlot after Harmony")
  
  ggsave(
    file.path(save_dir,"Harmony_elbowplot.png"),
    plot=p_elbow,
    width=10,
    height=6
  )
  
  ###########################################################
  # PCA / Harmony PLOT
  ###########################################################
  
  p_dim <- DimPlot(
    seu,
    reduction="harmony",
    group.by=batch_col
  )
  
  ggsave(
    file.path(save_dir,"Harmony_dimplot_batch.png"),
    plot=p_dim,
    width=10,
    height=6
  )
  
  ###########################################################
  # LISI AFTER CORRECTION
  ###########################################################
  
  message("Computing LISI after correction...")
  
  emb_int <- Embeddings(seu,"harmony")
  
  lisi_int <- compute_lisi(
    emb_int,
    meta,
    c(batch_col,"condition","source")
  )
  
  seu$emulsion_LISI_after <- lisi_int[[batch_col]]
  seu$condition_LISI_after <- lisi_int$condition
  seu$source_LISI_after <- lisi_int$source
  
  
  p_after <- FeaturePlot(
    seu,
    features=c("emulsion_LISI_after",
               "condition_LISI_after",
               "source_LISI_after"),
    cols=c("blue","yellow","red"),
    min.cutoff="q05",
    max.cutoff="q95"
  )
  
  ggsave(
    file.path(save_dir,
              "LISI_scores_after_correction.png"),
    plot=p_after,
    width=12,
    height=6
  )
  
  
  ###########################################################
  # LISI COMPARISON
  ###########################################################
  
  df <- rbind(
    data.frame(score=lisi_raw[[batch_col]],
               type="Raw",
               var="Batch"),
    
    data.frame(score=lisi_int[[batch_col]],
               type="Integrated",
               var="Batch"),
    
    data.frame(score=lisi_raw$condition,
               type="Raw",
               var="Condition"),
    
    data.frame(score=lisi_int$condition,
               type="Integrated",
               var="Condition")
  )
  
  p_comp <- ggplot(df,aes(type,score,fill=type))+
    geom_violin(trim=TRUE)+
    facet_wrap(~var,scales="free_y")+
    theme_classic()+
    NoLegend()+
    ggtitle("LISI comparison before vs after Harmony")
  
  ggsave(
    file.path(save_dir,
              "LISI_comparison.png"),
    plot=p_comp,
    width=12,
    height=6
  )
  
  ###########################################################
  # SAVE SEURAT
  ###########################################################
  
  if(save_seu){
    message("Saving updated Seurat object...")
    qs::qsave(seu,seu_file)
  }
  
  message("Harmony pipeline completed successfully")
}

##############
#--- MAIN ---#
##############

main <- function(){
  
  option_list <- list(
    
    make_option(
      c("-i","--input"),
      type="character",
      help="Seurat object (.qs)"
    ),
    
    make_option(
      c("-d","--dataset"),
      type="character",
      help="Dataset name"
    ),
    
    make_option(
      c("-b","--batch"),
      type="character",
      default="emulsion",
      help="Batch column"
    ),
    
    make_option(
      c("-o","--outdir"),
      type="character",
      default="03_result"
    ),
    
    make_option(
      c("-s","--save"),
      type="logical",
      default=TRUE
    )
  )
  
  parser <- OptionParser(option_list=option_list)
  opt <- parse_args(parser)
  
  run_harmony(
    opt$input,
    opt$dataset,
    opt$batch,
    opt$outdir,
    save_seu=opt$save
  )
  
}

if(sys.nframe()==0){
  main()
}