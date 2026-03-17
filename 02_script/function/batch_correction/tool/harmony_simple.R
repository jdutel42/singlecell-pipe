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

################
#--- COLOR ----#
################

pal <- c(
  "#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00",
  "#A65628","#F781BF","#999999","#66C2A5","#FC8D62",
  "#8DA0CB","#E78AC3","#A6D854","#FFD92F","#E5C494",
  "#1B9E77","#D95F02","#7570B3","#E7298A","#66A61E",
  "#E6AB02","#A6761D","#666666","#8DD3C7","#BEBADA",
  "#FB8072","#80B1D3","#FDB462","#B3DE69","#FCCDE5",
  "#D9D9D9","#BC80BD","#CCEBC5","#FFED6F","#1F78B4",
  "#33A02C","#FB9A99","#E31A1C","#FDBF6F","#CAB2D6",
  "#6A3D9A","#FFFF99","#B15928","#8E0152","#DE77AE",
  "#7FBC41","#4EB265","#4393C3","#2166AC","#B2182B",
  "#F4A582","#92C5DE","#FDDBC7","#D1E5F0","#542788"
)

##################
#--- FUNCTION ---#
##################

run_harmony <- function(seu_file,
                        dataset,
                        batch_col,
                        dims_use = 11,
                        theta = 3,
                        sigma = 0.05,
                        outdir,
                        save_seu = TRUE){
  
  message("Loading Seurat object...")
  seu <- qs::qread(seu_file)
  
  if(!inherits(seu,"Seurat")){
    stop("Input object must be Seurat")
  }
  
  save_dir <- file.path(outdir,dataset,"Integration","Final","Harmony")
  
  dir.create(
    file.path(save_dir),
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  
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
  
  
  # p1_before <- FeaturePlot(
  #   seu,
  #   features=c("emulsion_LISI_before",
  #              "condition_LISI_before",
  #              "source_LISI_before"),
  #   cols=c("blue","yellow","red"),
  #   min.cutoff="q05",
  #   max.cutoff="q95"
  # )
  # 
  # ggsave(
  #   file.path(save_dir,
  #             "LISI_scores_before_correction.png"),
  #   plot=p1_before,
  #   width=12,
  #   height=6,
  #   dpi=300
  # )
  # 
  # ggsave(
  #   file.path(save_dir,
  #             "LISI_scores_before_correction.pdf"),
  #   plot=p1_before,
  #   width=12,
  #   height=6
  # )
  # 
  # 
  # p2_before <- ggplot(seu@meta.data,
  #        aes(x = predicted.celltype.l2, y = condition_LISI)) +
  #   geom_boxplot(outlier.size = 0.3) +
  #   coord_flip() +
  #   theme_classic()
  # 
  # ggsave(
  #   filename = file.path(
  #     DIR_RESULT,
  #     "All_Tcells",
  #     "Harmony",
  #     "All_Tcell_RNA_Harmony_LISI_scores_by_celltype_before_correction.png"),
  #   plot = p2_before,
  #   width = 12,
  #   height = 6)
  
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
    file.path(save_dir,"harmony_convergence.png"),
    plot=p_converge,
    width=12,
    height=6,
    dpi=300
  )
  
  ggsave(
    file.path(save_dir,"harmony_convergence.pdf"),
    plot=p_converge,
    width=12,
    height=6
  )
  
  ###########################################################
  # ELBOW
  ###########################################################
  
  DefaultAssay(seu) <- "SCT"
  p_elbow <- ElbowPlot(
    seu,
    ndims=dims_use,
    reduction="harmony"
  ) + ggtitle("ElbowPlot after Harmony")
  
  ggsave(
    file.path(save_dir,"harmony_elbowplot.png"),
    plot=p_elbow,
    width=10,
    height=6,
    dpi=300
  )
  
  ggsave(
    file.path(save_dir,"harmony_elbowplot.pdf"),
    plot=p_elbow,
    width=10,
    height=6
  )
  
  ###########################################################
  # PCA / Harmony PLOT
  ###########################################################
  
  p_pca_harmony <- DimPlot(
    seu,
    reduction = "harmony",
    group.by = batch_col
  ) + 
    ggtitle("PCA DimPlot after Harmony - RNA assay") +
    scale_color_manual(values = pal)
  
  ggsave(
    file.path(save_dir,"harmony_dimplot_batch.png"),
    plot=p_pca_harmony,
    width=10,
    height=6,
    dpi=300
  )
  
  ggsave(
    file.path(save_dir,"harmony_dimplot_batch.pdf"),
    plot=p_pca_harmony,
    width=10,
    height=6
  )
  
  p_pca_harmony_wrap <- DimPlot(
    seu,
    reduction = "harmony",
    group.by = batch_col
  ) +
    ggtitle("PCA DimPlot after Harmony - RNA assay") +
    scale_color_manual(values = pal) +
    facet_wrap(as.formula(paste0("~", batch_col)), ncol = 6) +
    NoLegend()
  
  ggsave(
    file.path(save_dir,"harmony_dimplot_batch_facet.png"),
    plot=p_pca_harmony_wrap,
    width=12,
    height=8,
    dpi=300
  )
  
  ggsave(
    file.path(save_dir,"harmony_dimplot_batch_facet.pdf"),
    plot=p_pca_harmony_wrap,
    width=12,
    height=8
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
  
  
  # p_after <- FeaturePlot(
  #   seu,
  #   features=c("emulsion_LISI_after",
  #              "condition_LISI_after",
  #              "source_LISI_after"),
  #   cols=c("blue","yellow","red"),
  #   min.cutoff="q05",
  #   max.cutoff="q95"
  # )
  # 
  # ggsave(
  #   file.path(save_dir,
  #             "LISI_scores_after_correction.png"),
  #   plot=p_after,
  #   width=12,
  #   height=6
  # )
  
  
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
  
  p_comp_violin <- ggplot(df,aes(type,score,fill=type))+
    geom_violin(trim=TRUE)+
    facet_wrap(~var,scales="free_y")+
    theme_classic()+
    NoLegend()+
    ggtitle("LISI comparison before vs after Harmony") +
    labs(y = "LISI Score (Mixing/Integration)")
  
  ggsave(
    file.path(save_dir,
              "all_tcell_rna_harmony_lisi_scores_comparison_violin.png"),
    plot=p_comp_violin,
    width=12,
    height=6
  )
  
  ggsave(
    filename = file.path(
      save_dir,
      "all_tcell_rna_harmony_lisi_scores_comparison_violin.pdf"),
    plot = p_comp_violin,
    width = 12,
    height = 6)
  
  p_comp_box <- ggplot(df, aes(type, score, fill = type)) +
    geom_boxplot() +
    facet_wrap(~var, scales = "free_y") +
    theme_classic() +
    NoLegend() +
    ggtitle("LISI scores before and after Harmony correction") +
    labs(y = "LISI Score (Mixing/Integration)")
  
  ggsave(
    filename = file.path(
      save_dir,
      "all_tcell_rna_harmony_lisi_scores_comparison_boxplot.png"),
    plot = p_comp_box,
    width = 12,
    height = 6,
    dpi = 300)
  
  ggsave(
    filename = file.path(
      save_dir,
      "all_tcell_rna_harmony_lisi_scores_comparison_boxplot.pdf"),
    plot = p_comp_box,
    width = 12,
    height = 6
  )
  
  ###########################################################
  # SAVE SEURAT
  ###########################################################
  
  # ---- Save annotated Seurat object ----
  message("Save updated seurat object...")
  if (!is.null(save_seu) && save_seu != "") {
    qs::qsave(seu, file = save_seu)
    cat("Annotated Seurat object saved here:", save_seu, "\n")
  }
  message("Seurat object with library complexity QC annotations saved.")
  
  message("Library complexity QC completed successfully")
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
      c("-D","--dims"),
      type="integer",
      default=11,
      help="Number of dimensions to use for Harmony"
    ),
    
    make_option(
      c("-t","--theta"),
      type="numeric",
      default=3,
      help="Harmony theta parameter (diversity clustering penalty)"
    ),
    
     make_option(
      c("-S","--sigma"),
      type="numeric",
      default=0.2,
      help="Harmony sigma parameter (bandwidth for soft k-means clustering)"
    ),
    
    make_option(
      c("-o","--outdir"),
      type="character",
      default="03_result"
    ),
    
    make_option(
      c("-s","--save"),
      type="character",
      default=NULL,
      help="Path to save the Seurat object (qs format)"
    )
  )
  
  parser <- OptionParser(option_list=option_list)
  opt <- parse_args(parser)
  
  # ---- Logging ----
  message("\n========== Pipeline Configuration ==========")
  message("Input file : ", opt$input)
  message("Dataset    : ", opt$dataset)
  message("Batch col  : ", opt$batch_col)
  message("PCA dims   : ", opt$dims)
  message("Theta      : ", opt$theta)
  message("Sigma      : ", opt$sigma)
  message("Output dir : ", opt$outdir)
  message("Save Seurat : ", opt$save)
  message("===========================================\n")
  
  run_harmony(
    opt$input,
    opt$dataset,
    opt$batch,
    opt$dims,
    opt$theta,
    opt$sigma,
    opt$outdir,
    opt$save
  )
  
}

if(sys.nframe()==0){
  main()
}