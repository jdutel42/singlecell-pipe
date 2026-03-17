#!/usr/bin/env Rscript

# ============================================================
# Harmony Grid Search Integration
# ============================================================

suppressPackageStartupMessages({
  library(optparse)
  library(qs)
  library(Seurat)
  library(harmony)
})


##################
#--- FUNCTION ---#
##################

run_harmony <- function(seu_file_path, dataset_name, batch_col, outdir, theta, sigma, dims) {
  
  ########################
  # LOAD OBJECT
  ########################
  
  message("Loading Seurat object...")
  
  seu <- qs::qread(seu_file_path)
  
  if(!inherits(seu,"Seurat")){
    stop("Input must be a Seurat object")
  }
  
  # ------------------------------------------------------------
  # 2. VALIDATION
  # ------------------------------------------------------------
  
  if(!(batch_col %in% colnames(seu@meta.data))){
    stop("Batch column not found in metadata")
  }
  
  if(!"pca" %in% names(seu@reductions)){
    stop("PCA reduction not found. Run PCA before grid search.")
  }
  
  ########################
  # GRID PARAMETERS
  ########################
  
  theta_grid <- as.numeric(strsplit(theta,",")[[1]])
  sigma_grid <- as.numeric(strsplit(sigma,",")[[1]])
  
  ########################
  # OUTPUT DIR
  ########################
  
  save_dir <- file.path(
    outdir,
    dataset_name,
    "Integration",
    "Harmony_gridsearch"
  )
  
  dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)
  
  ########################
  # RUN GRIDSEARCH
  ########################
  
  DefaultAssay(seu) <- "SCT"
  
  for(theta in theta_grid){
    
    for(sigma in sigma_grid){
      
      message(
        paste0(
          "Running Harmony with theta=",theta,
          " sigma=",sigma
        )
      )
      
      seu_tmp <- RunHarmony(
        seu,
        group.by.vars = batch_col,
        reduction = "pca",
        max.iter = 20,
        dims.use = 1:dims,
        theta = theta,
        sigma = sigma,
        reduction.save = "harmony",
        verbose = TRUE
      )
      
      filename <- paste0(
        "harmony_theta",theta,
        "_sigma",sigma,
        ".qs"
      )
      
      qs::qsave(
        seu_tmp,
        file.path(save_dir,filename),
        preset = "fast"
      )
      
    }
    
  }

  message("Harmony grid search completed.")
  
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
      c("-i","--input"),
      type="character",
      help="Input Seurat object (.qs)"
    ),
    
    make_option(
      c("-b","--batch_col"),
      type="character",
      default="emulsion",
      help="Batch variable in metadata"
    ),
    
    make_option(
      c("-o","--outdir"),
      type="character",
      default="03_result",
      help="Output directory"
    ),
    
    make_option(
      c("-d","--dataset"),
      type="character",
      help="Dataset name"
    ),
    
    make_option(
      c("-t","--theta"),
      type="character",
      default="1,2,3,4",
      help="Theta grid (comma separated)"
    ),
    
    make_option(
      c("-s", "--sigma"),
      type="character",
      default="0.05,0.1,0.2",
      help="Sigma grid (comma separated)"
    ),
    
    make_option(
      c("-p", "--dims"),
      type="integer",
      default=30,
      help="Number of PCA dimensions"
    )
    
  )
  
  opt <- parse_args(OptionParser(option_list=option_list))


  
  # ---- Logging ----
  message("\n========== Pipeline Configuration ==========")
  message("Input file : ", opt$input)
  message("Dataset    : ", opt$dataset)
  message("Batch col  : ", opt$batch_col)
  message("Output dir : ", opt$outdir)
  message("Theta grid : ", opt$theta)
  message("Sigma grid : ", opt$sigma)
  message("PCA dims   : ", opt$dims)
  message("===========================================\n")
  
  # ---- Run QC Mitochondria ----
  run_harmony(opt$input, opt$dataset, opt$batch_col, opt$outdir, opt$theta, opt$sigma, opt$dims)
}



################
#--- SOURCE ---#
################

if (sys.nframe() == 0) {
  main()
}