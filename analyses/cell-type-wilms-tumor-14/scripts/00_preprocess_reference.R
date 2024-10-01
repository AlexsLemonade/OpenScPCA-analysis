#!/usr/bin/env Rscript
# parse arguments from command line

library(optparse)

option_list <- list(
  make_option(
    opt_str = c("--in_fetal_atlas"),
    type = "character",
    default = NULL,
    help = "Path to kidney atlas (fetal) in h5ad"
  ),
  make_option(
    opt_str = c("--out_fetal_atlas"),
    type = "character",
    default = NULL,
    help = "Path to converted seurat object for kidney atlas (fetal)"
  )
)

# Parse options
opt <- parse_args(OptionParser(option_list = option_list))

# make sure all input files exist
stopifnot(
  "in_fetal_atlas does not exist" = file.exists(opt$in_fetal_atlas)
)

# set up path to this analysis module
path_repo <- rprojroot::find_root(rprojroot::is_git_root)
path_anal <- file.path(path_repo,"analyses","cell-type-wilms-tumor-14") 

source(file = file.path(path_anal,"scripts","utils","00_preprocess_reference_functions.R"))

########### Prepare reference seurat obj ##########
prepare_fetal_atlas(in_fetal_atlas = opt$in_fetal_atlas,
                    out_fetal_atlas = opt$out_fetal_atlas,
                    use_exist = T)

