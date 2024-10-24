#!/usr/bin/env Rscript
# parse arguments from command line
library(optparse)
library(dplyr)

option_list <- list(
  make_option(
    opt_str = c("--reference"),
    type = "character",
    default = NULL,
    help = "Path to converted seurat object for reference database"
  ),
  make_option(
    opt_str = c("--metadata"),
    type = "character",
    default = NULL,
    help = "Path to cohort metadata"
  ),
  make_option(
    opt_str = c("--libraries"),
    type = "character",
    default = "all",
    help = "list of libraries to run, comma separated. Default is processing all libraries in this cohort."
  ),
  make_option(
    opt_str = c("--testing"),
    type = "logical",
    default = FALSE,
    action = "store_true",
    help = "Use this flag when running on test data"
  ),
  make_option(
    opt_str = c("--run_LogNormalize"),
    type = "logical",
    default = FALSE,
    action = "store_true",
    help = "Use this flag to generate results based on normalization method LogNormalize"
  ),
  make_option(
    opt_str = c("--run_SCT"),
    type = "logical",
    default = FALSE,
    action = "store_true",
    help = "Use this flag to generate results based on normalization method SCT"
  )
)

# Parse options
opt <- parse_args(OptionParser(option_list = option_list))
running_ci <- opt$testing
run_LogNormalize <- opt$run_LogNormalize
run_SCT <- opt$run_SCT

# make sure all input files exist
stopifnot(
  "reference does not exist" = file.exists(opt$reference),
  "metadata does not exist" = file.exists(opt$metadata)
)





path_repo <- rprojroot::find_root(rprojroot::is_git_root)
path_anal <- file.path(path_repo,"analyses","cell-type-wilms-tumor-14") 
# path_meta <- file.path(path_repo,"data","current","SCPCP000014","single_cell_metadata.tsv") # keep for debug
path_meta <- file.path(opt$metadata)
meta <- read.table(path_meta, sep = "\t", header = TRUE, stringsAsFactors = FALSE) 

if (opt$libraries == "all") {
  libraries <- meta$scpca_library_id
} else {
  libraries <- stringr::str_split_1(opt$libraries, ",")
}

# Use the default value for k.weight when working with real data
if (running_ci) {
  k_weight <- 15  # Smaller number when working with simulated data (n = 100 cells)
} else {
  k_weight <- 50  # (Current) Seurat default
}


# create output dirs
scratch_out_dir <- file.path(path_anal, "scratch", "01_anchor_transfer_seurat")
dir.create(scratch_out_dir, showWarnings = FALSE, recursive = TRUE)
results_out_dir <- file.path(path_anal, "results", "01_anchor_transfer_seurat")
dir.create(results_out_dir, showWarnings = FALSE, recursive = TRUE)
plots_out_dir <- file.path(path_anal, "plots", "01_anchor_transfer_seurat")
dir.create(results_out_dir, showWarnings = FALSE, recursive = TRUE)

source(file = file.path(path_anal,"scripts","utils","01_anchor_transfer_seurat_functions.R"))


########### Run anchor transfer ##########
# library = "SCPCL000850"
#ref_obj <- SeuratObject::LoadSeuratRds("scratch/00_preprocess_reference/kidneyatlas.rdsSeurat")
ref_obj <- SeuratObject::LoadSeuratRds(opt$reference)

if (run_LogNormalize) {
  used_assay <- "RNA"
  
  # make sure output folder exists
  dir.create(file.path(scratch_out_dir,used_assay), showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(results_out_dir,used_assay), showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(plots_out_dir,used_assay), showWarnings = FALSE, recursive = TRUE)
  
  purrr::walk(
    libraries,
    \(library) run_anchorTrans(path_anal = path_anal, 
                               scratch_out_dir = file.path(scratch_out_dir,used_assay), 
                               results_out_dir = file.path(results_out_dir,used_assay),
                               plots_out_dir = file.path(plots_out_dir,used_assay),
                               ref_obj = ref_obj, 
                               library = library, 
                               level = "compartment",
                               k_weight = k_weight,
                               unknown_cutoff = 0.5, ndims = 15,
                               obj_assay = used_assay)
  )
  
  purrr::walk(
    libraries,
    \(library) run_anchorTrans(path_anal = path_anal, 
                               scratch_out_dir = file.path(scratch_out_dir,used_assay), 
                               results_out_dir = file.path(results_out_dir,used_assay),
                               plots_out_dir = file.path(plots_out_dir,used_assay),
                               ref_obj = ref_obj, 
                               library = library, 
                               level = "celltype",
                               k_weight = k_weight,
                               unknown_cutoff = 0.5, ndims = 15,
                               obj_assay = used_assay)
  )
}

if (run_SCT) {
  used_assay <- "SCT"
  
  # make sure output folder exists
  dir.create(file.path(scratch_out_dir,used_assay), showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(results_out_dir,used_assay), showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(plots_out_dir,used_assay), showWarnings = FALSE, recursive = TRUE)
  
  purrr::walk(
    libraries,
    \(library) run_anchorTrans(path_anal = path_anal, 
                               scratch_out_dir = file.path(scratch_out_dir,used_assay), 
                               results_out_dir = file.path(results_out_dir,used_assay),
                               plots_out_dir = file.path(plots_out_dir,used_assay),
                               ref_obj = ref_obj, 
                               library = library, 
                               level = "compartment",
                               k_weight = k_weight,
                               unknown_cutoff = 0.5, ndims = 15,
                               obj_assay = used_assay)
  )
  
  purrr::walk(
    libraries,
    \(library) run_anchorTrans(path_anal = path_anal, 
                               scratch_out_dir = file.path(scratch_out_dir,used_assay), 
                               results_out_dir = file.path(results_out_dir,used_assay),
                               plots_out_dir = file.path(plots_out_dir,used_assay),
                               ref_obj = ref_obj, 
                               library = library, 
                               level = "celltype",
                               k_weight = k_weight,
                               unknown_cutoff = 0.5, ndims = 15,
                               obj_assay = used_assay)
  )
}
