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
    help = "list of libraries to run (starts with ), seperated by comma. Default is processing all libraries in this cohort."
  ),
  make_option(
    opt_str = c("--testing"),
    type = "logical",
    default = FALSE,
    action = "store_true",
    help = "Use this flag when running on test data"
  )
)

# Parse options
opt <- parse_args(OptionParser(option_list = option_list))
running_ci <- opt$testing

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
#ref_obj <- SeuratObject::LoadSeuratRds("scratch/00_preprocess_reference/kidneyatlas.rdsSeurat")
ref_obj <- SeuratObject::LoadSeuratRds(opt$reference)

purrr::walk(
  libraries,
  \(library) run_anchorTrans(path_anal = path_anal, 
                            scratch_out_dir = scratch_out_dir, 
                            results_out_dir = results_out_dir,
                            plots_out_dir = plots_out_dir,
                            ref_obj = ref_obj, 
                            library = library, 
                            level = "compartment",
                            k_weight = k_weight,
                            unknown_cutoff = 0.5, ndims = 20)
)

purrr::walk(
  libraries,
  \(library) run_anchorTrans(path_anal = path_anal, 
                            scratch_out_dir = scratch_out_dir, 
                            results_out_dir = results_out_dir,
                            plots_out_dir = plots_out_dir,
                            ref_obj = ref_obj, 
                            library = library, 
                            level = "celltype",
                            k_weight = k_weight,
                            unknown_cutoff = 0.5, ndims = 20)
)

