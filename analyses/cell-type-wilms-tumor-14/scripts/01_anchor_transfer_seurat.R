#!/usr/bin/env Rscript
# parse arguments from command line
args = commandArgs(trailingOnly=TRUE)
if (length(args) == 0) {
  path_repo <- rprojroot::find_root(rprojroot::is_git_root)
  nsample <- 10
} else if (length(args) == 1) {
  nsample <- as.numeric(args[1])
  path_repo <- rprojroot::find_root(rprojroot::is_git_root)
} else if (length(args) == 2) {
  nsample <- as.numeric(args[1])
  path_repo <- args[2]
}else {
  stop("Usage: Rscript 01_anchor_transfer_seurat.R <optional_num_of_test_samples> <opt_path_repo>", call.=FALSE)
}

library(dplyr)
library(Seurat)
library(ggpubr)

path_anal <- file.path(path_repo,"analyses","cell-type-wilms-tumor-14") 
path_proj <- file.path(path_repo,"data","current","SCPCP000014")
path_meta <- file.path(path_proj,"single_cell_metadata.tsv")
meta <- read.table(path_meta, sep = "\t", header = T, stringsAsFactors = F) %>%
  head(n = nsample)

# create output dirs
scratch_out_dir <- file.path(path_anal, "scratch", "01_anchor_transfer_seurat")
dir.create(scratch_out_dir, showWarnings = F, recursive = T)
results_out_dir <- file.path(path_anal, "results", "01_anchor_transfer_seurat")
dir.create(results_out_dir, showWarnings = F, recursive = T)

source(file = file.path(path_anal,"scripts","utils","01_anchor_transfer_seurat_functions.R"))

########### Prepare reference seurat obj ##########
prepare_fetal_atlas(scratch_out_dir = scratch_out_dir, use_exist = T)

########### Run anchor transfer ##########
ref_obj <- SeuratObject::LoadSeuratRds(file.path(scratch_out_dir, "kidneyatlas.rdsSeurat"))

purrr::walk(
  meta$scpca_sample_id,
  \(sample) run_anchorTrans(path_anal = path_anal, scratch_out_dir = scratch_out_dir, results_out_dir = results_out_dir,
                                     ref_obj = ref_obj, sample = sample,
                                     unknown_cutoff = 0.5)
)


