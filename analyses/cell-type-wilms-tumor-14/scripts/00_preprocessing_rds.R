#!/usr/bin/env Rscript
# parse arguments from command line
args = commandArgs(trailingOnly=TRUE)
if (length(args) == 0) {
  path_repo <- rprojroot::find_root(rprojroot::is_git_root)
} else if (length(args) == 1) {
  path_repo <- args[1]
} else {
  stop("Usage: Rscript 00_preprocessing_rds.R <optional_path_to_repo>", call.=FALSE)
}

library(dplyr)
library(Seurat)

path_anal <- file.path(path_repo,"analyses","cell-type-wilms-tumor-14") 
path_proj <- file.path(path_repo,"data","current","SCPCP000014")
path_meta <- file.path(path_proj,"single_cell_metadata.tsv")
meta <- read.table(path_meta, sep = "\t", header = T, stringsAsFactors = F)
db_proj <- file.path(path_repo,"data","current","results","doublet-detection","SCPCP000014")

########### per sample rds preprocessing ##########
source(file = file.path(path_anal,"scripts","utils","00_preprocessing_rds_functions.R"))

purrr::walk2(
  meta$scpca_sample_id,
  meta$scpca_library_id,
  \(sample, library) process_sce(sample, library, path_proj = path_proj, db_proj = db_proj, path_anal = path_anal, get_logcounts = TRUE)
)