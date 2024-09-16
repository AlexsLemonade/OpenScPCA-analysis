#!/usr/bin/env Rscript
# parse arguments from command line
args = commandArgs(trailingOnly=TRUE)
if (length(args) != 1) {
  stop("Usage: Rscript --vanilla 00_preprocessing_rds.R <path_to_repo>", call.=FALSE)
} else if (length(args) == 1) {
    path_repo <- args[1]

}

library(dplyr)
library(Seurat)

path_anal <- paste0(path_repo,"/analyses/cell-type-wilms-tumor-14") 
path_proj <- paste0(path_repo,"/data/current/SCPCP000014")
path_meta <- paste0(path_proj,"/single_cell_metadata.tsv")
meta <- read.table(path_meta, sep = "\t", header = T, stringsAsFactors = F)
db_proj <- paste0(path_repo,"/data/current/results/doublet-detection/SCPCP000014")

########### per sample rds preprocessing ##########
source(file = paste0(path_anal,"/scripts/utils/00_preprocessing_rds_functions.R"))

for (i in 1:nrow(meta)) {
  sample <- meta[i,2]
  library <- meta[i,3]
  process_sce(sample, library, path_proj = path_proj, db_proj = db_proj, path_anal = path_anal)
}