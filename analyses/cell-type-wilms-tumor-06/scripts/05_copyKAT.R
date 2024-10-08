#!/usr/bin/env Rscript

# Run `copyKAT` for one sample with or without a healthy reference
# copyKAT
#
# USAGE:
# Rscript copyKAT.R \
#   --sample_id SCPCS000194
#   --ncore 16
#

library(optparse)
library(Seurat)
library(copykat)

# Parse arguments --------------------------------------------------------------
# set up arguments
option_list <- list(
  make_option(
    opt_str = c("-s", "--sample_id"),
    type = "character",
    default = "SCPCS000205",
    help = "The sample_id of the sample to be used for inference of genomic copy number using copyKAT "
  ),
  make_option(
    opt_str = c("-c", "--n_core"),
    type = "integer",
    default = 16,
    help = "number of cores used to run copyKAT"
  )
  ,
  make_option(
    opt_str = c("-d", "--distance"),
    type = "integer",
    default = 16,
    help = "method used to calculate distance in copyKAT"
  )
)

opts <- parse_args(OptionParser(option_list = option_list))

# paths to data ----------------------------------------------------------------

# The base path for the OpenScPCA repository, found by its (hidden) .git directory
repository_base <- rprojroot::find_root(rprojroot::is_git_root)
# The path to this module
module_base <- file.path(repository_base, "analyses", "cell-type-wilms-tumor-06")
# Path to the result directory
result_dir <- file.path(module_base, "results", opts$sample_id)
# path to output copykat object
copykat_output_obj_noref <- file.path(result_dir,  "copyKAT", "noref", "05_final-copykat.rds")
copykat_output_obj_ref <- file.path(result_dir, "copyKAT", "ref", "05_final-copykat.rds")


# Read in data -----------------------------------------------------------------
srat <- readRDS(
  file.path(result_dir,  paste0("02b-fetal_kidney_label-transfer_",  opts$sample_id, ".Rds"))
)

# Extract raw counts -----------------------------------------------------------
exp.rawdata <- GetAssayData(object = srat, assay = "RNA", layer = "counts")

# Extract normal cells ---------------------------------------------------------
normal_cell <- WhichCells(object = srat, expression = fetal_kidney_predicted.compartment %in% c("endothelium", "immune"))

# Run copyKAT ------------------------------------------------------------------
# change working directory of the script to the result directory
# this ensures copykat files get saved to the right location
# there is no option to specify an output directory when running copykat
dir.create(file.path(result_dir,  "05_copyKAT", "noref", opts$distance), recursive = TRUE)
setwd(file.path(result_dir,  "05_copyKAT", "noref", opts$distance))
copykat.noref <- copykat(rawmat=exp.rawdata, 
                         sam.name=opts$sample_id, 
                         distance=opts$distance, 
                         norm.cell.names="", 
                         genome="hg20",
                         n.cores= opts$n_core, 
                         id.type = "E",
                         plot.genes = FALSE,
                         output.seg = FALSE,
                         KS.cut = 0.05)

dir.create(file.path(result_dir,  "05_copyKAT", "ref", opts$distance), recursive = TRUE)
setwd(file.path(result_dir,  "05_copyKAT", "ref", opts$distance))
copykat.ref <- copykat(rawmat=exp.rawdata, 
                       sam.name=opts$sample_id, 
                       distance=opts$distance, 
                       norm.cell.names=normal_cell, 
                       genome="hg20",
                       n.cores= opts$n_core, 
                       id.type = "E",
                       plot.genes = FALSE,
                       output.seg = FALSE,
                       KS.cut = 0.05
                        )



