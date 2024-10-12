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
library(fs)
library(jpeg)
library(png)

# Parse arguments --------------------------------------------------------------
# set up arguments
option_list <- list(
  make_option(
    opt_str = c("-s", "--sample_id"),
    type = "character",
    default = "SCPCS000179",
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
    type = "character",
    default = "euclidean",
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


# Create directories to save the results of copykat with/without reference using opts$distance
dir.create(file.path(result_dir,  "05_copyKAT", "noref", opts$distance), recursive = TRUE)
dir.create(file.path(result_dir,  "05_copyKAT", "ref", opts$distance), recursive = TRUE)

# define scratch directory for tempory saving the output of copykat
scratch_dir <- file.path(module_base, "scratch", opts$sample_id)
dir.create(scratch_dir, recursive = TRUE)

# path for copykat rds output
name_ref <- glue::glue("05_copykat_",opts$sample_id,"_ref_distance-", opts$distance, ".rds")
name_ref_full <- file.path(result_dir,  "05_copyKAT", "ref", opts$distance, name_ref)
name_no_ref <- glue::glue("05_copykat_",opts$sample_id,"_noref_distance-", opts$distance, ".rds")
name_no_ref_full <- file.path(result_dir,  "05_copyKAT", "noref", opts$distance, name_no_ref)

# path to scratch and final heatmap file to copy over
jpeg_file <- glue::glue(opts$sample_id,"_copykat_heatmap.jpeg")
scratch_jpeg <- file.path(scratch_dir, jpeg_file)
output_jpeg_ref <- file.path(result_dir,  "05_copyKAT", "ref", opts$distance, glue::glue("05_copykat_",opts$sample_id,"_ref_distance-", opts$distance, "_copykat_heatmap.png"))
output_jpeg_noref <- file.path(result_dir,  "05_copyKAT", "noref", opts$distance, glue::glue("05_copykat_",opts$sample_id,"_noref_distance-", opts$distance, "_copykat_heatmap.png"))

# path to scratch and final .txt prediction file to copy over
prediction_file <- glue::glue(opts$sample_id,"_copykat_prediction.txt")
scratch_prediction <- file.path(scratch_dir, prediction_file)
output_prediction_ref <- file.path(result_dir,  "05_copyKAT", "ref", opts$distance, glue::glue("05_copykat_",opts$sample_id,"_ref_distance-", opts$distance, "_copykat_prediction.txt"))
output_prediction_noref <- file.path(result_dir,  "05_copyKAT", "noref", opts$distance, glue::glue("05_copykat_",opts$sample_id,"_noref_distance-", opts$distance, "_copykat_prediction.txt"))

# path to scratch and final .txt CNA file to copy over
CNA_file <- glue::glue(opts$sample_id,"_copykat_CNA_results.txt")
scratch_CNA <- file.path(scratch_dir, CNA_file)
output_CNA_ref <- file.path(result_dir,  "05_copyKAT", "ref", opts$distance, glue::glue("05_copykat_",opts$sample_id,"_ref_distance-", opts$distance, "_copykat_CNA_results.txt"))
output_CNA_noref <- file.path(result_dir,  "05_copyKAT", "noref", opts$distance, glue::glue("05_copykat_",opts$sample_id,"_noref_distance-", opts$distance, "_copykat_CNA_results.txt"))


# change working directory of the script to the scratch directory
# this ensures copykat files get saved to the right location
# there is no option to specify an output directory when running copykat
setwd(scratch_dir)

# Read in data -----------------------------------------------------------------
srat <- readRDS(
  file.path(result_dir,  paste0("02b-fetal_kidney_label-transfer_",  opts$sample_id, ".Rds"))
)

# Extract raw counts -----------------------------------------------------------
exp.rawdata <- GetAssayData(object = srat, assay = "RNA", layer = "counts")

# Extract normal cells ---------------------------------------------------------
normal_cell <- WhichCells(object = srat, expression = fetal_kidney_predicted.compartment %in% c("endothelium", "immune"))

# Run copyKAT without reference ------------------------------------------------

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

# Save copykat output without reference ----------------------------------------

saveRDS(copykat.noref, name_no_ref_full)

img <- readJPEG(scratch_jpeg)
writePNG(img, target = output_jpeg_noref)

fs::file_copy(scratch_prediction, output_prediction_noref, overwrite = TRUE)
fs::file_copy(scratch_CNA, output_CNA_noref, overwrite = TRUE)

# Run copyKAT with reference  --------------------------------------------------

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

# Save copykat output with reference -------------------------------------------

saveRDS(copykat.ref,name_ref_full)

img <- readJPEG(scratch_jpeg)
writePNG(img, target = output_jpeg_ref)

fs::file_copy(scratch_prediction, output_prediction_ref, overwrite = TRUE)
fs::file_copy(scratch_CNA, output_CNA_ref, overwrite = TRUE)
