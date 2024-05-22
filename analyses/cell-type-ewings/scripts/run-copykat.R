#!/usr/bin/env Rscript

# This script is used to run CopyKAT on a given ScPCA library 
# The output will be found in `outdir`/`results_prefix` and all files will be labeled with the `library_id`

project_root <- here::here()
renv::load(project_root)

library(optparse)
library(copykat)
library(SingleCellExperiment)

option_list <- list(
  make_option(
    opt_str = c("--sce_file"),
    type = "character",
    help = "Path to RDS file containing a processed SingleCellExperiment object from scpca-nf"
  ),
  make_option(
    opt_str = c("--reference_cell_file"),
    type = "character",
    default = NULL,
    help = "Optional file that contains a table of cell barcodes. 
      Any cells with Normal in the reference_cell_class column will be used as a reference for CopyKAT"
  ),
  make_option(
    opt_str = c("--results_dir"),
    type = "character",
    default = NULL,
    help = "Full path to folder to save final CopyKAT results"
  ),
  make_option(
    opt_str = c("--scratch_dir"),
    type = "character",
    default = file.path(project_root, "scratch", "copykat"),
    help = "path to folder where all intermediate CopyKAT results live"
  ),
  make_option(
    opt_str = c("-t", "--threads"),
    type = "integer",
    default = 1,
    help = "Number of multiprocessing threads to use."
  )
)

# Parse options
opt <- parse_args(OptionParser(option_list = option_list))

# Set up -----------------------------------------------------------------------

# make sure path to sce file exists
stopifnot("sce_file does not exist" = file.exists(opt$sce_file))

# check that results dir was provided 
stopifnot("Must provide a --results_dir to save results" = !is.null(opt$results_dir))

# read in sce file
sce <- readr::read_rds(opt$sce_file)

# get library id to use as prefix for all output file names 
library_id <- metadata(sce)$library_id

# contstruct and create output folder if not already present 
fs::dir_create(opt$results_dir)

# define scratch directory 
scratch_dir <- file.path(opt$scratch_dir, library_id)

# path to output copykat object 
copykat_output_obj <- file.path(opt$results_dir, glue::glue("{library_id}_final-copykat.rds"))

# path to scratch and final heatmap file to copy over 
jpeg_file <- glue::glue("{library_id}_copykat_heatmap.jpeg")
scratch_jpeg <- file.path(scratch_dir, jpeg_file)
output_jpeg <- file.path(opt$results_dir, jpeg_file)

# change working directory of the script to the scratch directory
# this ensures copykat files get saved to the right location
# there is no option to specify an output directory when running copykat
setwd(scratch_dir)

# Run CopyKAT ------------------------------------------------------------------

# if normal cells exist then read in and set option for running copyKAT
if(!is.null(opt$reference_cell_file)){
  
  # make sure normal cells file exists 
  stopifnot("reference file does not exist" = file.exists(opt$reference_cell_file))
  normal_cells <- readr::read_tsv(opt$reference_cell_file) |> 
    dplyr::filter(reference_cell_class == "Normal") |> 
    dplyr::pull(barcodes)
  
  # check that all normal cells are in colnames of sce 
  stopifnot("All barcodes in the normal_cells file are not found in the sce object" = 
              all(normal_cells %in% colnames(sce)))
  
} else {
  # otherwise set normal cells to default 
  normal_cells <- ""
}

# Run copyKat without normal cells 
copykat_result <- copykat(
  rawmat = as.matrix(counts(sce)),
  id.type = "E",
  sam.name = library_id,
  norm.cell.names = normal_cells,
  plot.genes = FALSE,
  output.seg = FALSE,
  n.cores = opt$threads
)

# Save outputs 
readr::write_rds(copykat_result, copykat_output_obj)

# copy over png file 
fs::file_copy(scratch_jpeg, output_jpeg, overwrite = TRUE)

