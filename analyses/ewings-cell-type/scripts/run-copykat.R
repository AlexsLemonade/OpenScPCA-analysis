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
    opt_str = c("--normal_cells"),
    type = "character",
    default = NULL,
    help = "Optional file that contains a list of cell barcodes that correspond to normal cells and will be used as a reference"
  ),
  make_option(
    opt_str = c("--results_dir"),
    type = "character",
    default = NULL,
    help = "Folder within `--copykat_results_prefix/library_id` to save copyKAT results"
  ),
  make_option(
    opt_str = c("--copykat_results_prefix"),
    type = "character",
    default = file.path(project_root, "results", "copykat"),
    help = "path to folder where all copykat results live"
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
stopifnot("Must provide a --results_dir to save results" = is.null(opt$results_dir))

# read in sce file
sce <- readr::read_rds(opt$sce_file)

# get library id to use as prefix for all output file names 
library_id <- metadata(sce)$library_id

# contstruct and create output folder if not already present 
output_dir <- file.path(opt$copykat_results_prefix, library_id, opt$results_dir) 
fs::dir_create(output_dir)

# change working directory of the script to the output directory
# this ensures copykat files get saved to the right location
# there is no option to specify an output directory when running copykat
setwd(output_dir)

# Run CopyKAT ------------------------------------------------------------------

# if normal cells exist then read in and set option for running copyKAT
if(!is.null(opt$normal_cells)){
  
  # make sure normal cells file exists 
  stopifnot("normal_cells file does not exist" = file.exists(opt$normal_cells))
  normal_cells <- readLines(opt$normal_cells)
  
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

