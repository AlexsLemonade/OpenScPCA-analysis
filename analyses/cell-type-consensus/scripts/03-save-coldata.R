#!/usr/bin/env Rscript

# This script is used to grab the colData from a SCE object and save it as a TSV file

library(optparse)

option_list <- list(
  make_option(
    opt_str = c("--sce_file"),
    type = "character",
    help = "Path to RDS file containing a processed SingleCellExperiment object from scpca-nf"
  ),
  make_option(
    opt_str = c("--output_file"),
    type = "character",
    help = "Path to file where colData will be saved, must end in `.tsv`"
  )
)

# Parse options
opt <- parse_args(OptionParser(option_list = option_list))

# Set up -----------------------------------------------------------------------

# make sure input files exist
stopifnot(
  "sce file does not exist" = file.exists(opt$sce_file)
)

# load SCE
suppressPackageStartupMessages({
  library(SingleCellExperiment)
})

# Extract colData --------------------------------------------------------------

# read in sce 
sce <- readr::read_rds(opt$sce_file)

# extract ids 
library_id <- metadata(sce)$library_id
# account for multiplexed libraries that have multiple samples 
# for now just combine sample ids into a single string and don't worry about demultiplexing 
sample_id <- metadata(sce)$sample_id |> 
  paste0(collapse = ";")
project_id <- metadata(sce)$project_id

# check if cell line since cell lines don't have any cell type assignments 
# account for having more than one sample and a list of sample types
# all sample types should be the same theoretically 
is_cell_line <- all(metadata(sce)$sample_type == "cell line")

# only create and write table for non-cell line samples
if(!is_cell_line){
  
  # get df with ids, barcodes, and cell type assignments
  celltype_df <- colData(sce) |> 
    as.data.frame() |> 
    dplyr::mutate(
      project_id = project_id,
      sample_id = sample_id,
      library_id = library_id
    ) |> 
    dplyr::select(
      project_id,
      sample_id,
      library_id,
      barcodes, 
      contains("celltype") # get both singler and cellassign with ontology 
    )
  
  # save tsv 
  readr::write_tsv(celltype_df, opt$output_file)
  
}

