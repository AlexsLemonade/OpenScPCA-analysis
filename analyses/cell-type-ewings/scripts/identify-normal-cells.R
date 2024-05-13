#!/usr/bin/env Rscript

# This script is used to generate a list of normal cells for a given ScPCA libray 
# The normal cell barcodes can then be used as input for any downstream analysis that require a 
# list of normal cells (e.g., CopyKAT)

project_root <- here::here()
renv::load(project_root)

library(optparse)
library(SingleCellExperiment)

option_list <- list(
  make_option(
    opt_str = c("--sce_file"),
    type = "character",
    help = "Path to RDS file containing a processed SingleCellExperiment object from scpca-nf"
  ),
  make_option(
    opt_str = c("--singler_normal_cells"),
    type = "character",
    default = NULL,
    help = "Comma separated list of normal cell types to pull from SingleR annotations.
      If both singler_normal_cells and cellassign_normal_cells are used,
      the intersection of cells from both will be saved as normal cells."
  ),
  make_option(
    opt_str = c("--cellassign_normal_cells"),
    type = "character",
    default = NULL,
    help = "Comma separated list of normal cell types to pull from CellAssign annotations. 
      If both singler_normal_cells and cellassign_normal_cells are used,
      the intersection of cells from both will be saved as normal cells."
  ),
  make_option(
    opt_str = c("--output_dir"),
    type = "character",
    default = file.path(project_root, "references", "normal_cell_lists"),
    help = "path to folder to save normal cell file"
  ),
  make_option(
    opt_str = c("--output_filename"),
    type = "character",
    help = "filename to use for saved list of normal cells. Will be saved inside `--output_dir/library_id`.
      Must end in `.txt`"
  )
)

# Parse options
opt <- parse_args(OptionParser(option_list = option_list))

# Set up -----------------------------------------------------------------------

# make sure path to sce file exists
stopifnot("sce_file does not exist" = file.exists(opt$sce_file))

# make sure output filename has the correct extension 
stopifnot("output_filename must have .txt extension" = stringr::str_detect(opt$output_filename, ".txt$"))

# read in sce file
sce <- readr::read_rds(opt$sce_file)

# get library id 
library_id <- metadata(sce)$library_id

# create output directory if not already present
output_dir <- file.path(opt$output_dir, library_id)
fs::dir_create(output_dir)
normal_cell_file <- file.path(output_dir, opt$output_filename)


# Create normal cell vector ----------------------------------------------------

coldata_df <- colData(sce) |>
  as.data.frame()

# get list of singleR cell types
if(!is.null(opt$singler_normal_cells)){
  singler_normal_cells <- stringr::str_split(opt$singler_normal_cells, pattern = ",") |>
    unlist()
  if(!all(singler_normal_cells %in% unique(sce$singler_celltype_annotation))){
    stop("--singler_normal_cells must be present in the singler_celltype_annotation column of the SCE object.")
  }
  
  singler_cells <- coldata_df |>
    dplyr::filter(singler_celltype_annotation %in% singler_normal_cells) |>
    dplyr::pull(barcodes)
  
} else {
  singler_cells <- c()
}

# get list of CelllAssign cell types
if(!is.null(opt$cellassign_normal_cells)){
  cellassign_normal_cells <- stringr::str_split(opt$cellassign_normal_cells, pattern = ",") |>
    unlist()
  if(!all(cellassign_normal_cells %in% unique(sce$cellassign_celltype_annotation))){
    stop("--cellassign_normal_cells must be present in the cellassign_celltype_annotation column of the SCE object.")
  }
  
  cellassign_cells <- coldata_df |>
    dplyr::filter(cellassign_celltype_annotation %in% cellassign_normal_cells) |>
    dplyr::pull(barcodes)
} else {
  cellassign_cells <- c()
}

# create vector of cells that should be used as normal reference
all_normal_cells <- intersect(singler_cells, cellassign_cells)

# export list of normal cells 
writeLines(all_normal_cells, normal_cell_file)


