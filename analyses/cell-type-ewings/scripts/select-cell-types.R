#!/usr/bin/env Rscript

# This script is used to generate a table of reference cells for a given ScPCA libray that are annotated by 
# SingleR and/or CellAssign to have a specific label
# The output is a table that contains columns for barcodes, reference cell class, 
# and then columns for each of the original cell type annotations for that cell 
# The `reference_cell_class` column will contain either tumor or normal 

project_root <- here::here()
renv::load(project_root)

suppressPackageStartupMessages({
  library(optparse)
  library(SingleCellExperiment)
})

option_list <- list(
  make_option(
    opt_str = c("--sce_file"),
    type = "character",
    help = "Path to RDS file containing a processed SingleCellExperiment object from scpca-nf"
  ),
  make_option(
    opt_str = c("--normal_cells"),
    type = "character",
    default = "",
    help = "Comma separated list of cell types to pull from annotations present and save as normal cells."
  ),
  make_option(
    opt_str = c("--tumor_cells"),
    type = "character",
    default = NULL,
    help = "Comma separated list of cell types to pull from annotations present and save as tumor cells."
  ),
  make_option(
    opt_str = c("--output_filename"),
    type = "character",
    help = "Full path to file to store saved table of cell barcodes."
  )
)

# Parse options
opt <- parse_args(OptionParser(option_list = option_list))

# Set up -----------------------------------------------------------------------

stopifnot(
  # make sure path to sce file exists
  "sce_file does not exist" = file.exists(opt$sce_file)),

  # check that either normal or tumor cell types are provided 
  "Either normal_cells or tumor_cells must be provided" = 
    !is.null(opt$normal_cells) || !is.null(opt$tumor_cells)
)       

# read in sce file
sce <- readr::read_rds(opt$sce_file)

# create output directory if not already present
output_dir <- dir(opt$output_filename)
fs::dir_create(output_dir)

# Create table of reference cells ----------------------------------------------

# possible annotation columns 
annotation_columns <- c(
  "singler_celltype_annotation",
  "cellassign_celltype_annotation"
)

# get dataframe with barcodes and associated annotation 
coldata_df <- colData(sce) |>
  as.data.frame() |> 
  dplyr::select(barcodes, all_of(annotation_columns)) |> 
  tidyr::pivot_longer(all_of(annotation_columns),
                      names_to = "celltype_method",
                      values_to = "celltype")

# get list of normal and tumor cell types to save as refernences 
if(!is.null(opt$normal_cells)){
  normal_cell_types <- stringr::str_split_1(opt$normal_cells, pattern = ",")
} else {
  normal_cell_types <- ""
}

if(!is.null(opt$tumor_cells)){
  tumor_cell_types <- stringr::str_split(opt$tumor_cells, pattern = ",") |>
    unlist()
} else {
  tumor_cell_types <- ""
}

# label cells as either tumor or normal based on which cell type was annotated 
ref_table_df <- coldata_df |> 
  dplyr::mutate(reference_cell_class = dplyr::case_when(
    celltype %in% normal_cell_types ~ "Normal",
    celltype %in% tumor_cell_types ~ "Tumor",
    .default = "NA"
  )) |> 
  # remove any that weren't tumor or normal 
  dplyr::filter(reference_cell_class != "NA") |> 
  # get one row per barcode and one column for each annotation type 
  tidyr::pivot_wider(values_from = celltype,
                     names_from = celltype_method)

# export table of reference cells 
readr::write_tsv(ref_table_df, opt$output_filename)


