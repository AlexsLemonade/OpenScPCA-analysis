#!/usr/bin/env Rscript
# This script creates normal references for on specified reference cell groups
# References are created by pooling relevant cell types across all project samples, excluding PDX's.
# We create three references:
# - `all-normal`: All normal cells, which here includes immune, endothelial, epithelial, and adipocytes
# - `immune`: All immune cells
# - `endo`: All endothelial cells
# For more information about how references were determined, see the exploratory notebook in:
# ../../exploratory-notebooks/SCPCP000004/01_nb-consensus-cell-types.Rmd

suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(optparse)
})

# source functions
module_dir <- here::here()
source(
  file.path(module_dir, "scripts", "utils.R")
)

# Define a helper function to split strings
split_input_str <- function(x) {
  strsplit(x, ",")[[1]] |>
    stringr::str_trim()
}

option_list <- list(
  make_option(
    opt_str = "--merged_sce_file",
    type = "character",
    default = "",
    help = "Path to the merged SingleCellExperiment object for the given project"
  ),
  make_option(
    opt_str = "--reference_celltype_group_file",
    type = "character",
    default = "",
    help = "Path to file mapping consensus cell types to reference cell type groups"
  ),
  make_option(
    opt_str = "--reference_groups",
    type = "character",
    default = "",
    help = "Comma-separated list of reference group names to build references for.
      Use the string 'normal' to include a combined set of normal cell types, which can be further specifed with the `normal_reference_groups` argument."
  ),
  make_option(
    opt_str = "--normal_reference_groups",
    type = "character",
    default = "",
    help = "Optional comma-separated list of reference group names to include in the 'normal' reference; only used when 'normal' is provided to the `reference_groups` argument. 
    If not provided, all normal cell types will be included in the 'normal' grouping."
  ),
  make_option(
    opt_str = "--reference_output_dir",
    type = "character",
    default = "",
    help = "Path to directory where normal reference files, named as `{reference group}.rds`,will be saved"
  ),
  make_option(
    opt_str = "--exclude_libraries",
    type = "character",
    default = "",
    help = "Comma-separated list of any libraries which should be excluded from the references"
  )
)

# Parse options ----------------
opts <- parse_args(OptionParser(option_list = option_list))

# check input files and paths
stopifnot(
  "merged_sce_file does not exist" = file.exists(opts$merged_sce_file),
  "reference_celltype_group_file does not exist" = file.exists(opts$reference_celltype_group_file)
)

# read input data -----------
celltype_group_df <- readr::read_tsv(opts$reference_celltype_group_file)
merged_sce <- readRDS(opts$merged_sce_file)

# Ensure is_xenograft and is_cell_line columns are logical
# context: https://github.com/AlexsLemonade/scpcaTools/issues/302
merged_sce$is_xenograft <- as.logical(merged_sce$is_xenograft)
merged_sce$is_cell_line <- as.logical(merged_sce$is_cell_line)

# Define reference groups -------
reference_groups <- split_input_str(opts$reference_groups)
normal_reference_groups <- split_input_str(opts$normal_reference_groups)

# prepare data frame of all possible cells to consider -----------
exclude_libraries <- split_input_str(opts$exclude_libraries)     

cell_df <- colData(merged_sce) |>
  as.data.frame() |>
  # remove PDX and cell lines
  # remove any excluded ids passed in
  dplyr::filter(
    !is_xenograft,
    !is_cell_line, 
    !(library_id %in% exclude_libraries) 
  ) |>
  dplyr::select(cell_id, consensus_annotation = consensus_celltype_annotation) |>
  dplyr::inner_join(celltype_group_df) 

# create and export a reference SCE for each reference group -----
reference_groups |>
  purrr::walk(
    \(ref_group) {
      
      output_file <- file.path(
        opts$reference_output_dir,
        glue::glue("{ref_group}.rds")
      )
      
      # define cell groups to consider
      if (ref_group == "normal") {
        # use just provided normal types, or use all normal types
        if (length(normal_reference_groups) == 0) {
          all_groups <- uniue(cell_df$reference_cell_group)
        } else {
          all_groups <- normal_reference_groups
        }
      } else {
        all_groups <- ref_group
      }
      
      cell_ids <- cell_df |>
        dplyr::filter(reference_cell_group %in% all_groups) |>
        dplyr::pull(cell_id)
      
      # Subset and export
      merged_sce[, cell_ids] |>
        clean_sce() |>
        readr::write_rds(output_file, compress = "gz")
    })
