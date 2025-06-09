#!/usr/bin/env Rscript
# This script creates normal references for use with Neuroblastoma (SCPCP000004) samples.
# References are created by pooling relevant cell types across all project samples, excluding PDX's.
# We create three references:
# - `all-normal`: All normal cells
# - `immune`: All immune cells
# - `endo`: All endothelial immune cells
# For more information about how references were determined, see the exploratory notebook in:
# ../../exploratory-notebooks/SCPCP000004/01_nb-consensus-cell-types.Rmd

suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(optparse)
})

# source functions
source(here::here("scripts", "build-normal-reference", "utils.R"))

option_list <- list(
  make_option(
    opt_str = "--merged_sce_file",
    type = "character",
    default = "",
    help = "Path to the merged SingleCellExperiment object"
  ),
  make_option(
    opt_str = "--celltype_dir",
    type = "character",
    default = "",
    help = "Path to consensus cell type annotations for project SCPCP000004"
  ),
  make_option(
    opt_str = "--metadata_file",
    type = "character",
    default = "",
    help = "Path to single_cell_metadata.tsv file for project SCPCP000004"
  ),
  make_option(
    opt_str = "--celltype_group_file",
    type = "character",
    default = "",
    help = "Path to file mapping consensus cell types to reference cell type groups"
  ),
  make_option(
    opt_str = "--reference_normal_rds",
    type = "character",
    help = "Path to output RDS file to save an SCE file to use as a normal reference with all normal cells"
  ),
  make_option(
    opt_str = "--reference_immune_rds",
    type = "character",
    help = "Path to output RDS file to save an SCE file to use as a normal reference with all immune cells"
  ),
  make_option(
    opt_str = "--reference_endo_rds",
    type = "character",
    help = "Path to output RDS file to save an SCE file to use as a normal reference with all endothelial cells"
  )
)



# Parse options ----------------
opts <- parse_args(OptionParser(option_list = option_list))

# check input files and paths
stopifnot(
  "merged_sce_file does not exist" = file.exists(opts$merged_sce_file),
  "metadata_file does not exist" = file.exists(opts$metadata_file),
  "celltype_group_file does not exist" = file.exists(opts$celltype_group_file),
  "celltype_dir does not exist" = dir.exists(opts$celltype_dir)
)

# define, check consensus cell type files
celltype_files <- list.files(
  path = opts$celltype_dir,
  pattern = "_processed_consensus-cell-types\\.tsv\\.gz$",
  recursive = TRUE,
  full.names = TRUE
) |>
  # set names as library id
  purrr::set_names(
    \(x) {
      stringr::str_split_i(basename(x), "_", 1)
    }
  )
stopifnot("Could not find celltype files" = length(celltype_files) > 0)

# read input data -----------
celltype_group_df <- readr::read_tsv(opts$celltype_group_file)

merged_sce <- readRDS(opts$merged_sce_file)

# determine non-PDX libraries to keep
tumor_library_ids <- readr::read_tsv(opts$metadata_file) |>
  dplyr::filter(!is_xenograft) |>
  dplyr::pull(scpca_library_id)

# Prepare data frame of all cell types
celltype_df <- celltype_files |>
  # consider only non-PDX
  purrr::keep_at(tumor_library_ids) |>
  purrr::map(readr::read_tsv) |>
  purrr::list_rbind() |>
  dplyr::mutate(sce_cell_id = glue::glue("{library_id}-{barcodes}")) |>
  dplyr::select(sce_cell_id, consensus_annotation) |>
  dplyr::inner_join(celltype_group_df)


# Define vectors of cell ids for each reference -------

# All-normal reference
# as determined in ../../exploratory-notebooks/SCPCP000004/01_nb-consensus-cell-types.Rmd,
# use the following groups in `all-normal`: endothelial, epithelial, adipocytes, immune
normal_groups <- c("endothelial", "epithelial", "adipocytes", "immune")
all_normal_ids <- celltype_df |>
  dplyr::filter(reference_cell_group %in% normal_groups) |>
  dplyr::pull(sce_cell_id)

# Immune reference
immune_ids <- celltype_df |>
  dplyr::filter(reference_cell_group == "immune") |>
  dplyr::pull(sce_cell_id)

# Endothelial reference
endo_ids <- celltype_df |>
  dplyr::filter(reference_cell_group == "endothelial") |>
  dplyr::pull(sce_cell_id)

# Subset the SCE to create references -------------

# first create the all-normal, and then the immune and endo as subsets from that one
normal_reference <- merged_sce[, all_normal_ids] |>
  clean_sce()

# clean up for space
rm(merged_sce)
gc()

# Ensure that the consensus cell type is recorded in reference SCE if not present
if (!("consensus_annotation" %in% colnames(colData(normal_reference)))) {
  colData(normal_reference) <- consensus_to_coldata(colData(normal_reference), celltype_df)
}

immune_reference <- normal_reference[, immune_ids]
endo_reference <- normal_reference[, endo_ids]

# Export references ---------
readr::write_rds(normal_reference, opts$reference_normal_rds, compress = "gz")
readr::write_rds(immune_reference, opts$reference_immune_rds, compress = "gz")
readr::write_rds(endo_reference, opts$reference_endo_rds, compress = "gz")
