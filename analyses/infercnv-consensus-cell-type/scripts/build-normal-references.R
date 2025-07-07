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
  file.path(module_dir, "scripts", "build-normal-reference", "utils.R")
)

option_list <- list(
  make_option(
    opt_str = "--merged_sce_file",
    type = "character",
    default = "../../data/current/results/merge-sce/SCPCP000004/SCPCP000004_merged.rds",
    help = "Path to the merged SingleCellExperiment object for the given project"
  ),
  make_option(
    opt_str = "--metadata_file",
    type = "character",
    default = "../../data/current/SCPCP000004/single_cell_metadata.tsv",
    help = "Path to single_cell_metadata.tsv file for the given project"
  ),
  make_option(
    opt_str = "--reference_groups",
    type = "character",
    default = "",
    help = "Comma-separated list of reference group names to build references for.
      Use the string 'all-normal' to include a combined set of normal cell types, which can be further specifed with the `normal_reference_groups` argument."
  ),
  make_option(
    opt_str = "--normal_reference_groups",
    type = "character",
    default = "",
    help = "If 'all-normal' is provided to the `reference_groups` arguments, this argument can be used to specify which cell types to include in that category."
  ),
  make_option(
    opt_str = "--reference_celltype_group_file",
    type = "character",
    default = file.path(module_dir, "references", "reference-cell-groups.tsv"),
    help = "Path to file mapping consensus cell types to reference cell type groups"
  ),
  make_option(
    opt_str = "--exclude_samples",
    type = "character",
    default = "",
    help = "Comma-separated list of any samples which should be excluded from the references"
  ),
  make_option(
    opt_str = "--reference_output_dir",
    type = "character",
    help = "Path to directory where normal reference files, named as `{reference cell group}.rds`,will be saved"
  )
)



# Parse options ----------------
opts <- parse_args(OptionParser(option_list = option_list))

# check input files and paths
stopifnot(
  "merged_sce_file does not exist" = file.exists(opts$merged_sce_file),
  "metadata_file does not exist" = file.exists(opts$metadata_file),
  "reference_celltype_group_file does not exist" = file.exists(opts$reference_celltype_group_file)
)

# read input data -----------
celltype_group_df <- readr::read_tsv(opts$reference_celltype_group_file)

merged_sce <- readRDS(opts$merged_sce_file)

# keep only tumor libraries
tumor_library_ids <- readr::read_tsv(opts$metadata_file) |>
  dplyr::filter(!is_xenograft, !is_cell_line) |>
  dplyr::pull(scpca_library_id)

celltype_df <- merged_sce |>
  colData() |>
  as.data.frame() |>
  # consider only non-PDX
  dplyr::filter(library_id %in% tumor_library_ids) |>
  head()
dplyr::select(barcodes, library_id, consensus_annotation = consensus_celltype_annotation) |>
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
