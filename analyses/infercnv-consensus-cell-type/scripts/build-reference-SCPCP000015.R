#!/usr/bin/env Rscript
# This script creates normal references for use with Ewing sarcoma (SCPCP000015) samples.
# References are created by pooling relevant cell types across all project samples.
# We create three references:
# - `endo`: All endothelial cells
# - `immune`: All immune cells
# - `endo-immune`: All endothelial and immune cells
# References exclude any cells from PDX samples or for which the `cell-type-ewings` OpenScPCA-analysis module labeled as "tumor"

suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(optparse)
})

# source functions
source(here::here("scripts", "utils.R"))

option_list <- list(
  make_option(
    opt_str = "--merged_sce_file",
    type = "character",
    default = "",
    help = "Path to the merged SingleCellExperiment object"
  ),
  make_option(
    opt_str = "--cell_type_ewings_dir",
    type = "character",
    default = "",
    help = "Path to directory containing results from the `cell-type-ewings` module"
  ),
  make_option(
    opt_str = "--metadata_file",
    type = "character",
    default = "",
    help = "Path to single_cell_metadata.tsv file for project SCPCP000015"
  ),
  make_option(
    opt_str = "--reference_endo",
    type = "character",
    help = "Path to output RDS file to save an SCE file to use as a normal reference with all Ewing endothelial cells"
  ),
  make_option(
    opt_str = "--reference_immune",
    type = "character",
    help = "Path to output RDS file to save an SCE file to use as a normal reference with all Ewing immune cells"
  ),
  make_option(
    opt_str = "--reference_endo_immune",
    type = "character",
    help = "Path to output RDS file to save an SCE file to use as a normal reference with all Ewing endothelial and immune cells"
  ),
  make_option(
    opt_str = "--reference_tsv",
    type = "character",
    help = "Path to output TSV file to save cell types included in each reference",
  ),
  make_option(
    opt_str = "--immune_ref_url",
    type = "character",
    default = "https://raw.githubusercontent.com/AlexsLemonade/OpenScPCA-analysis/refs/heads/main/analyses/cell-type-consensus/references/consensus-immune-cell-types.tsv",
    help = "URL of the OpenScPCA consensus cell type immune cells"
  )
)

# Parse options ----------------
opts <- parse_args(OptionParser(option_list = option_list))

# check input directories
stopifnot(
  "merged_sce_file does not exist" = file.exists(opts$merged_sce_file),
  "metadata_file does not exist" = file.exists(opts$metadata_file),
  "cell_type_ewings_dir does not exist" = dir.exists(opts$cell_type_ewings_dir)
)


# Paths -----------------

# find and check the cell-type-ewings files
celltype_files <- list.files(
  path = opts$cell_type_ewings_dir,
  pattern = "_ewing-celltype-assignments\\.tsv$",
  recursive = TRUE,
  full.names = TRUE
) |>
  # name by library_id
  purrr::set_names(
    \(x) {
      stringr::str_split_i(basename(x), "_", 1)
    }
  )

stopifnot(
  "Could not find celltype files" = length(celltype_files) > 0
)


# get library ids to include, which are only the "Tumor" samples
# we'll use this to exclude PDX samples
include_libraries <- readr::read_tsv(opts$metadata_file) |>
  dplyr::filter(sample_type == "Tumor") |>
  dplyr::pull(scpca_library_id)

merged_sce <- readRDS(opts$merged_sce_file)


# Define immune cells -------------
immune_celltypes <- readr::read_tsv(opts$immune_ref_url) |>
  dplyr::pull(consensus_annotation)

# Define data frames of cells to include in each reference ----------

# All consensus annotations
consensus_df <- celltype_files |>
  # only read in libraries of interest
  purrr::keep_at(include_libraries) |>
  purrr::map(readr::read_tsv) |>
  purrr::list_rbind() |>
  dplyr::mutate(sce_cell_id = glue::glue("{library_id}-{barcodes}")) |>
  dplyr::select(sce_cell_id, ewing_annotation, consensus_annotation)

# Pull out immune cells not labeled as tumor
immune_cell_ids <- consensus_df |>
  dplyr::filter(
    consensus_annotation %in% immune_celltypes,
    !(stringr::str_detect(ewing_annotation, "tumor"))
  ) |>
  dplyr::pull(sce_cell_id)

# Subset to the endothelial cells only
endo_celltypes <- c(
  "endothelial cell",
  "blood vessel endothelial cell",
  "microvascular endothelial cell",
  "pericyte"
)
endo_cell_ids <- consensus_df |>
  dplyr::filter(
    consensus_annotation %in% endo_celltypes,
    !(stringr::str_detect(ewing_annotation, "tumor"))
  ) |>
  dplyr::pull(sce_cell_id)


# Subset the SCEs to create references -------------

# First make the combined endo+immune reference
combined_reference <- merged_sce[, c(immune_cell_ids, endo_cell_ids)] |>
  # remove unneeded slots
  clean_sce()

# Ensure that the consensus cell type is recorded in reference SCE if not present
if (!("consensus_annotation" %in% colnames(colData(combined_reference)))) {
  colData(combined_reference) <- consensus_to_coldata(colData(combined_reference), consensus_df)
}

# Create the immune reference
immune_reference <- combined_reference[, immune_cell_ids]

# Create the endo reference
endo_reference <- combined_reference[, endo_cell_ids]


# Export references ---------
readr::write_rds(immune_reference, opts$reference_immune, compress = "gz")
readr::write_rds(endo_reference, opts$reference_endo, compress = "gz")
readr::write_rds(combined_reference, opts$reference_endo_immune, compress = "gz")

# Export TSV of reference cell types ---------------
dplyr::bind_rows(
  data.frame(
    reference_name = "immune",
    consensus_celltype = immune_celltypes
  ),
  data.frame(
    reference_name = "endo",
    consensus_celltype = endo_celltypes
  ),
  data.frame(
    reference_name = "endo-immune",
    consensus_celltype = c(endo_celltypes, immune_celltypes)
  )
) |>
  readr::write_tsv(opts$reference_tsv)
