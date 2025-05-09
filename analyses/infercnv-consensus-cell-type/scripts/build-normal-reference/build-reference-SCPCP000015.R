#!/usr/bin/env Rscript
# This script creates normal references for use with Ewing sarcoma (SCPCP000015) samples
# We create three references:
# - `endo`: All endothelial cells
# - `immune`: All immune cells
# - `endo-immune`: All endothelial and immune cells
# References exclude any cells which the `cell-type-ewings` OpenScPCA-analysis module labeled as "tumor"

suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(optparse)
})

option_list <- list(
  make_option(
    opt_str = "--merged_sce_file",
    type = "character",
    help = "Path to the merged SingleCellExperiment object"
  ),
  make_option(
    opt_str = "--cell_type_ewings_dir",
    type = "character",
    help = "Path to directory containing results from the `cell-type-ewings` module"
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
  ), make_option(
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
  "cell_type_ewings_dir does not exist" = dir.exists(opts$cell_type_ewings_dir)
)

# Paths -----------------

# find and check the cell-type-ewings files
celltype_files <- list.files(
  path = opts$cell_type_ewings_dir,
  pattern = "_ewing-celltype-assignments\\.tsv$",
  recursive = TRUE,
  full.names = TRUE
)
stopifnot(
  "Could not find celltype files" = length(celltype_files) > 0
)

merged_sce <- readRDS(opts$merged_sce_file)


# Define immune cells -------------
immune_celltypes <- readr::read_tsv(opts$immune_ref_url) |>
  dplyr::pull(consensus_annotation)

# Define data frames of cells to include in each reference ----------

# All consensus annotations
consensus_df <- celltype_files |>
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
endo_cell_types <- c(
  "endothelial cell",
  "blood vessel endothelial cell",
  "microvascular endothelial cell"
)
endo_cell_ids <- consensus_df |>
  dplyr::filter(
    consensus_annotation %in% endo_cell_types,
    !(stringr::str_detect(ewing_annotation, "tumor"))
  ) |>
  dplyr::pull(sce_cell_id)


# Subset the SCEs to create references -------------

# Define helper function to remove unneeded slots
# from a reference SCE to save space
clean_sce <- function(sce) {
  logcounts(sce) <- NULL
  assay(sce, "spliced") <- NULL
  reducedDim(sce, "PCA") <- NULL
  reducedDim(sce, "UMAP") <- NULL

  # ensure the counts matrix is sparse
  counts(sce) <- as(counts(sce), "CsparseMatrix")

  return(sce)
}


# First make the combined endo+immune reference
combined_reference <- merged_sce[, c(immune_cell_ids, endo_cell_ids)] |>
  # remove unneeded slots
  clean_sce()


# Ensure that the consensus cell type is recorded in reference SCE if not present
if (!("consensus_annotation" %in% colnames(colData(combined_reference)))) {
  colData(combined_reference) <- colData(combined_reference) |>
    as.data.frame() |>
    # temporarily make the rownames a column so we can join consensus
    tibble::rownames_to_column(var = "sce_cell_id") |>
    dplyr::left_join(consensus_df, by = "sce_cell_id") |>
    dplyr::select(-sce_cell_id) |>
    # make it a DataFrame again
    DataFrame(row.names = rownames(colData(combined_reference)))
}


# Create the immune reference
immune_reference <- combined_reference[, immune_cell_ids]

# Create the endo reference
endo_reference <- combined_reference[, endo_cell_ids]


# Export references ---------
readr::write_rds(immune_reference, opts$reference_immune, compress = "gz")
readr::write_rds(endo_reference, opts$reference_endo, compress = "gz")
readr::write_rds(combined_reference, opts$reference_endo_immune, compress = "gz")
