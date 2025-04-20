#!/usr/bin/env Rscript
# This script creates normal references for use with Ewing sarcoma (SCPCP000015) samples
# Cells to include in the immune-related references were determined in ../exploratory-notebooks/01_ewings-consensus-cell-types.Rmd

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
    opt_str = "--reference_immune",
    type = "character",
    help = "Path to output RDS file to save an SCE file to use as a normal reference with all Ewing immune cells"
  ),
  make_option(
    opt_str = "--reference_immune_subset",
    type = "character",
    help = "Path to output RDS file to save an SCE file to use as a normal reference with a subset Ewing immune cells, specifically macrophages and T cell types"
  ),
  make_option(
    opt_str = "--reference_endo",
    type = "character",
    help = "Path to output RDS file to save an SCE file to use as a normal reference with all Ewing endothelial cells"
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
immune_cells_df <- consensus_df |>
  dplyr::filter(
    consensus_annotation %in% immune_celltypes,
    !(stringr::str_detect(ewing_annotation, "tumor"))
  )

# Subset further to only macrophage and T-cell cell types
immune_subset_cells_df <- immune_cells_df |>
  dplyr::filter(
    consensus_annotation == "macrophage" |
      stringr::str_detect(consensus_annotation, "T cell")
  )

# Subset to the endothelial cells only
endo_cell_types <- c(
  "endothelial cell",
  "blood vessel endothelial cell",
  "microvascular endothelial cell"
)
endo_cells_df <- consensus_df |>
  dplyr::filter(
    consensus_annotation %in% endo_cell_types,
    !(stringr::str_detect(ewing_annotation, "tumor"))
  )

# Subset the SCEs to create references -------------

# Define helper function to remove unneeded slots
# from a reference SCE to save space
clean_sce <- function(sce) {
  logcounts(sce) <- NULL
  assay(sce, "spliced") <- NULL
  reducedDim(sce, "PCA") <- NULL
  reducedDim(sce, "UMAP") <- NULL

  return(sce)
}

# first the reference with all immune cells
all_immune_reference <- merged_sce[, immune_cells_df$sce_cell_id]

# Ensure that the consensus cell type is recorded in reference SCE if not present
if (!("consensus_annotation" %in% colnames(colData(all_immune_reference)))) {
  colData(all_immune_reference) <- colData(all_immune_reference) |>
    as.data.frame() |>
    # temporarily make the rownames a column so we can join consensus
    tibble::rownames_to_column(var = "sce_cell_id") |>
    dplyr::left_join(consensus_df, by = "sce_cell_id") |>
    dplyr::select(-sce_cell_id) |>
    # make it a DataFrame again
    DataFrame(row.names = rownames(colData(all_immune_reference)))
}

# remove some SCE slots to save space
all_immune_reference <- clean_sce(all_immune_reference)

# subset it to create the second reference
subset_immune_reference <- all_immune_reference[, immune_subset_cells_df$sce_cell_id]

# finally, the endo reference
endo_reference <- merged_sce[, endo_cells_df$sce_cell_id] |>
  clean_sce()

# Export references ---------
readr::write_rds(all_immune_reference, opts$reference_immune, compress = "gz")
readr::write_rds(subset_immune_reference, opts$reference_immune_subset, compress = "gz")
readr::write_rds(endo_reference, opts$reference_endo, compress = "gz")
