#!/usr/bin/env Rscript
#
# This script exports an SCE and AnnData version of a given NBAtlas Seurat object
# We retain only the raw counts, normalized counts, and cell metadata in the converted objects
# In addition, only cells present in the provided `cell_id_file` are included in the final objects

library(optparse)

option_list <- list(
  make_option(
    opt_str = c("--nbatlas_file"),
    type = "character",
    default = "",
    help = "Path to Seurat version of an NBAtlas object"
  ),
  make_option(
    opt_str = c("--tumor_metadata_file"),
    type = "character",
    default = "",
    help = "Path to RDS file with data frame containing NBAtlas tumor metadata."
  ),
  make_option(
    opt_str = c("--cell_id_file"),
    type = "character",
    default = "",
    help = "Path to text file with cell ids present in the full atlas dated `20241203`. Any cell ids not present in this list will be excluded from the converted reference."
  ),
  make_option(
    opt_str = c("--sce_file"),
    type = "character",
    help = "Path to output RDS file to hold an SCE version of the NBAtlas object"
  ),
  make_option(
    opt_str = c("--anndata_file"),
    type = "character",
    help = "Path to output H5AD file to hold an AnnData version of the NBAtlas object"
  ),
  make_option(
    opt_str = c("--testing"),
    action = "store_true",
    default = FALSE,
    help = "Use this flag when running in CI to output a smaller version of NBAtlas for testing"
  ),
  make_option(
    opt_str = c("--seed"),
    type = "integer",
    default = 2025,
    help = "Random seed used to subset NBAtlas when --testing is specified"
  )
)

# Parse options and check arguments
opts <- parse_args(OptionParser(option_list = option_list))

stopifnot(
  "nbatlas_file does not exist" = file.exists(opts$nbatlas_file),
  "tumor_metadata_file does not exist" = file.exists(opts$tumor_metadata_file),
  "cell_id_file does not exist" = file.exists(opts$cell_id_file),
  "sce_file was not provided" = !is.null(opts$sce_file),
  "anndata_file was not provided" = !is.null(opts$anndata_file)
)

# load the bigger libraries after passing checks
suppressPackageStartupMessages({
  library(Seurat)
  library(SingleCellExperiment)
  library(zellkonverter)
})
set.seed(opts$seed)

# read input files and determine relevant cell ids
nbatlas_seurat <- readRDS(opts$nbatlas_file)
tumor_cells <- readRDS(opts$tumor_metadata_file) |>
  rownames()
all_cell_ids <- readr::read_lines(opts$cell_id_file)

# keep only cells that are present in in `all_ids`
nbatlas_seurat <- subset(nbatlas_seurat, cells = all_cell_ids)

# if testing, subset to fewer cells: keep 5% of each label
if (opts$testing) {
  keep_cells <- nbatlas_seurat@meta.data |>
    tibble::rownames_to_column("cell_id") |>
    dplyr::group_by(Cell_type) |>
    dplyr::sample_frac(0.05) |>
    dplyr::pull(cell_id)

  nbatlas_seurat <- subset(nbatlas_seurat, cells = keep_cells)
}

# convert to SCE, using `as` to avoid CI error
# note that Seurat gives some deprecation warnings, but this works
nbatlas_sce <- as.SingleCellExperiment(nbatlas_seurat)

# remove Seurat file to save space
rm(nbatlas_seurat)
gc()

# Update SCE innards:
# - remove reducedDim for space
# - add `cell_id` and `in_tumor_zoom` columns to colData

reducedDim(nbatlas_sce) <- NULL

colData(nbatlas_sce) <- colData(nbatlas_sce) |>
  as.data.frame() |>
  dplyr::mutate(
    cell_id = rownames(colData(nbatlas_sce)),
    in_tumor_zoom = cell_id %in% tumor_cells
  ) |>
  DataFrame(row.names = rownames(colData(nbatlas_sce)))


# export reformatted NBAtlas objects, first SCE then AnnData
readr::write_rds(
  nbatlas_sce,
  opts$sce_file,
  compress = "gz"
)

zellkonverter::writeH5AD(
  nbatlas_sce,
  opts$anndata_file,
  X_name = "counts",
  compression = "gzip"
)
