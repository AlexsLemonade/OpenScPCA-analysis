#!/usr/bin/env Rscript
#
# This script exports an SCE and AnnData version of a given NBAtlas Seurat object
# The SCE object retains the raw counts, normalized counts, and cell metadata
# The AnnData object retains only the raw counts for the top 2000 high-variance genes, and cell metadata
#  This allows for a smaller object export and lower memory usage during SCVI/SCANVI training

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
    opt_str = c("--sce_file"),
    type = "character",
    help = "Path to output RDS file to hold an SCE version of the NBAtlas object. If not provided, the SCE object will not be saved."
  ),
  make_option(
    opt_str = c("--anndata_file"),
    type = "character",
    help = "Path to output H5AD file to hold an AnnData version of the NBAtlas object. If not provided, the AnnData object will not be saved."
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
  "tumor_metadata_file does not exist" = file.exists(opts$tumor_metadata_file)
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

# if testing, subset to fewer cells: keep 5% of each label
# Note that we use the `Cell_type_wImmuneZoomAnnot` for annotation but the
#  `Cell_type` label here to subset, since `Cell_type_wImmuneZoomAnnot` isn't in
#  the NBAtlas subset which is used in testing
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


# export reformatted NBAtlas objects if requested
# if (!is.null(opts$sce_file)) {
#   readr::write_rds(
#     nbatlas_sce,
#     opts$sce_file,
#     compress = "gz"
#   )
# }

if (!is.null(opts$anndata_file)) {
  # We can pare this object down substantially to save space:
  # - subset to top 2000 HVGs (batch-aware)
  # - remove logcounts

  gene_var <- scran::modelGeneVar(
    nbatlas_sce,
    block = nbatlas_sce$Sample
  )
  hv_genes <- scran::getTopHVGs(gene_var, n = 2000)

  nbatlas_sce <- nbatlas_sce[hv_genes,]
  logcounts(nbatlas_sce) <- NULL

  zellkonverter::writeH5AD(
    nbatlas_sce,
    opts$anndata_file,
    X_name = "counts",
    compression = "gzip"
  )
}
