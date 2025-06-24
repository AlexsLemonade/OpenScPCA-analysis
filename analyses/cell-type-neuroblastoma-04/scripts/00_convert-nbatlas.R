#!/usr/bin/env Rscript
# This script exports an SCE and AnnData version of a given NBAtlas Seurat object
# We retain only the raw counts, normalized counts, and cell metadata in the converted objects

library(optparse)

option_list <- list(
  make_option(
    opt_str = c("--seurat_file"),
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
  "seurat_file does not exist" = file.exists(opts$seurat_file),
  "tumor_metadata_file does not exist" = file.exists(opts$tumor_metadata_file),
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

# read input files and determine tumor cell ids
nbatlas_seurat <- readRDS(opts$seurat_file)
tumor_cells <- readRDS(opts$tumor_metadata_file) |>
  rownames()

# convert Seurat to SCE object directly, to save space in the final object
nbatlas_sce <- SingleCellExperiment(
  assays = list(
    counts = nbatlas_seurat[["RNA"]]$counts,
    logcounts = nbatlas_seurat[["RNA"]]$data
  )
)

# add in colData, including updated cell labels with `neuroendocrine-tumor` label for tumor cells
colData(nbatlas_sce) <- nbatlas_seurat@meta.data |>
  dplyr::mutate(
    cell_id = rownames(colData(nbatlas_sce)),
    NBAtlas_label = ifelse(
      cell_id %in% tumor_cells,
      "Neuroendocrine-tumor",
      Cell_type
    )
  ) |>
  DataFrame(row.names = rownames(colData(nbatlas_sce)))

# remove Seurat file to save space
rm(nbatlas_seurat)
gc()

# if testing, subset the SCE to fewer cells: keep 5% of each label
if (opts$testing) {
  keep_cells <- colData(nbatlas_sce) |>
    as.data.frame() |>
    dplyr::group_by(NBAtlas_label) |>
    dplyr::sample_frac(0.05) |>
    dplyr::pull(cell_id)

  nbatlas_sce <- nbatlas_sce[, keep_cells]
}

# export reformatted NBAtlas objects: SCE and AnnData
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
