#!/usr/bin/env Rscript
#
# This script exports an SCE and AnnData version of a given NBAtlas Seurat object
# The SCE object retains the raw counts, normalized counts, and cell metadata
# The AnnData object retains only the raw counts for the top 2000 high-variance genes, and cell metadata
#  This allows for a smaller object export and lower memory usage during SCVI/SCANVI training
# In addition, a text file with the top 2000 high-variance genes is exported
#
# During processing, one piece of metadata in the object is further updated:
# The `Platform` value for the Costa2022 Study should be `10x_v3.1` and not `10x_v3`
# See this issue discussion for context:
# https://github.com/AlexsLemonade/OpenScPCA-analysis/pull/1231#discussion_r2226070913

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
    opt_str = c("--nbatlas_hvg_file"),
    type = "character",
    help = "Path to output text file to save top 2000 HVGs of the NBAtlas object. This is only exported if the anndata_file is provided."
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
# - update Costa2022 Platform to 10X_v3.1

reducedDims(nbatlas_sce) <- NULL

colData(nbatlas_sce) <- colData(nbatlas_sce) |>
  as.data.frame() |>
  dplyr::mutate(
    cell_id = rownames(colData(nbatlas_sce)),
    in_tumor_zoom = cell_id %in% tumor_cells,
    Platform = ifelse(
      Study == "Costa2022",
      "10X_v3.1",
      Platform)
  ) |>
  DataFrame(row.names = rownames(colData(nbatlas_sce)))


# export reformatted NBAtlas objects if requested
if (!is.null(opts$sce_file)) {
  readr::write_rds(
    nbatlas_sce,
    opts$sce_file,
    compress = "gz"
  )
}

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

  # export the AnnData object
  zellkonverter::writeH5AD(
    nbatlas_sce,
    opts$anndata_file,
    X_name = "counts",
    compression = "gzip"
  )

  # export text file with the HVGs
  readr::write_lines(hv_genes, opts$nbatlas_hvg_file)
}
