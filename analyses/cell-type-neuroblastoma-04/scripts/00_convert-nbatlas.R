#!/usr/bin/env Rscript
# This script exports an SCE and AnnData version of a given NBAtlas Seurat object
# We retain only the raw counts and cell metadata in the converted objects, since normalized will not be used

library(optparse)

option_list <- list(
  make_option(
    opt_str = c("--seurat_file"),
    type = "character",
    default = "",
    help = "Path to Seurat version of an NBAtlas object"
  ),
  make_option(
    opt_str = c("--tumor_cells_file"),
    type = "character",
    default = "",
    help = "Path to text file listing all tumor cells in the full NBAtlas."
  ),
  make_option(
    opt_str = c("--sce_file"),
    type = "character",
    default = "",
    help = "Path to output RDS file to hold an SCE version of the NBAtlas object"
  ),
  make_option(
    opt_str = c("--anndata_file"),
    type = "character",
    default = "",
    help = "Path to output H5AD file to hold an AnnData version of the NBAtlas object"
  )
)

# Parse options and check arguments
opts <- parse_args(OptionParser(option_list = option_list))

stopifnot(
  "seurat_file does not exist" = file.exists(opts$seurat_file),
  "tumor_cells_file does not exist" = file.exists(opts$tumor_cells_file)
)

# load the bigger libraries after passing checks
suppressPackageStartupMessages({
  library(Seurat)
  library(SingleCellExperiment)
  library(zellkonverter)
})

# read input files
nbatlas_seurat <- readRDS(opts$seurat_file)
tumor_cells <- readLines(opts$tumor_cells_file)

# convert Seurat to SCE object directly, to save space in the final object
nbatlas_sce <- SingleCellExperiment(
  assays = list(counts = nbatlas_seurat[["RNA"]]$counts)
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
