#!/usr/bin/env Rscript
#
# This script prepares the merged SCPCP000004 SCE object for input to scANVI:
# - unneeded slots are removed to save space/memory
# - rownames are converted to gene symbols
# - the object is subset to the NBAtlas HVGs
# - fields to match NBAtlas covariate encoding are added to the colData
# - the object is exported as an AnnData object

library(optparse)

option_list <- list(
  make_option(
    opt_str = c("--merged_sce_file"),
    type = "character",
    default = "",
    help = "Path to the SCPCP000004 merged SCE file"
  ),
  make_option(
    opt_str = c("--nbatlas_hvg_file"),
    type = "character",
    default = "",
    help = "Path to text file with top 2000 HVGs of the NBAtlas object, as gene symbols"
  ),
  make_option(
    opt_str = c("--prepared_anndata_file"),
    type = "character",
    default = "",
    help = "Path to output the updated merged AnnData file to use as scANVI query input"
  )
)

# Parse options and check arguments
opts <- parse_args(OptionParser(option_list = option_list))

stopifnot(
  "merged_sce_file does not exist" = file.exists(opts$merged_sce_file),
  "nbatlas_hvg_file does not exist" = file.exists(opts$nbatlas_hvg_file)
)

# load the bigger libraries after passing checks
suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(zellkonverter)
})


# read merged sce
merged_sce <- readRDS(opts$merged_sce_file)

# read gene symbols
hv_gene_symbols <- readr::read_lines(opts$nbatlas_hvg_file)

# remove assays, reducedDims, metadata we won't need
assay(merged_sce, "logcounts") <- NULL
assay(merged_sce, "spliced") <- NULL
reducedDims(merged_sce) <- NULL
metadata(merged_sce) <- list() # also removing since some of this interferes with anndata conversion


# convert to gene symbols
merged_sce <- rOpenScPCA::sce_to_symbols(merged_sce, reference = "sce")

# subset sce
# note that only 1975 genes are present out of 2000, which is 98-99%
# scanvi warns if there is less than 80% overlap, so this is fine. sources:
# https://github.com/scverse/scvi-tools/blob/70564c397b789943230b900500c557f31905d91b/src/scvi/model/base/_archesmixin.py#L40
# https://github.com/scverse/scvi-tools/blob/70564c397b789943230b900500c557f31905d91b/src/scvi/model/base/_archesmixin.py#L479 intersecting_genes <- intersect(rownames(merged_sce), hv_gene_symbols) # 1975 genes; 99% present which is fine per scanvi
intersecting_genes <- intersect(rownames(merged_sce), hv_gene_symbols) # 1975 genes
merged_sce <- merged_sce[intersecting_genes, ]


# add colData columns that match NBAtlas naming
colData(merged_sce) <- colData(merged_sce) |>
  as.data.frame() |>
  dplyr::mutate(
    Sample = library_id,
    Assay = ifelse(suspension_type == "cell", "single-cell", "single-nucleus"),
    Platform = ifelse(stringr::str_detect(tech_version, "10Xv3"), "10X_v3", "10X_v2")
  ) |>
  # recode NAs to support anndata conversion
  # source: https://github.com/AlexsLemonade/scpcaTools/blob/d0fe377284aaa1b4b0647374060e5c699b4c3a48/R/sce_to_anndata.R#L78
  dplyr::mutate(
    dplyr::across(dplyr::where(\(x) all(is.na(x))), as.logical)
  ) |>
  DataFrame(row.names = colnames(merged_sce))


# export as an AnnData object
zellkonverter::writeH5AD(
  merged_sce,
  opts$prepared_anndata_file,
  X_name = "counts",
  compression = "gzip"
)
