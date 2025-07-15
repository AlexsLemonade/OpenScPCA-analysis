#!/usr/bin/env Rscript
#
# This script trains a SingleR model from a given NBAtlas object
# The script additionally requires an input SCE object to determine genes to restrict to

suppressWarnings({
  suppressPackageStartupMessages({
    library(optparse)
    library(SingleCellExperiment)
  })
})


option_list <- list(
  make_option(
    opt_str = c("--nbatlas_sce"),
    type = "character",
    default = "",
    help = "Path to an NBAtlas object in SCE format"
  ),
  make_option(
    opt_str = c("--sce_file"),
    type = "character",
    default = "",
    help = "Path to an ScPCA SCE object for determining genes to restrict to"
  ),
  make_option(
    opt_str = c("--singler_model_file"),
    type = "character",
    default = "",
    help = "Path to RDS file to save trained SingleR model"
  ),
  make_option(
    opt_str = c("--aggregate_reference"),
    action = "store_true",
    default = FALSE,
    help = "Whether to aggregate the reference using `SingleR::aggregateReference()` before training. Default: FALSE."
  ),
  make_option(
    opt_str = c("--separate_tumor"),
    action = "store_true",
    default = FALSE,
    help = "Whether to use a separate `Neuroendocrine-tumor` category for cells present in the NBAtlas tumor zoom. Default: TRUE"
  ),
  make_option(
    opt_str = c("--filter_genes"),
    action = "store_true",
    default = FALSE,
    help = "Whether to remove mitochondrial and ribosomal genes from the reference before training the model. Default: FALSE"
  ),
  make_option(
    opt_str = c("--testing"),
    action = "store_true",
    default = FALSE,
    help = "Whether this script is being run on test data. If TRUE, the `Cell_type` column is used for annotation."
  ),
  make_option(
    opt_str = c("--threads"),
    type = "integer",
    default = 4,
    help = "Number of threads for SingleR to use"
  ),
  make_option(
    opt_str = c("--seed"),
    type = "integer",
    default = 2025,
    help = "Random seed"
  )
)

# Parse options and check arguments
opts <- parse_args(OptionParser(option_list = option_list))
stopifnot(
  "nbatlas_sce does not exist" = file.exists(opts$nbatlas_sce),
  "sce_file does not exist" = file.exists(opts$sce_file)
)
set.seed(opts$seed)

if (opts$threads == 1) {
  bp_param <- BiocParallel::SerialParam()
} else {
  bp_param <- BiocParallel::MulticoreParam(opts$threads)
}

# Read atlas
nbatlas_sce <- readRDS(opts$nbatlas_sce)

# Define restrict vector for model training
sce_rowdata <- readRDS(opts$sce_file) |>
  rowData()
restrict_genes <- intersect(
  sce_rowdata$gene_symbol,
  rownames(nbatlas_sce)
)

# if we are filtering genes, remove those from the restrict_genes list
# this way they will not be included in the model
if (opts$filter_genes) {
  remove_genes <- c(
    grep("^MT-", restrict_genes, value = TRUE),
    grep("^RP[SL]", restrict_genes, value = TRUE)
  )
  restrict_genes <- restrict_genes[!(restrict_genes %in% remove_genes)]
}

# In testing, we use the `Cell_type` column
# Otherwise, we use the `Cell_type_wImmuneZoomAnnot` column which has finer-grained immune annotations
if (opts$testing) {
  celltype_column <- "Cell_type"
} else {
  celltype_column <- "Cell_type_wImmuneZoomAnnot"
}

# Define vector of cell labels, considering "tumor" if opts$separate_tumor is TRUE
if (opts$separate_tumor) {
  celltype_label <- ifelse(
    nbatlas_sce$in_tumor_zoom,
    "Neuroendocrine-tumor",
    colData(nbatlas_sce)[[celltype_column]]
  )
} else {
  celltype_label <- colData(nbatlas_sce)[[celltype_column]]
}

# Create and export an aggregated version of the reference
nbatlas_trained <- SingleR::trainSingleR(
  ref = nbatlas_sce,
  labels = celltype_label, # specify vector we created above
  # note the aggregated references are also fairly sparse,
  # so this is appropriate for either type
  de.method = "wilcox",
  restrict = restrict_genes,
  # aggregate as specified
  aggr.ref = opts$aggregate_reference,
  BPPARAM = bp_param
)
readr::write_rds(nbatlas_trained, opts$singler_model_file)
