#!/usr/bin/env Rscript

# This script converts the datasets used for benchmarking into SCE and AnnData objects,
# mirroring ScPCA dataset format.
# The original data files were obtained from Zenodo and are named <dataset name>.rds.
# Each is an RDS object containing a list of two items:
#  [[1]] A raw counts matrix
#  [[2]] A vector of doublet/singlet calls for each barcode
# For each dataset, we normalize counts and calculate PCA, but we perform no additional filtering.
# The exported files are named `<dataset name>-<sce/anndata>.<rds/h5ad>.
# Each has a variable `ground_truth_doublets` representing the singlet/doublet calls:
#   - In the SCE file, this is in the colData slot
#   - In the AnnData file, this is in the obs slot

# Load libraries
library(SingleCellExperiment)
library(optparse)


option_list <- list(
  make_option(
    "--dataset_name",
    type = "character",
    help = "Name of dataset to process."
  ),
  make_option(
    c("--input_dir"),
    type = "character",
    help = "Directory containing data files to format."
  ),
  make_option(
    c("--output_dir"),
    type = "character",
    help = "Directory to save formatted files to."
  )
)
opts <- parse_args(OptionParser(option_list = option_list))
fs::dir_create(opts$output_dir)

input_file <- file.path(opts$input_dir, glue::glue("{opts$dataset_name}.rds"))
if (!file.exists(input_file)) {
  stop(
    glue::glue("Input file could not be found at: `{input_file}`.")
  )
}

output_sce_file <- file.path(opts$output_dir, glue::glue("{opts$dataset_name}.rds"))
output_anndata_file <- file.path(opts$output_dir, glue::glue("{opts$dataset_name}.h5ad"))

dat <- readRDS(input_file)
mat <- dat[[1]] # raw counts matrix
calls <- dat[[2]] # "singlet" or "doublet"

# Create, process, and export SCE ------------
sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = mat))
sce$ground_truth_doublets <- calls

# normalization
qclust <- scran::quickCluster(sce)
sce <- scran::computeSumFactors(sce, clusters = qclust) |>
  scuttle::logNormCounts()

# PCA
gene_variance <- scran::modelGeneVar(sce)
var_genes <- scran::getTopHVGs(gene_variance, n = 2000)
sce <- scater::runPCA(
  sce,
  ncomponents = 20,
  subset_row = var_genes
)

# export
readr::write_rds(sce, output_sce_file)

# Export AnnData version -----------
zellkonverter::writeH5AD(sce, output_anndata_file)
