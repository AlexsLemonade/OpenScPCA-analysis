#!/usr/bin/env Rscript

# Run `copyKAT` for one sample with or without a healthy reference
#
# USAGE:
# Rscript copyKAT.R \
#   --sample_id SCPCS000194 \
#   --use_reference <ref, noref> \
#   --distance <spearman, euclidean> \
#   --ncore 8
#
# Additional optional arguments include:
# --score_threshold: Annotation prediction score threshold to use when identifying cells to use in reference
# --seed: Integer to set the random seed

library(optparse)
library(Seurat)
library(copykat)

# Parse arguments --------------------------------------------------------------
# set up arguments
option_list <- list(
  make_option(
    opt_str = c("-s", "--sample_id"),
    type = "character",
    default = "SCPCS000179",
    help = "The sample_id of the sample to be used for inference of genomic copy number using copyKAT "
  ),
  make_option(
    opt_str = c("-c", "--n_core"),
    type = "integer",
    default = 16,
    help = "number of cores used to run copyKAT"
  ),
  make_option(
    opt_str = c("-d", "--distance"),
    type = "character",
    default = "euclidean",
    help = "method used to calculate distance in copyKAT"
  ),
  make_option(
    opt_str = c("-r", "--use_reference"),
    type = "character",
    default = "ref",
    help = "either to run copyKAT with or without reference normal cells"
  ),
  make_option(
    opt_str = c("-t", "--score_threshold"),
    type = "numeric",
    default = 0.85,
    help = "Threshold prediction score from label transfer to consider a normal cell in the reference"
  ),
  make_option(
    opt_str = "--seed",
    type = "integer",
    default = 12345,
    help = "random seed to set"
  )
)

opts <- parse_args(OptionParser(option_list = option_list))

# Note that the version of copyKAT used here overrides random seeds, so while we set one, it isn't used:
# https://github.com/navinlabcode/copykat/blob/d7d6569ae9e30bf774908301af312f626de4cbd5/R/copykat.R#L33
set.seed(opts$seed)

# paths to data ----------------------------------------------------------------

# The base path for the OpenScPCA repository, found by its (hidden) .git directory
repository_base <- rprojroot::find_root(rprojroot::is_git_root)
# The path to this module
module_base <- file.path(repository_base, "analyses", "cell-type-wilms-tumor-06")
# Path to the result directory
result_dir <- file.path(module_base, "results", opts$sample_id)


# Create directories to save the results of copykat with/without reference using opts$distance
output_dir <- file.path(result_dir, "05_copyKAT", opts$use_reference, opts$distance)
fs::dir_create(output_dir)

# Define output file names
output_rds <- file.path(
  output_dir,
  glue::glue("05_copykat_{opts$sample_id}_{opts$use_reference}_distance-{opts$distance}.rds")
)

output_heatmap_file <- file.path(
  output_dir,
  glue::glue("05_copykat_{opts$sample_id}_{opts$use_reference}_distance-{opts$distance}_copykat_heatmap.jpg")
)

output_prediction_file <- file.path(
  output_dir,
  glue::glue("05_copykat_{opts$sample_id}_{opts$use_reference}_distance-{opts$distance}_copykat_prediction.txt")
)

output_cna_file <- file.path(
  output_dir,
  glue::glue("05_copykat_{opts$sample_id}_{opts$use_reference}_distance-{opts$distance}_copykat_CNA_results.txt")
)

# Read in data -----------------------------------------------------------------
srat <- readRDS(
  file.path(result_dir, paste0("02b-fetal_kidney_label-transfer_", opts$sample_id, ".Rds"))
)

# Extract raw counts -----------------------------------------------------------
exp.rawdata <- GetAssayData(object = srat, assay = "RNA", layer = "counts")

# Define normal cells, if reference should be used  ----------------------------
if (opts$use_reference == "ref") {
  normal_cells <- WhichCells(
    object = srat,
    expression = fetal_kidney_predicted.compartment %in% c("endothelium", "immune") &
      fetal_kidney_predicted.compartment.score > opts$score_threshold
  )
} else {
  normal_cells <- ""
}


# Run copyKAT without reference ------------------------------------------------

copykat.ref <- copykat(
  rawmat = exp.rawdata,
  sam.name = opts$sample_id,
  distance = opts$distance,
  norm.cell.names = normal_cells,
  genome = "hg20",
  n.cores = opts$n_core,
  id.type = "E",
  plot.genes = FALSE,
  output.seg = FALSE,
  KS.cut = 0.05
)

# Clean up and copykat output reference ----------------------------------------

saveRDS(copykat.ref, output_rds)

# move output files to final location
fs::file_move(glue::glue("{opts$sample_id}_copykat_prediction.txt"), output_prediction_file)
fs::file_move(glue::glue("{opts$sample_id}_copykat_heatmap.jpeg"), output_heatmap_file)
fs::file_move(glue::glue("{opts$sample_id}_copykat_CNA_results.txt"), output_cna_file)

# remove extra files we don't need to save
fs::file_delete(glue::glue("{opts$sample_id}_copykat_CNA_raw_results_gene_by_cell.txt"))
fs::file_delete(glue::glue("{opts$sample_id}_copykat_clustering_results.rds"))
