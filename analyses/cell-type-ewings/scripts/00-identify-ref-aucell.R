#!/usr/bin/env Rscript

# This script is used to identify a reference threshold to use for classifying cells with AUCell
# The reference SCE is run through AUCell using the marker gene list 
# The auc threshold output by AUCell for this sample is returned as the output

project_root <- here::here()
renv::load(project_root)

library(optparse)
library(SingleCellExperiment)

option_list <- list(
  make_option(
    opt_str = c("--sce_file"),
    type = "character",
    help = "Path to RDS file containing a processed SingleCellExperiment object from scpca-nf 
      to be used to run AUCell and determine an auc threshold."
  ),
  make_option(
    opt_str = c("--marker_genes_file"),
    type = "character",
    default = file.path(project_root, "references", "tumor-marker-genes.tsv"),
    help = "Path to file containing all marker genes used to run AUCell. 
      Must contain `cell_type` and `ensembl_gene_id` columns.
      Any marker genes where `cell_type` is 'Tumor' will be used."
  ),
  make_option(
    opt_str = c("--output_file"),
    type = "character",
    help = "Path to file where AUC threshold will be stored."
  )
)

# Parse options
opt <- parse_args(OptionParser(option_list = option_list))

# Set up -----------------------------------------------------------------------

# make sure input files exist 
stopifnot(
  "sce file does not exist" = file.exists(opt$sce_file),
  "marker genes file does not exist" = file.exists(opt$marker_genes_file)
)

# read in SCE
sce <- readr::read_rds(opt$sce_file)

# get list of marker genes to use with AUCell
marker_genes <- readr::read_tsv(opt$marker_genes_file, show_col_types = FALSE) |> 
  # account for genes being from multiple sources
  dplyr::select(cell_type, ensembl_gene_id) |> 
  dplyr::distinct() |> 
  dplyr::filter(cell_type == "tumor") |> 
  dplyr::pull(ensembl_gene_id)

# turn it into a gene set 
marker_gene_set <- marker_genes |> 
  GSEABase::GeneSet(setName = "markers")

# AUCell -----------------------------------------------------------------------

# run AUCell
auc_results <- AUCell::AUCell_run(counts(sce), marker_gene_set)

# assign cells and identify auc thresholds
# don't output plots 
auc_assignments <- AUCell::AUCell_exploreThresholds(
  auc_results, 
  assign = TRUE, 
  plotHist = FALSE
  )

# pull out threshold used for assigning cells from ref sample 
ref_threshold <- auc_assignments$markers$aucThr$selected

# save threshold to output file 
writeLines(ref_threshold, opt$output_file)

