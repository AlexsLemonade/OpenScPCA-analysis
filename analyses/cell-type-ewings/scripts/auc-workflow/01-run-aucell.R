#!/usr/bin/env Rscript

# This script is used to classify tumor and normal cells using AUCell with a marker gene list 
# The threshold for defining cells as tumor cells should be determined separately by running 
# `00-identify-ref-aucell.R`

project_root <- here::here()

library(optparse)

option_list <- list(
  make_option(
    opt_str = c("--sce_file"),
    type = "character",
    help = "Path to RDS file containing a processed SingleCellExperiment object from scpca-nf 
      to be used to run AUCell."
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
    opt_str = c("--auc_threshold"),
    type = "double",
    default = NULL,
    help = "auc threshold as calculated by AUCell to use for identifying cells as tumor cells.
      If no threshold is provided, will use the threshold determined by AUCell for the input SCE."
  ), 
  make_option(
    opt_str = c("--return_auc"),
    action = "store_true",
    default = FALSE, 
    help = "If used, the AUC value used to classify tumor cells will be printed to stdout"
  ),
  make_option(
    opt_str = c("--output_file"),
    type = "character",
    help = "Path to file where results will be saved. 
      This will be a TSV file containing three columns, `barcodes`, `auc`, `auc_classification`."
  ),
  make_option(
    opt_str = c("-t", "--threads"),
    type = "integer",
    default = 1,
    help = "Number of multiprocessing threads to use."
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

# load SCE 
suppressPackageStartupMessages({
  library(SingleCellExperiment)
})


# set up multiprocessing params
if (opt$threads > 1) {
  bp_param <- BiocParallel::MulticoreParam(opt$threads)
} else {
  bp_param <- BiocParallel::SerialParam()
}

# make sure directory exists for writing output
output_dir <- dir(opt$output_file)
fs::dir_create(output_dir)

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
auc_results <- AUCell::AUCell_run(counts(sce), marker_gene_set, BPPARAM = bp_param)

# get threshold 
if(is.null(opt$auc_threshold)){
  # assign cells and identify auc thresholds
  # don't output plots 
  auc_assignments <- AUCell::AUCell_exploreThresholds(
    auc_results, 
    assign = TRUE, 
    plotHist = FALSE
  )
  
  # pull out threshold used for assigning cells 
  auc_threshold <- auc_assignments$markers$aucThr$selected
  
} else {
  auc_threshold <- opt$auc_threshold
}

# create data frame with auc for each cell and classification
auc_df <- auc_results@assays@data$AUC |> 
  as.data.frame() |> 
  tidyr::pivot_longer(everything(), 
                      names_to = "barcodes", 
                      values_to = "auc") |> 
  dplyr::mutate(auc_classification = dplyr::if_else(auc >= auc_threshold, "Tumor", "Normal"))

# if return auc is true, print to stdout
if(opt$return_auc){
  write(auc_threshold, stdout())
}

# export results as table
readr::write_tsv(auc_df, opt$output_file)
