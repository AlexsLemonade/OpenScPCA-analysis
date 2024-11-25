#!/usr/bin/env Rscript

# This script performs clustering on a single library across a set of clustering parameters and algorithms
# Clustering results are saved as a TSV file 

suppressPackageStartupMessages({
  library(optparse)
  library(SingleCellExperiment)
  library(rOpenScPCA)
})

option_list <- list(
  make_option(
    opt_str = c("--sce_file"),
    type = "character",
    help = "Path to processed SCE file to use for calculating clusters."
  ),
  make_option(
    opt_str = c("--output_file"),
    type = "character",
    help = "Path to TSV file to save clustering results."
  ),
  make_option(
    opt_str = c("--louvain"),
    action = "store_true",
    default = FALSE,
    help = "Indicate whether or not to perform Louvain clustering."
  ), 
  make_option(
    opt_str = c("--leiden_mod"),
    action = "store_true",
    default = FALSE,
    help = "Indicate whether or not to perform Leiden clustering with the modularity objective function."
  ), 
  make_option(
    opt_str = c("--leiden_cpm"),
    action = "store_true",
    default = FALSE,
    help = "Indicate whether or not to perform Leiden clustering with the CPM objective function."
  ), 
  make_option(
    opt_str = c("--nn_range"),
    default = "5,10,15,20,25,30,35,40",
    type = "character", 
    help = "Comma separated list with a range of values to use for nearest neighbors parameter."
  ), 
  make_option(
    opt_str = c("--louvain_res_range"),
    default = ".5,1,1.5",
    type = "character",
    help = "Comma separated list with a range of values to use for resolution with the Louvain clustering algorithm.
      Only used with the `--louvain` flag."
  ),
  make_option(
    opt_str = c("--mod_res_range"),
    default = ".5,1,1.5",
    type = "character",
    help = "Comma separated list with a range of values to use for resolution with the Leiden clustering algorithm and modularity objective function.
      Only used with the `--leiden_mod` flag."
  ),
  make_option(
    opt_str = c("--cpm_res_range"),
    default = "0.001,0.005,0.01",
    type = "character",
    help = "Range of values to use for resolution with the Leiden clustering algorithm and CPM objective function.
      Only used with the `--leiden_cpm` flag."
  ),
  make_option(
    opt_str = c("-t", "--threads"),
    type = "integer",
    default = 4,
    help = "Number of multiprocessing threads to use."
  ),
  make_option(
    opt_str = c("--seed"),
    type = "integer",
    default = 2024,
    help = "A random seed for reproducibility."
  )
)

# Parse options
opt <- parse_args(OptionParser(option_list = option_list))

# Set up -----------------------------------------------------------------------

# set seed
set.seed(opt$seed)

# make sure inputs exist
stopifnot(
  "SCE file does not exist" = file.exists(opt$sce_file),
  "Must indicate at least one clustering algorithm to use, one of `--louvain`, `--leiden_mod`, or `--leiden_cpm`" = opt$louvain && opt$leiden_cpm && opt$leiden_cpm
)

# create output directory if it doesn't exist
output_dir <- dirname(opt$output_file)
fs::dir_create(output_dir)


# read in input sce files
sce <- readr::read_rds(opt$sce_file)

# get nn and resolution as numbers
# define a helper function to read in comma separated lists
split_numbers <- function(param) {
  param |>
    stringr::str_split_1(param, ",") |>
    as.numeric()
 }

nn_range <- split_numbers(opt$nn_range)
louvain_res_range <- split_numbers(opt$louvain_res_range)
mod_res_range <- split_numbers(opt$mod_res_range)
cpm_res_range <- split_numbers(opt$cpm_res_range)

# Clustering -------------------------------------------------------------------

# make a list of clustering options
cluster_opt_list <- list()
if(opt$louvain){
  cluster_opt_list$louvain <- list(
    algorithm = "louvain",
    nn_range = nn_range, 
    res_range = louvain_res_range
  )
}

if(opt$leiden_mod){
  cluster_opt_list$leiden_mod <- list(
    algorithm = "leiden",
    objective_function = "modularity", 
    nn_range = nn_range, 
    res_range = mod_res_range
  )
} 

if(opt$leiden_cpm){
  cluster_opt_list$leiden_cpm <- list(
    algorithm = "leiden",
    objective_function = "CPM", 
    nn_range = nn_range, 
    res_range = cpm_res_range
  )
}


cluster_results_df <- cluster_opt_list |> 
  purrr::map(\(opt_list){
    rOpenScPCA::sweep_clusters(
      sce,
      algorithm = opt_list$algorithm,
      objective_function = opt_list$objective_function,
      nn = opt_list$nn_range, 
      resolution = opt_list$res_range,
      threads = opt$threads,
      seed = opt$seed
    ) |>
      dplyr::bind_rows()
  }) |> 
  dplyr::bind_rows(.id = "cluster_method")
  

# export results 
readr::write_tsv(cluster_results_df, opt$output_file)
