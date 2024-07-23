#!/usr/bin/env Rscript

# This script is used to calculate gene set scores for a set of EWS-FLI1 gene sets from MsigDB
# ZHANG_TARGETS_OF_EWSR1_FLI1_FUSION - https://www.gsea-msigdb.org/gsea/msigdb/human/geneset/ZHANG_TARGETS_OF_EWSR1_FLI1_FUSION.html
# RIGGI_EWING_SARCOMA_PROGENITOR_UP - https://www.gsea-msigdb.org/gsea/msigdb/human/geneset/RIGGI_EWING_SARCOMA_PROGENITOR_UP.html?ex=1
# SILIGAN_TARGETS_OF_EWS_FLI1_FUSION_DN https://www.gsea-msigdb.org/gsea/msigdb/cards/SILIGAN_TARGETS_OF_EWS_FLI1_FUSION_DN

# output is a TSV where each cell is a row and each gene set score is a column 
# calculates the sum, mean, z-scaled(sum), and z-scaled(mean) for each geneset 

project_root <- here::here()

suppressPackageStartupMessages({
  library(optparse)
  library(SingleCellExperiment)
})

option_list <- list(
  make_option(
    opt_str = c("--sce_file"),
    type = "character",
    help = "Path to RDS file containing a processed SingleCellExperiment object from scpca-nf"
  ),
  make_option(
    opt_str = c("--output_file"),
    type = "character",
    default = NULL,
    help = "Full path to TSV file to save gene set scores"
  )
)

# Parse options
opt <- parse_args(OptionParser(option_list = option_list))

# Define function for calculating gene set score -------------------------------

# function to calculate sum and mean of all genes in a given geneset 
# genes_df is a dataframe with one row per gene and a column `gs_name`
# All genes that have `gs_name == geneset_name` are included in the calculation
calculate_score <- function(sce, genes_df, geneset_name){
  
  # pull out genes that belong to specified geneset 
  genes <- genes_df |> 
    dplyr::filter(gs_name == geneset_name,
                  ensembl_gene %in% rownames(sce)) |> 
    dplyr::pull(ensembl_gene)
  
  # get expression data for all genes in geneset
  gene_exp <- logcounts(sce[genes, ]) |>
    as.matrix() |>
    t() |>
    as.data.frame()
  
  # create data frame with sum, mean, and scaled values 
  scores_df <- data.frame(
    sum = rowSums(gene_exp),
    mean = rowMeans(gene_exp),
    row.names = rownames(gene_exp)
  ) |> 
    dplyr::mutate(
      # need to use as.numeric to get rid of [,1] which will cause an error when trying to save the tsv 
      scaled_sum = as.numeric(scale(sum)),
      scaled_mean = as.numeric(scale(mean))
    ) 
  
  # use geneset name to set column names 
  # only keep the first name of the dataset
  geneset_id <- stringr::word(geneset_name, sep = "_")
  colnames(scores_df) <- glue::glue("{colnames(scores_df)}-{geneset_id}")
  
  return(scores_df)
  
}

# Set up -----------------------------------------------------------------------

# make sure path to sce file exists and output file is provided 
stopifnot("sce_file does not exist" = file.exists(opt$sce_file), 
          "Must provide a --output_file to save results" = !is.null(opt$output_file))

# read in sce file
sce <- readr::read_rds(opt$sce_file)

# library ID 
library_id <- metadata(sce)$library_id

# contstruct and create output folder if not already present 
results_dir <- dirname(opt$output_file)
fs::dir_create(results_dir)

# Grab gene sets ---------------------------------------------------------------

# ews gene set list 
ews_gene_sets <- c(
  "ZHANG_TARGETS_OF_EWSR1_FLI1_FUSION",
  "RIGGI_EWING_SARCOMA_PROGENITOR_UP",
  "SILIGAN_TARGETS_OF_EWS_FLI1_FUSION_DN"
)

# pull gene sets from msigbdr 
# all gene sets are part of C2, CGP
genes_df <- msigdbr::msigdbr(species= "Homo sapiens", 
                             category = "C2",
                             subcategory = "CGP") |> 
  # only keep relevant gene sets 
  dplyr::filter(gs_name %in% ews_gene_sets)


# Calculate gene set score -----------------------------------------------------

# get dataframe with scores for all gene sets 
# each column is a geneset score and each row is a cell 
all_scores_df <- ews_gene_sets |> 
  purrr::map(\(geneset_name) {calculate_score(sce, genes_df, geneset_name)}) |> 
  dplyr::bind_cols() |> 
  tibble::rownames_to_column("barcodes")

# Export to TSV 
readr::write_tsv(all_scores_df, opt$output_file)


