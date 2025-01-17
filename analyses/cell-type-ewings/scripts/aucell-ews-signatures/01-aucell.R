#!/usr/bin/env Rscript

# This script is used to run `AUCell` on a single SCE object for a set of marker gene sets 
# gene sets used are custom gene sets and a set of Ewing specific gene sets from MsigDB
# the results are exported as a single TSV file with the following columns: 
# `gene_set`, `barcodes`, `auc`, and `auc_threshold`

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
    opt_str = c("--custom_geneset_dir"),
    type = "character",
    default = file.path(project_root, "references", "gene_signatures"), 
    help = "Path to directory containing custom gene sets to use with AUCell. 
      All TSV files in the directory will be used and must contain the `ensembl_gene_id` column.
      File names will be used as the name of the gene set."
  ),
  make_option(
    opt_str = c("--max_rank_threshold"),
    type = "double",
    default = 0.01,
    help = "Percentage of genes detected to set as the `aucMaxRank`. 
      Must be a number between 0 and 1."
  ),
  make_option(
    opt_str = c("--output_file"),
    type = "character",
    help = "Path to file where results will be saved"
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
    default = 2025,
    help = "A random seed for reproducibility."
  )
)

# Parse options
opt <- parse_args(OptionParser(option_list = option_list))

# Set up -----------------------------------------------------------------------

# make sure input files exist
stopifnot(
  "sce file does not exist" = file.exists(opt$sce_file),
  "max_rank_threshold must be between 0 and 1" = (opt$max_rank_threshold <= 1 & opt$max_rank_threshold > 0)
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

# check that max rank is number between 0-1

# make sure directory exists for writing output
output_dir <- dir(opt$output_file)
fs::dir_create(output_dir)

# read in SCE
sce <- readr::read_rds(opt$sce_file)

# remove genes that are not detected from SCE object 
genes_to_remove <- rowData(sce)$detected > 0 
filtered_sce <- sce[genes_to_remove , ]

# Prep gene sets ---------------------------------------------------------------

# get list of files with custom gene sets to use 
gene_set_files <- list.files(opt$custom_geneset_dir, recursive = TRUE, full.names = TRUE)
# get names of gene sets using name of the files 
gene_set_names <- stringr::str_replace(basename(gene_set_files), ".tsv", "")

# read in custom gene sets
custom_genes_list <- gene_set_files |> 
  purrr::set_names(gene_set_names) |> 
  purrr::map(\(file) {
    gene_ids <- readr::read_tsv(file) |> 
      dplyr::pull(ensembl_gene_id) |>
      unique()
  })

# grab msigdb gene sets 
ews_gene_sets <- c(
  "staege" = "STAEGE_EWING_FAMILY_TUMOR",
  "miyagawa_up" = "MIYAGAWA_TARGETS_OF_EWSR1_ETS_FUSIONS_UP",
  "miyagawa_down" = "MIYAGAWA_TARGETS_OF_EWSR1_ETS_FUSIONS_DN",
  "zhang"= "ZHANG_TARGETS_OF_EWSR1_FLI1_FUSION", 
  "riggi_up" = "RIGGI_EWING_SARCOMA_PROGENITOR_UP",
  "riggi_down" = "RIGGI_EWING_SARCOMA_PROGENITOR_DN",
  "kinsey_up" = "KINSEY_TARGETS_OF_EWSR1_FLII_FUSION_UP",
  "kinsey_down"= "KINSEY_TARGETS_OF_EWSR1_FLII_FUSION_DN"
)

# pull gene sets from msigbdr
# all gene sets are part of C2, CGP
msig_genes_df <- msigdbr::msigdbr(
  species = "Homo sapiens",
  category = "C2",
  subcategory = "CGP"
) |>
  # only keep relevant gene sets
  dplyr::filter(gs_name %in% ews_gene_sets)

# create named list of genes in each gene set
msig_genes_list <- ews_gene_sets |>
  purrr::map(\(name){
    genes <- msig_genes_df |>
      dplyr::filter(gs_name == name) |>
      dplyr::pull(ensembl_gene) |> 
      unique()
  })

# combine custom and msig
all_genes_list <- c(custom_genes_list, msig_genes_list)

# build GeneSetCollection for AUCell
collection <- all_genes_list |>
  purrr::imap(\(genes, name) GSEABase::GeneSet(genes, setName = name)) |> 
  GSEABase::GeneSetCollection()

# Run AUCell -------------------------------------------------------------------

# run AUCell
counts_mtx <- counts(sce)
max_rank <- ceiling(opt$max_rank_threshold*nrow(counts_mtx))

auc_results <- AUCell::AUCell_run(
  counts_mtx, 
  collection, 
  aucMaxRank = max_rank,
  #BPPARAM = bp_param
)

# Get threshold ----------------------------------------------------------------

# get auc threshold for each geneset 
auc_thresholds <- AUCell::AUCell_exploreThresholds(
  auc_results,
  assign = TRUE,
  plotHist = FALSE
) |> 
  # extract select auc threshold 
  purrr::map_dbl(\(results){
    results$aucThr$selected
  })

# put into a data frame for easy joining with all auc values 
threshold_df <- data.frame(
  gene_set = names(auc_thresholds),
  auc_threshold = auc_thresholds
)

# Combine and export results ---------------------------------------------------

# create data frame with auc for each cell and each geneset 
auc_df <- auc_results@assays@data$AUC |>
  as.data.frame() |>
  tibble::rownames_to_column("gene_set") |> 
  tidyr::pivot_longer(!"gene_set",
                      names_to = "barcodes",
                      values_to = "auc"
  ) |> 
  # add in threshold column 
  dplyr::left_join(threshold_df, by = "gene_set")

# export results as table
readr::write_tsv(auc_df, opt$output_file)
