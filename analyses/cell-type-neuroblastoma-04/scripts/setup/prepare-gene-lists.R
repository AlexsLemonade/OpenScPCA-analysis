#!/usr/bin/env Rscript
#
# This script was used to prepare several TSVs with marker genes to use to validation cell type annotations:
# - A TSV of marker genes from NBAtlas (Bonine et al (2024))
#   - Only gene symbols for which we have known corresponding Ensembl IDs in ScPCA objects are retained
# - A TSV of marker genes from the `cell-type-consensus` module validation genes
#
# To run this script, manually download the Table S2 and Table S5 excel spreadsheets from <https://doi.org/10.1016/j.celrep.2024.114804>.
# Run the script from this directory as follows:
# Rscript prepare-gene-lists.R --table_s2_excel_file <path to Table S2 excel file> --table_s5_excel_file <path to Table S5 excel file>

library(optparse)
module_dir <- here::here()
renv::load(module_dir)

option_list <- list(
  make_option(
    opt_str = c("--table_s2_excel_file"),
    type = "character",
    default = "scratch/table_S2.xlsx",
    help = "Path to Table S2 excel spreadsheet downloaded from <https://doi.org/10.1016/j.celrep.2024.114804>"
  ),
  make_option(
    opt_str = c("--table_s5_excel_file"),
    type = "character",
    default = "scratch/table_S5.xlsx",
    help = "Path to Table S5 excel spreadsheet downloaded from <https://doi.org/10.1016/j.celrep.2024.114804>"
  ),
  make_option(
    opt_str = c("--consensus_marker_genes_url"),
    type = "character",
    default = "https://raw.githubusercontent.com/AlexsLemonade/OpenScPCA-analysis/refs/heads/main/analyses/cell-type-consensus/references/validation-markers.tsv",
    help = "URL to cell-type-consensus validation markers"
  ),
  make_option(
    opt_str = c("--nbatlas_marker_gene_file"),
    type = "character",
    default = file.path(module_dir, "references", "nbatlas-marker-genes.tsv"),
    help = "Path to output TSV file to store NBAtlas marker genes"
  ),
  make_option(
    opt_str = c("--consensus_marker_gene_file"),
    type = "character",
    default = file.path(module_dir, "references", "consensus-marker-genes.tsv"),
    help = "Path to output TSV file to store marker genes taken from cell-type-consensus validation markers"
  )
)

# Parse options and check arguments
opts <- parse_args(OptionParser(option_list = option_list))
stopifnot(
  "table_s2_excel_file does not exist" = file.exists(opts$table_s2_excel_file),
  "table_s5_excel_file does not exist" = file.exists(opts$table_s5_excel_file),
  "nbatlas_marker_gene_file was not provided" = !is.null(opts$nbatlas_marker_gene_file),
  "consensus_marker_gene_file was not provided" = !is.null(opts$consensus_marker_gene_file)
)

#### NBAtlas marker genes -----------------------------------


s2_keep_sheets <- c("Endothelial", "Fibroblast", "NE", "RBCs", "Schwann", "Stromal other") # only keep non-immune cells
s5_discard_sheets <- c("Doublets") # remove the Doublets genes


# Find all sheet names and limit to those we want; each sheet name is a cell type
s2_sheet_names <- readxl::excel_sheets(opts$table_s2_excel_file) |>
  purrr::set_names()
s2_sheet_names <- s2_sheet_names[s2_keep_sheets]

s5_sheet_names <- readxl::excel_sheets(opts$table_s5_excel_file) |>
  purrr::set_names()
s5_sheet_names <- s5_sheet_names[!(s5_sheet_names %in% s5_discard_sheets)]


scpca_gene_table <- rOpenScPCA::scpca_gene_reference |>
  dplyr::select(gene_ids, gene_symbol_scpca)


# Map over relevant sheets in both spreadsheet to get genes


# First, define a helper function to prepare an excel sheet for joining
# sheet_name: The name of the sheet to read in
# excel_file_path: Path to excel file with the given sheet
prepare_excel_sheet <- function(sheet_name, excel_file_path) {
  readxl::read_excel(excel_file_path, sheet = sheet_name) |>
    # keep only genes we have
    dplyr::inner_join(scpca_gene_table, by = c("gene" = "gene_symbol_scpca")) |>
    # keep pval and LFC so we can rank/subset later as needed
    dplyr::select(ensembl_gene_id = gene_ids, gene_symbol = gene, p_val_adj, avg_log2FC)
}


s2_df <- s2_sheet_names |>
  purrr::map(
    prepare_excel_sheet,
    opts$table_s2_excel_file
  ) |>
  purrr::list_rbind(names_to = "NBAtlas_label")

s5_df <- s5_sheet_names |>
  purrr::map(
    prepare_excel_sheet,
    opts$table_s5_excel_file
  ) |>
  purrr::list_rbind(names_to = "NBAtlas_label")

# Combine data frames and finalize processing


nbatlas_markers_df <- s2_df |>
  dplyr::bind_rows(s5_df) |>
  dplyr::mutate(
    # indicate if up or down-regulated, since there are positive and negative LFCs
    direction = ifelse(avg_log2FC > 0, "up", "down"),
    source = "NBAtlas",
    # recode some labels to match what is in the objects
    NBAtlas_label = dplyr::case_when(
      NBAtlas_label == "NE" ~ "Neuroendocrine",
      NBAtlas_label == "Cycling" ~ "Immune cycling",
      NBAtlas_label == "TOX-KIT NK cell" ~ "TOX2+/KIT+ NK cell",
      .default = NBAtlas_label)
  )

# export
readr::write_tsv(nbatlas_markers_df, opts$nbatlas_marker_gene_file)

#### cell-type-consensus marker genes --------------------------

validation_df <- readr::read_tsv(opts$consensus_marker_genes_url)

consensus_markers_df <- validation_df |>
  # keep all annotations and decide which to use for validation once results are in
  dplyr::select(
    validation_group_annotation,
    ensembl_gene_id,
    gene_symbol,
    gene_observed_count
  ) |>
  dplyr::mutate(
    # all of our genes are up-regulated in the given cell type
    direction = "up",
    source = "cell-type-consensus"
  )

# export
readr::write_tsv(consensus_markers_df, opts$consensus_marker_gene_file)
