#!/usr/bin/env Rscript
#
# This script was used to prepare several TSVs with marker genes to use to validation cell type annotations:
# - A TSV of marker genes from NBAtlas (Bonine et al (2024))
#   - Only gene symbols for which we have known corresponding Ensembl IDs in ScPCA objects are retained
# - A TSV of marker genes from the `cell-type-consensus` module validation genes
#
# To run this script, manually download the Table S2 excel spreadsheet from <https://doi.org/10.1016/j.celrep.2024.114804>.
# Run the script from this directory as follows:
# Rscript prepare-gene-lists.R --excel_file <path to excel file>

library(optparse)
module_dir <- here::here()
renv::load(module_dir)

option_list <- list(
  make_option(
    opt_str = c("--excel_file"),
    type = "character",
    default = "",
    help = "Path to Table S2 excel spreadsheet downloaded from <https://doi.org/10.1016/j.celrep.2024.114804>"
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
    default = here::here("references/nbatlas-marker-genes.tsv"),
    help = "Path to output TSV file to store NBAtlas marker genes"
  ),
  make_option(
    opt_str = c("--consensus_marker_gene_file"),
    type = "character",
    default = here::here("references/consensus-marker-genes.tsv"),
    help = "Path to output TSV file to store marker genes taken from cell-type-consensus validation markers"
  )
)

# Parse options and check arguments
opts <- parse_args(OptionParser(option_list = option_list))
stopifnot(
  "excel_file does not exist" = file.exists(opts$excel_file),
  "nbatlas_marker_gene_file was not provided" = !is.null(opts$nbatlas_marker_gene_file),
  "consensus_marker_gene_file was not provided" = !is.null(opts$consensus_marker_gene_file)
)

#### NBAtlas marker genes -----------------------------------

# Find all sheet names; each sheet name is a cell type
sheet_names <- readxl::excel_sheets(opts$excel_file)
names(sheet_names) <- sheet_names

scpca_gene_table <- rOpenScPCA::scpca_gene_reference |>
  dplyr::select(gene_ids, gene_symbol_scpca)

# Map over sheets to get genes, keeping only those we have in ScPCA
nbatlas_markers_df <- sheet_names |>
  purrr::map(
    \(x) {
      readxl::read_excel(opts$excel_file, sheet = x) |>
        # keep only genes we have
        dplyr::inner_join(scpca_gene_table, by = c("gene" = "gene_symbol_scpca")) |>
        # keeping pval and LFC will allow us to rank/subset later as needed
        dplyr::select(ensembl_gene_id = gene_ids, gene_symbol = gene, p_val_adj, avg_log2FC)
    }
  ) |>
  purrr::list_rbind(names_to = "NBAtlas_label") |>
  dplyr::mutate(
    # indicate if up or down-regulated, since there are positive and negative LFCs
    direction = ifelse(avg_log2FC > 0, "up", "down"),
    source = "NBAtlas",
    # recode NE -> Neuroendocrine
    mutate(NBAtlas_label = ifelse(NBAtlas_label == "NE", "Neuroendocrine", NBAtlas_label))
  )

# export
readr::write_tsv(nbatlas_markers_df, opts$nbatlas_marker_gene_file)

#### cell-type-consensus marker genes --------------------------

# map of names to match NBAtlas labels with our consensus labels (excluding NE)
# consensus = nbatlas
celltype_map <- c(
  "B cell" = "B cell",
  "endothelial cell" = "Endothelial",
  "fibroblast" = "Fibroblast",
  "myeloid" = "Myeloid",
  "natural killer cell" = "NK cell",
  "dendritic cell" = "pDC",
  "plasma cell" = "Plasma",
  "stromal cell" = "Stromal other", # may not be a great match, but can't hurt
  "T cell" = "T cell"
)

validation_df <- readr::read_tsv(opts$consensus_marker_genes_url)

consensus_markers_df <- validation_df |>
  dplyr::filter(validation_group_annotation %in% names(celltype_map)) |>
  dplyr::select(NBAtlas_label = validation_group_annotation, ensembl_gene_id, gene_symbol) |>
  dplyr::mutate(
    NBAtlas_label = dplyr::recode(NBAtlas_label, !!!celltype_map),
    # all of our genes are up-regulated in the given cell type
    direction = "up",
    source = "cell-type-consensus"
  )

# export
readr::write_tsv(consensus_markers_df, opts$consensus_marker_gene_file)
