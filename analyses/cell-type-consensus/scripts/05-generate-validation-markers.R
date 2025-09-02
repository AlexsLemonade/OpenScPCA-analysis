#!/usr/bin/env Rscript

# This script is used to create two tables of marker genes to help validate the consensus cell type assignments 
# markers are obtained from CellMarker2.0 for Human 
# Prior to running this script, the CellMarker file must be downloaded using 00-download-cellmarker-ref.sh

# the outputs are the `consensus-markers.tsv` and `validation-markerts.tsv` file s
# that contains a list of the top 10 marker genes for each cell ontology ID, one for consensus cell types and one for broad validation groups
# top marker genes are based on the frequency that marker gene shows up in the total tissues that cell type was found in 
# for some cell types the marker gene list is longer than 10 and includes any genes with equivalent percentage to the 10th gene in the list 

# Function for formatting top marker gene table --------------------------------
# takes a table of marker genes and list of ontology Ids and creates a dataframe with 
# one row per marker gene per ontology Id 
# also includes how often that gene is found as a marker gene for all included ids

create_top_markers_table <- function(
    ontology_ids, #list of ids to include marker genes for
    ontology_column, # prefix for columns that include annotation and ontology Ids for cell types
    annotation_column, 
    celltype_groups_df = consensus_groups_df, # dataframe with ontology column and annotation column
    marker_df = cell_marker_df # data frame with marker genes from cellmarker 
) {
  
  # select unique cell types from 
  celltypes_df <- celltype_groups_df |> 
    dplyr::select({{ontology_column}}, {{annotation_column}}) |> 
    unique()
  
  # filter markers to only those for provided ontology IDs
  filtered_markers_df <- marker_df |> 
    dplyr::filter(cellontology_id %in% ontology_ids)
  
  # get total number of tissues for each cell type 
  total_tissues_df <- filtered_markers_df |> 
    dplyr::select(cellontology_id, tissue_type) |> 
    unique() |>
    dplyr::count(cellontology_id, name = "celltype_total_tissues")
  
  # create data frame with top 10 marker genes per cell type 
  # ranking is based on the percentage of tissues the gene shows up in for a given cell type  
  top_markers_df <- filtered_markers_df |> 
    # first count how many tissues each gene shows up in for each cell type 
    dplyr::count(cellontology_id, gene_symbol, sort = TRUE, name = "number_of_tissues") |> 
    # add total possible tissue types 
    dplyr::left_join(total_tissues_df) |> 
    dplyr::mutate(
      # get the total percent of tissues that have that marker gene in that cell type 
      percent_tissues = round((number_of_tissues/celltype_total_tissues) * 100, 2)
    ) |> 
    # get the top 10 for each group, if there's a tie all will be saved
    dplyr::group_by(cellontology_id) |> 
    dplyr::top_n(10, percent_tissues) |> 
    dplyr::ungroup() |> 
    # add column for total number of times that gene is observed in the table
    dplyr::add_count(gene_symbol, name = "gene_observed_count") |> 
    # bring in ensembl ids
    dplyr::left_join(mapped_ids, by = c("gene_symbol")) |>
    # match to human readable name for ontology group 
    dplyr::left_join(celltypes_df, by = setNames(ontology_column, "cellontology_id")) |> 
    dplyr::relocate(
      rlang::sym(annotation_column), ensembl_gene_id, .after = cellontology_id
    ) |> 
    # arrange genes by cell type and percentage of tissues 
    dplyr::arrange( !!rlang::sym(annotation_column), desc(percent_tissues)) |> 
    # rename cellontology_id column
    dplyr::rename( !!rlang::sym(ontology_column) := "cellontology_id")
  
  return(top_markers_df)
  
}

# Load libraries ---------------------------------------------------------------

# need for mapping library ids 
library(AnnotationDbi)
library(org.Hs.eg.db)

# File paths -------------------------------------------------------------------
module_base <- rprojroot::find_root(rprojroot::is_renv_project)

# cell marker file 
cellmarker_file <- file.path(module_base, "references", "Cell_marker_Human.xlsx")

# cell consensus ref file 
consensus_validation_file <- file.path("references", "consensus-validation-groups.tsv")

# define output file 
validation_group_markers_file <- file.path("references", "validation-markers.tsv")
consensus_celltype_markers_file <- file.path("references", "consensus-markers.tsv")


# Define order of groups -------------------------------------------------------

# set validation group order 
# individual consensus cell types will also use the order from these broader groups 
# for example, all B cell subtypes will be listed first
validation_group_order <- c(
  # lymphocytes
  "B cell",
  "plasma cell",
  "T cell",
  "innate lymphoid cell",
  # myeloid
  "dendritic cell",
  "macrophage",
  "monocyte",
  "myeloid cell",
  "natural killer cell",
  # hpsc 
  "hematopoietic precursor cell",
  # stem cell/ development
  "embryonic cell (metazoa)",
  "stem cell",
  # stromal cells
  "adipocyte",
  "chondrocyte",
  "cardiocyte",
  "endothelial cell",
  "endocrine cell",
  "epithelial cell",
  "fibroblast",
  "kidney cell",
  "lung secretory cell",
  "melanocyte",
  "myofibroblast cell",
  "muscle cell",
  "pericyte",
  "stromal cell",
  # neural cells
  "neuron",
  "astrocyte",
  "glial cell",
  # other 
  "chemoreceptor cell",
  "ciliated cell",
  "histamine secreting cell",
  "mesangial cell",
  "pigment cell",
  "retinal cell",
  "surfactant secreting cell"
)

# Prep references --------------------------------------------------------------

# Read in consensus and validation groups
consensus_groups_df <- readr::read_tsv(consensus_validation_file) |> 
  # set factors and order of validation groups 
  dplyr::mutate(validation_group_annotation = factor(validation_group_annotation, levels = validation_group_order)) |> 
  dplyr::arrange(validation_group_annotation) |> 
  # set factors and order of consensus cell types 
  dplyr::mutate(consensus_annotation = factor(consensus_annotation, levels = consensus_annotation))

# get list of all ontology IDs that we care about for groups 
group_ontology_ids <- consensus_groups_df |>
  dplyr::pull(validation_group_ontology) |> 
  unique()

consensus_ontology_ids <- consensus_groups_df |> 
  dplyr::pull(consensus_ontology) |> 
  unique()

all_ids <- c(group_ontology_ids, consensus_ontology_ids) |> unique()

# read in cell marker and filter to relevant ontology IDs
cell_marker_df <- readxl::read_xlsx(cellmarker_file) |> 
  dplyr::mutate(
    cellontology_id = stringr::str_replace(cellontology_id, "_", ":")
  ) |> 
  dplyr::filter(cellontology_id %in% all_ids) |> 
  # rename for later joining
  dplyr::select(cellontology_id, gene_symbol = Symbol, tissue_type) |> 
  # get rid of any thing that doesn't have a gene symbol 
  tidyr::drop_na() |> 
  unique() # only want one instance of each cell type/gene combo per tissue type 

# get ensembl ids for genes 
mapped_ids <-AnnotationDbi::mapIds(
  org.Hs.eg.db,
  keys = unique(cell_marker_df$gene_symbol),
  keytype = "SYMBOL",
  column = "ENSEMBL",
  multiVals = "first"
) |>
  tibble::enframe(name = 'gene_symbol',
                  value = 'ensembl_gene_id')

# Get top markers -------------------------------------------------------------

create_top_markers_table(
  ontology_ids = group_ontology_ids,
  annotation_column = "validation_group_annotation",
  ontology_column = "validation_group_ontology",
  celltype_groups_df = consensus_groups_df,
  marker_df = cell_marker_df
) |> 
  readr::write_tsv(validation_group_markers_file)

create_top_markers_table(
  ontology_ids = consensus_ontology_ids,
  annotation_column = "consensus_annotation",
  ontology_column = "consensus_ontology",
  celltype_groups_df = consensus_groups_df,
  marker_df = cell_marker_df
) |> 
  readr::write_tsv(consensus_celltype_markers_file)
