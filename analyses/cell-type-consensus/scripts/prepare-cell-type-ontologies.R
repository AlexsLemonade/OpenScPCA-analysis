#!/usr/bin/env Rscript

# This script is used to match cell types present in PanglaoDB to cell ontology IDs
# The output is a TSV file that contains `ontology_id`, `human_readable_value`, and `panglao_cell_type`

# If the output file already exists any cell types with assigned ontology IDs will not be modified 
# Only cell types without ontology IDs will be used when matching to avoid any NA replacements 

# File paths -------------------------------------------------------------------
module_base <- rprojroot::find_root(rprojroot::is_renv_project)

# read in original ref file 
ref_file <- file.path(module_base, "references", "PanglaoDB_markers_2020-03-27.tsv")
ref_df <- readr::read_tsv(ref_file)

# cell ontology ref file 
ontology_output_file <- file.path(module_base, "references", "panglao-cell-type-ontologies.tsv")

# Prep ontology terms ----------------------------------------------------------

# get uberon ontology terms and ids
ol <- rols::Ontologies()
cell_ontology <- ol[["cl"]]
terms <- rols::Terms(cell_ontology)
labels <- rols::termLabel(terms) 

# data frame of id and human readable value 
label_df <- data.frame(
  ontology_id = names(labels),
  human_readable_value = labels
)

# Get existing ontology assignments --------------------------------------------

# check if ontology ID assignments already exist 
if(file.exists(ontology_output_file)){
  
  # remove any cell types without an assigned ontology
  existing_ontology_df <- readr::read_tsv(ontology_output_file) |> 
    dplyr::filter(!is.na(ontology_id))
  
  # get a list of cell types with existing ontology IDs 
  existing_ontology_cell_types <- existing_ontology_df$panglao_cell_type
  
} else {
  # if no output file just make sure these variables exist 
  existing_ontology_df <- NULL
  existing_ontology_cell_types <- NULL
}

# Match ontology terms ---------------------------------------------------------

# make human readable values for all cell types present in Panglao reference file 
cell_type_df <- ref_df |>
  dplyr::select(panglao_cell_type = `cell type`) |>
  # remove any NAs
  tidyr::drop_na(panglao_cell_type) |>
  # add a column for joining
  # generally cell type terms are in lower case and singular
  dplyr::mutate(human_readable_value = tolower(panglao_cell_type) |> 
                  # make everything singular
                  # everything is either cells or plural version of cell type with an added s 
                  stringr::str_replace("s$", "") |> 
                  # make sure B and T stay capitalized for B and T cell
                  stringr::str_replace("^b ", "B ") |> 
                  stringr::str_replace("^t ", "T ")) |>
  unique()

# join ontology terms with cell types in reference file
cell_type_terms_df <- cell_type_df |>
  # remove any cell types that already have an assigned ontology ID
  dplyr::filter(!(panglao_cell_type %in% existing_ontology_cell_types)) |> 
  # join remaining ones with ontology terms 
  dplyr::left_join(label_df, by = "human_readable_value") |>
  # make sure ontology id and human readable columns are first
  dplyr::relocate(ontology_id, human_readable_value) |>
  dplyr::mutate(human_readable_value = ifelse(is.na(ontology_id),
                                              NA_character_,
                                              human_readable_value)) |> 
  # add back existing ontologies 
  dplyr::bind_rows(existing_ontology_df) |> 
  dplyr::arrange(human_readable_value)

# export to ontology tsv
readr::write_tsv(cell_type_terms_df, ontology_output_file)
