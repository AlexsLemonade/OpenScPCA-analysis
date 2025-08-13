#!/usr/bin/env Rscript

# This script is used to match cell types present in the SCimilarity model to cell ontology IDs
# The output is a TSV file that contains `scimilarity_celltype_ontology`, `scimilarity_celltype_annotation`, and `human_readable_value`

# running this script requires downloading and unzipping the model from zenodo (https://zenodo.org/records/10685499)
# the path to the unzipped model directory should be provided with the `--model_dir` argument.
# default is `models/model_v1.1`

library(optparse)
project_root <- here::here()

option_list <- list(
  make_option(
    opt_str = c("--model_annotations_file"),
    type = "character",
    default = file.path(project_root, "models", "model_v1.1", "annotation", "reference_labels.tsv"),
    help = "Path to file containing annotations for a SCimilarity model"
  ),
  make_option(
    opt_str = c("--missing_ontology_tsv"),
    type = "character",
    default = file.path(project_root, "references", "scimilarity-missing-ontology-assignments.tsv"),
    help = "Path to TSV file with the human readable value for any scimilarity cell type annotations that
      do not directly match to cell ontology IDs"
  ),
  make_option(
    opt_str = c("--output_ontology_tsv"),
    type = "character",
    default = file.path(project_root, "references", "scimilarity-mapped-ontologies.tsv"),
    help = "Path to output TSV file with SCimilarity annotations and ontology IDs"
  )
)

opts <- parse_args(OptionParser(option_list = option_list))
# Set up -----------------------------------------------------------------------

# check that annotations exist
stopifnot(
  "--model_dir is missing annotation/reference_labels.tsv" = file.exists(opts$model_annotations_file),
  "--missing_ontology_tsv does not exist" = file.exists(opts$missing_ontology_tsv)
)

# read in labels
scimilarity_labels_df <- readr::read_tsv(opts$model_annotations_file, col_names = "scimilarity_celltype_annotation") |>
  unique()

# read in missing values
missing_df <- readr::read_tsv(opts$missing_ontology_tsv)

# Prep ontology terms ----------------------------------------------------------

# get uberon ontology terms and ids
ol <- rols::Ontologies()
cell_ontology <- ol[["cl"]]
terms <- rols::Terms(cell_ontology)
labels <- rols::termLabel(terms)

# data frame of id and human readable value
ontology_labels_df <- data.frame(
  ontology_id = names(labels),
  human_readable_value = labels
)

# Add ontology IDs to annotation labels ----------------------------------------

labels_df <- scimilarity_labels_df |>
  dplyr::left_join(missing_df, by = "scimilarity_celltype_annotation") |>
  dplyr::mutate(human_readable_value = dplyr::if_else(is.na(human_readable_value),
                                                      scimilarity_celltype_annotation,
                                                      human_readable_value)) |>
  dplyr::left_join(ontology_labels_df, by = "human_readable_value") |>
  dplyr::rename(scimilarity_celltype_ontology = ontology_id)

# check that all labels have an ontology id
stopifnot(
  "Ontology IDs cannot be found for all cell type annotations" =
    sum(is.na(labels_df$scimilarity_celltype_ontology)) == 0
)

# export
readr::write_tsv(labels_df, opts$output_ontology_tsv)

