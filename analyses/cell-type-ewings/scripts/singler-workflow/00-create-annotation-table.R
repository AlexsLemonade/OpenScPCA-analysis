#!/usr/bin/env Rscript

# This script is used to create the reference TSV file with tumor cell annotations 
# for all samples that will be used to create a SingleR reference 
# here the reference samples being used are: SCPCL000822, SCPCL000824, SCPCL001114

# all annotation tables were created in validation notebooks found in 
# `exploratory_analysis/annotation_notebooks`

project_root <- here::here()

# Set up -----------------------------------------------------------------------

# Define file paths 
# directory with annotation tables for each library 
annotation_dir <- file.path(project_root, "results", "annotation_tables")

# annotation files 
ref_library_ids <- c("SCPCL000822", "SCPCL000824", "SCPCL001114")
ref_file_names <- glue::glue("{ref_library_ids}_tumor-classifications.tsv.gz")

file_pattern <- paste0(ref_file_names, collapse = "|")
annotation_files <- list.files(annotation_dir, pattern = file_pattern, recursive = TRUE, full.names = TRUE)

# output file 
output_annotation_file <- file.path(annotation_dir, "reference-tumor-cell-annotations.tsv")

# columns to use for tumor/ normal annotations 
annotation_columns <- c(
  "SCPCL000822" = "combined_cnv_classification", # consensus between copyKAT and InferCNV
  "SCPCL000824" = "consensus_classification", # marker gene expression cut off using SCPCL000822 + AUCell
  "SCPCL001114" = "consensus_classification" # marker gene expression cut off using SCPCL000822 + AUCell
)

# Create and export data frame -------------------------------------------------

annotation_df <- purrr::map2(annotation_files, annotation_columns, 
                             \(file, column){
                               
                               # select only the cell barcode and the desired classification
                               df <- readr::read_tsv(file) |> 
                                 dplyr::select(cell_barcode, 
                                               tumor_cell_classification = {{ column }}) |> 
                                 dplyr::mutate(
                                   # document classification method in saved table 
                                   method = stringr::str_remove(column, "_classification")
                                 )
                               
                               return(df)
                               
                             }) |> 
  # add library id as a column
  purrr::set_names(ref_library_ids) |> 
  dplyr::bind_rows(.id = "library_id")

readr::write_tsv(annotation_df, output_annotation_file)

