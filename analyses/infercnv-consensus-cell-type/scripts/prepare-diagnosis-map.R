#!/usr/bin/env Rscript

# This script was manually run to create the file `../references/broad-diagnosis-map.tsv`
# This script maps broad diagnosis categories to individual diagnoses in ScPCA that we plan to run inferCNV on

# Read diagnosis map -----------
diagnosis_grouping_url <- "https://raw.githubusercontent.com/AlexsLemonade/scpca-paper-figures/refs/heads/main/sample-info/diagnosis-groupings.tsv"
diagnosis_map <- readr::read_tsv(diagnosis_grouping_url)

# Prepare diagnosis map ---------
updated_diagnosis_map <- diagnosis_map |>
  # remove non-cancer
  dplyr::filter(ontology_id != "PATO:0000461") |>
  # split up leukemias based on cell types in their ontology name
  # set the "Other solid tumors" to have their ontology diagnosis be the group
  dplyr::mutate(
    diagnosis_group = dplyr::case_when(
      diagnosis_group == "Leukemia" & stringr::str_detect(human_readable_value, "T-cell") ~ "T-cell leukemia",
      diagnosis_group == "Leukemia" & stringr::str_detect(human_readable_value, "B-cell") ~ "B-cell leukemia",
      diagnosis_group == "Leukemia" & stringr::str_detect(human_readable_value, "myeloid") ~ "Acute myeloid leukemia",
      diagnosis_group == "Leukemia" & stringr::str_detect(human_readable_value, "Mixed") ~ "Mixed phenotype acute leukemia",
      diagnosis_group == "Other solid tumors" ~ human_readable_value,
      .default = diagnosis_group
    )) |>
  dplyr::arrange(diagnosis_group) |>
  # but keep leukemias together
  dplyr::arrange(stringr::str_detect(diagnosis_group, "leukemia"))
  
  

# Export diagnosis map ------------
readr::write_tsv(updated_diagnosis_map, "../references/broad-diagnosis-map.tsv")
