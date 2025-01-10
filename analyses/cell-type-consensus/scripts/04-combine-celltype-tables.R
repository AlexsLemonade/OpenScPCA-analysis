#!/usr/bin/env Rscript

# This script is used to combine all TSV files containing cell types into a single TSV file 
# The output TSV file will include the following added columns: 
# panglao_ontology: CL term assigned to panglao term
# panglao_annotation: human readable value associated with the CL term for panglao term 
# blueprint_annotation_fine: Fine-grained annotation from blueprint associated with singler_celltype_ontology 
# consensus_annotation: human readable name associated with the consensus label 
# consensus_ontology: CL ontology term for the consensus cell type 

project_root <- rprojroot::find_root(rprojroot::has_dir(".github"))

library(optparse)

option_list <- list(
  make_option(
    opt_str = c("--celltype_tsv_dir"),
    type = "character",
    help = "Path to directory containing TSV files with cell type annotations from single samples.
      All TSV files in this directory will be combined into a single file."
  ),
  make_option(
    opt_str = c("--panglao_ref_file"),
    default = file.path(project_root, "references", "panglao-cell-type-ontologies.tsv"),
    type = "character", 
    help = "Path to file with panglao assignments and associated cell ontology ids"
  ),
  make_option(
    opt_str = c("--consensus_ref_file"),
    default = file.path(project_root, "references", "consensus-cell-type-reference.tsv"),
    type = "character",
    help = "Path to file containing the reference for assigning consensus cell type labels"
  ), 
  make_option(
    opt_str = c("--output_file"),
    type = "character",
    help = "Path to file where combined TSV file will be saved. 
      File name must end in either `.tsv` or `.tsv.gz` to save a compressed TSV file"
  )
)

# Parse options
opt <- parse_args(OptionParser(option_list = option_list))

# Prep ref files ---------------------------------------------------------------

# make sure reference files exist 
stopifnot(
  "panglao reference file does not exist" = file.exists(opt$panglao_ref_file),
  "cell type consensus reference file does not exist" = file.exists(opt$consensus_ref_file),
  "output file must end in `.tsv` or `.tsv.gz`" = stringr::str_detect(opt$output_file, ".tsv|.tsv.gz")
)

# read in ref files 
# change names for panglao ref to match what's in the consensus file 
panglao_ref_df <- readr::read_tsv(opt$panglao_ref_file) |>
  dplyr::rename(
    panglao_ontology = ontology_id,
    panglao_annotation = human_readable_value,
    original_panglao_name = panglao_cell_type
  )

consensus_ref_df <- readr::read_tsv(opt$consensus_ref_file) |>
  # select columns to use for joining and consensus assigmments 
  dplyr::select(
    panglao_ontology, 
    original_panglao_name,
    blueprint_ontology,
    consensus_annotation,
    consensus_ontology
  )

# grab singler ref from celldex
blueprint_ref <- celldex::BlueprintEncodeData()

# get ontologies and human readable name into data frame for blueprint
# in scpca-nf we don't include the fine label so this lets us add it in 
blueprint_df <- data.frame(
  blueprint_ontology = blueprint_ref$label.ont,
  blueprint_annotation_fine = blueprint_ref$label.fine
) |>
  unique() |> 
  tidyr::drop_na()

# get list of all TSV files 
all_files <- list.files(path = opt$celltype_tsv_dir,
                        pattern = "*.tsv", 
                        full.names = TRUE) 

# read in TSV files and combine into a single df 
all_cells_df <- all_files |> 
  purrr::map(readr::read_tsv) |> 
  dplyr::bind_rows() |> 
  # add columns for panglao ontology and consensus
  # first add panglao ontology 
  dplyr::left_join(panglao_ref_df, by = c("cellassign_celltype_annotation" = "original_panglao_name")) |>
  # now add in all the blueprint columns
  dplyr::left_join(blueprint_df, by = c("singler_celltype_ontology" = "blueprint_ontology")) |> 
  # then add consensus labels
  dplyr::left_join(consensus_ref_df, 
                   by = c("singler_celltype_ontology" = "blueprint_ontology",
                          "cellassign_celltype_annotation" = "original_panglao_name",
                          "panglao_ontology")) |>
  # use unknown for NA annotation but keep ontology ID as NA
  dplyr::mutate(consensus_annotation = dplyr::if_else(is.na(consensus_annotation), "Unknown", consensus_annotation))

# export file 
readr::write_tsv(all_cells_df, opt$output_file)
