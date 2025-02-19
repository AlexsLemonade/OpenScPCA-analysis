#!/usr/bin/env Rscript

# This script is used to grab the existing cell type annotations from a processed SCE object
# and assign consensus cell types based on previous annotations from SingleR and CellAssign

# The previous annotations and consensus label will be saved to a TSV file with the following annotation columns:
# singler_celltype_ontology: Original ontology label from SingleR
# singler_celltype_annotation: Original cell type name from SingleR (blueprint main label)
# cellassign_celltype_annotation: Original cell type name from CellAssign (panglao name)
# panglao_ontology: CL term assigned to panglao term
# panglao_annotation: human readable value associated with the CL term for panglao term
# blueprint_annotation_cl: human readable value associated with the CL term for singler_celltype_ontology
# consensus_annotation: human readable name associated with the consensus label
# consensus_ontology: CL ontology term for the consensus cell type

# An additional TSV file containing the expression for a set of marker genes is saved
# the `logcounts` for all genes in `ensembl_gene_id` column of the `marker_gene_file` will be saved to the output file

# output TSV columns
# library_id
# barcodes
# ensembl_gene_id
# gene_expression

library(optparse)

option_list <- list(
  make_option(
    opt_str = c("--sce_file"),
    type = "character",
    help = "Path to RDS file containing a processed SingleCellExperiment object from scpca-nf"
  ),
  make_option(
    opt_str = c("--blueprint_ref_file"),
    type = "character",
    help = "Path to file with BlueprintEncodeData cell ontology IDs and associated cell ontology names"
  ), 
  make_option(
    opt_str = c("--panglao_ref_file"),
    type = "character",
    help = "Path to file with panglao assignments and associated cell ontology ids"
  ),
  make_option(
    opt_str = c("--consensus_ref_file"),
    type = "character",
    help = "Path to file containing the reference for assigning consensus cell type labels"
  ),
  make_option(
    opt_str = c("--marker_gene_file"),
    type = "character",
    help = "Path to file containing a table of marker genes.
     The logcounts for all genes in the `ensembl_gene_id` column that are expressed in the input SCE will be saved." 
  ),
  make_option(
    opt_str = c("--consensus_output_file"),
    type = "character",
    help = "Path to file where TSV file containing consensus cell types will be saved.
      File name must end in either `.tsv` or `.tsv.gz` to save a compressed TSV file"
  ),
  make_option(
    opt_str = c("--gene_exp_output_file"),
    type = "character",
    help = "Path to file where TSV file containing marker gene expression will be saved.
      File name must end in either `.tsv` or `.tsv.gz` to save a compressed TSV file"
  )
)

# Parse options
opt <- parse_args(OptionParser(option_list = option_list))

# Set up -----------------------------------------------------------------------

# make sure input files exist
stopifnot(
  "sce file does not exist" = file.exists(opt$sce_file),
  "blueprint reference file does not exist" = file.exists(opt$blueprint_ref_file),
  "panglao reference file does not exist" = file.exists(opt$panglao_ref_file),
  "cell type consensus reference file does not exist" = file.exists(opt$consensus_ref_file),
  "marker gene file does not exist" = file.exists(opt$marker_gene_file),
  "consensus output file must end in `.tsv` or `.tsv.gz`" = stringr::str_ends(opt$consensus_output_file, "\\.tsv(\\.gz)?"),
  "gene expression output file must end in `.tsv` or `.tsv.gz`" = stringr::str_ends(opt$gene_exp_output_file, "\\.tsv(\\.gz)?")
)

# load SCE
suppressPackageStartupMessages({
  library(SingleCellExperiment)
})

# read in file
markers_df <- readr::read_tsv(opt$marker_gene_file)

#check that ensembl gene id is present
stopifnot(
  "ensembl_gene_id column is missing from marker gene file" = "ensembl_gene_id" %in% colnames(markers_df)
)

# Extract colData --------------------------------------------------------------

# read in sce 
sce <- readr::read_rds(opt$sce_file)

# extract ids 
library_id <- metadata(sce)$library_id
# account for multiplexed libraries that have multiple samples 
# for now just combine sample ids into a single string and don't worry about demultiplexing 
sample_id <- metadata(sce)$sample_id |> 
  paste0(collapse = ";")
project_id <- metadata(sce)$project_id

# check if cell line since cell lines don't have any cell type assignments
# account for having more than one sample and a list of sample types
# all sample types should be the same theoretically
sample_type <- unique(metadata(sce)$sample_type)
is_cell_line <- "cell line" %in% sample_type

# grab coldata
coldata_df <- colData(sce) |>
  as.data.frame() |>
  # add unique sample/library information
  dplyr::mutate(
    project_id = project_id,
    sample_id = sample_id,
    library_id = library_id,
    # add in sample type to make sure we don't assign consensus cell types to cell lines
    sample_type = sample_type
  )

# only select sample info and cell type info, we don't need the rest of the coldata
# if sample is cell line, fill in celltype columns with NA
if (is_cell_line) {
  celltype_df <- coldata_df |>
    dplyr::select(
      project_id,
      sample_id,
      library_id,
      barcodes,
      sample_type
    ) |>
    dplyr::mutate(
      singler_celltype_ontology = NA,
      singler_celltype_annotation = NA,
      cellassign_celltype_annotation = NA
    )
} else {
  # otherwise select the cell type columns
  celltype_df <- coldata_df |>
    dplyr::select(
      project_id,
      sample_id,
      library_id,
      barcodes,
      sample_type,
      contains("celltype") # get both singler and cellassign with ontology
    )
}


# Prep references --------------------------------------------------------------

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

# grab blueprint ref 
blueprint_df <- readr::read_tsv(opt$blueprint_ref_file)

# Create consensus TSV ---------------------------------------------------------

all_assignments_df <- celltype_df |> 
  # add columns for panglao ontology and consensus
  # first add panglao ontology
  dplyr::left_join(panglao_ref_df, by = c("cellassign_celltype_annotation" = "original_panglao_name")) |>
  # now add in all the blueprint columns
  dplyr::left_join(blueprint_df, by = c("singler_celltype_ontology" = "blueprint_ontology")) |>
  # then add consensus labels
  dplyr::left_join(
    consensus_ref_df,
    by = c(
      "singler_celltype_ontology" = "blueprint_ontology",
      "cellassign_celltype_annotation" = "original_panglao_name",
      "panglao_ontology"
    )
  ) |>
  # use unknown for NA annotation but keep ontology ID as NA
  # if the sample type is cell line, keep as NA
  dplyr::mutate(consensus_annotation = dplyr::if_else(is.na(consensus_annotation) & (sample_type != "cell line"), "Unknown", consensus_annotation))

# export file
readr::write_tsv(all_assignments_df, opt$consensus_output_file)

# Save gene expression file ----------------------------------------------------


# list of all marker genes
all_markers <- markers_df |> 
  dplyr::pull(ensembl_gene_id) |> 
  unique()

# we only care about if that gene is expressed otherwise we won't waste memory and include it
expressed_genes <- rowData(sce) |> 
  as.data.frame() |>
  dplyr::filter(detected > 0) |> 
  dplyr::pull(gene_ids)

# get markers that are expressed
expressed_markers <- intersect(all_markers, expressed_genes)

# if no markers are found, fill in columns with NA
if(length(expressed_markers) == 0)(
  
  gene_exp_df <- data.frame(
    barcodes = rownames(sce),
    library_id = library_id, 
    ensembl_gene_id = NA,
    gene_symbol = NA,
    gene_expression = NA, 
    validation_group_annotation = NA,
    validation_group_ontology = NA
  )
  
) else {
  
  # get logcounts from sce for expressed genes 
  gene_exp_df <- scuttle::makePerCellDF(
    sce, 
    features = expressed_markers, 
    assay.type = "logcounts",
    use.dimred = FALSE
  ) |>
    tidyr::pivot_longer(starts_with("ENSG"), names_to = "ensembl_gene_id", values_to = "logcounts") |> 
    # only pull out relevant columns
    dplyr::select(barcodes, ensembl_gene_id, logcounts) |> 
    # add library id column
    dplyr::mutate(library_id = library_id) |> 
    dplyr::relocate(library_id, .before = 0)
  
}

# export the file 
readr::write_tsv(gene_exp_df, opt$gene_exp_output_file)
