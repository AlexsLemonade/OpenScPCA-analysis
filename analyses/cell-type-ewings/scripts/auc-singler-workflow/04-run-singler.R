#!/usr/bin/env Rscript

# This script runs SingleR on a processed SCE object
# All tumor cells present in a combined tumor reference is used
# as a reference along with BlueprintEncodeData from celldex
# Any cells that are from libraries from the same participant being annotated are removed from the reference


project_root <- here::here()

suppressPackageStartupMessages({
  library(optparse)
  library(SingleCellExperiment)
})

option_list <- list(
  make_option(
    opt_str = c("--sce_file"),
    type = "character",
    help = "Path to processed SCE file to annotate."
  ),
  make_option(
    opt_str = c("--tumor_reference_file"),
    type = "character",
    default = file.path(project_root, "scratch", "tumor-ref-singler.rds"),
    help = "Path to RDS file with merged object containing all tumor cells.
      The `ref_tumor_label` column in the `colData` will be used for the reference labels."
  ),
  make_option(
    opt_str = c("--output_file"),
    type = "character",
    help = "Path to TSV file to save annotations as obtained from `SingleR` for input SCE file."
  ),
  make_option(
    opt_str = c("--scratch_dir"),
    default = "scratch",
    type = "character",
    help = "Path to scratch directory to save celldex references. "
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
    default = 2024,
    help = "A random seed for reproducibility."
  )
)

# Parse options
opt <- parse_args(OptionParser(option_list = option_list))

# Set up -----------------------------------------------------------------------

# set seed
set.seed(opt$seed)

# make sure input files exist
stopifnot(
  "SCE file does not exist" = file.exists(opt$sce_file),
  "Tumor reference file does not exist" = file.exists(opt$tumor_reference_file)
)

# create output directory if it doesn't exist
output_dir <- dirname(opt$output_file)
fs::dir_create(output_dir)

# set up multiprocessing params
if (opt$threads > 1) {
  bp_param <- BiocParallel::MulticoreParam(opt$threads)
} else {
  bp_param <- BiocParallel::SerialParam()
}

# define reference files, if these exist we'll use these rather than downloading
blueprint_file <- file.path(opt$scratch_dir, "blueprint-celldex.rds")
hpca_file <- file.path(opt$scratch_dir, "hpca-celldex.rds")
cl_ont_file <- file.path(opt$scratch_dir, "cl-ont-ref.tsv")


# read in files
sce <- readr::read_rds(opt$sce_file)
tumor_ref <- readr::read_rds(opt$tumor_reference_file)

# Prep references --------------------------------------------------------------

# grab celldex references
if (!file.exists(blueprint_file)) {
  blueprint_ref <- celldex::BlueprintEncodeData(ensembl = TRUE)
  readr::read_rds(blueprint_ref, blueprint_file)
} else {
  blueprint_ref <- readr::read_rds(blueprint_file)
}

if (!file.exists(hpca_file)) {
  hpca_ref <- celldex::HumanPrimaryCellAtlasData(ensembl = TRUE)
  readr::write_tsv(hpca_ref, hpca_file)
} else {
  hpca_ref <- readr::read_rds(hpca_file)
}

# pull out library and participant id
library_id <- metadata(sce)$library_id
participant_id <- metadata(sce)$sample_metadata$participant_id

# pull out original tumor cells to combine with singler results later
sce_tumor_cells <- tumor_ref$barcodes[which(tumor_ref$library_id == library_id)]

# remove tumor cells from participant being annotated from reference
cells_to_keep <- tumor_ref$participant_id != participant_id
filtered_tumor_ref <- tumor_ref[, cells_to_keep]

# we don't need the original tumor ref anymore
rm(tumor_ref)

# Run SingleR ------------------------------------------------------------------

# run SingleR on remaining cells with tumor ref and blueprint
singler_results <- SingleR::SingleR(
  test = sce,
  ref = list(
    Blueprint = blueprint_ref,
    HPCA = hpca_ref,
    tumor_ref = filtered_tumor_ref
  ),
  labels = list(blueprint_ref$label.ont, hpca_ref$label.ont, filtered_tumor_ref$ref_tumor_label),
  BPPARAM = bp_param
)

# get ontology labels from ontoProc if we don't already have them stored
if (!file.exists(cl_ont_file)) {
  # get ontology labels
  cl_ont <- ontoProc::getOnto("cellOnto")
  cl_df <- data.frame(
    annotation = cl_ont$name, # human readable name
    ontology = names(cl_ont$name) # CL ID
  )

  readr::write_tsv(cl_df, cl_ont_file)
} else {
  cl_df <- readr::read_tsv(cl_ont_file)
}

# grab results from data frame and add human readable values for normal ont labels
results_df <- data.frame(
  barcodes = rownames(singler_results),
  ontology = singler_results$pruned.labels
) |>
  # add full human readable labels as a new column
  dplyr::left_join(cl_df, by = c("ontology")) |>
  # if its not found in blueprint its tumor so keep the original annotation
  dplyr::mutate(
    annotation = dplyr::coalesce(annotation, ontology),
    aucell_annotation = dplyr::if_else(barcodes %in% sce_tumor_cells,
      "tumor",
      "normal"
    )
  ) |>
  # specify which columns are from SingleR
  dplyr::rename(
    singler_annotation = annotation,
    singler_ontology = ontology
  )

# export tsv file with barcodes, aucell_annotation, singler_annotation and singler_ontology columns
readr::write_tsv(results_df, opt$output_file)
