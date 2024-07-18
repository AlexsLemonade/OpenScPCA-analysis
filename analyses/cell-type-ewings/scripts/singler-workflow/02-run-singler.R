#!/usr/bin/env Rscript

# This script runs SingleR on a processed SCE object 
# All tumor cells present in a combined tumor reference is used 
# as a reference along with BlueprintEncodeData from celldex 
# Any cells that are from the library being annotated are removed from the reference 


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
    opt_str = c("--output_annotations_file"),
    type = "character",
    help = "Path to TSV file to save annotations as obtained from `SingleR` for input SCE file."
  ), 
  make_option(
    opt_str = c("-t", "--threads"),
    type = "integer",
    default = 4,
    help = "Number of multiprocessing threads to use."
  )
)

# Parse options
opt <- parse_args(OptionParser(option_list = option_list))

# Set up -----------------------------------------------------------------------

# make sure input files exist 
stopifnot(
  "SCE file does not exist" = file.exists(opt$sce_file),
  "Tumor reference file does not exist" = file.exists(opt$tumor_reference_file)
)

# create output directory if it doesn't exist
output_dir <- dirname(opt$output_annotations_file)
fs::dir_create(output_dir)

# set up multiprocessing params
if (opt$threads > 1) {
  bp_param <- BiocParallel::MulticoreParam(opt$threads)
} else {
  bp_param <- BiocParallel::SerialParam()
}

# read in files 
sce <- readr::read_rds(opt$sce_file)
tumor_ref <- readr::read_rds(opt$tumor_reference_file)

# Prep references --------------------------------------------------------------

# grab celldex reference
blueprint_ref <- celldex::BlueprintEncodeData(ensembl = TRUE)

# remove tumor cells from SCE object for annotation 
library_id <- metadata(sce)$library_id
sce_tumor_cells <- tumor_ref$barcodes[which(tumor_ref$library_id == library_id)]
filtered_sce <- sce[ , !(colnames(sce) %in% sce_tumor_cells)]

# remove tumor cells from library being annotated from reference 
cells_to_keep <- tumor_ref$library_id != library_id
filtered_tumor_ref <- tumor_ref[, cells_to_keep]

# we don't need the original tumor ref anymore
rm(tumor_ref)

# Run SingleR ------------------------------------------------------------------

# run SingleR on remaining cells with tumor ref and blueprint
singler_results <- SingleR::SingleR(
  test = filtered_sce,
  ref = list(Blueprint = blueprint_ref,
             tumor_ref = filtered_tumor_ref),
  labels = list(blueprint_ref$label.ont, filtered_tumor_ref$ref_tumor_label),
  BPPARAM = bp_param
)

# get ontology labels
cl_ont <- ontoProc::getOnto("cellOnto")
cl_df <- data.frame(
  annotation = cl_ont$name, # CL ID
  ontology = names(cl_ont$name) # human readable name
)

# grab results from data frame and add human readable values for normal ont labels 
results_df <- data.frame(
  barcodes = rownames(singler_results),
  ontology = singler_results$pruned.labels
) |>
  # replace ont labels with full human readable labels
  dplyr::left_join(cl_df, by = c("ontology")) |>
  # if its not found in blueprint its tumor so keep the original annotation 
  dplyr::mutate(annotation = dplyr::if_else(is.na(annotation), ontology, annotation))

# combine back with tumor cells 
tumor_annotations_df <- data.frame(barcodes = sce_tumor_cells) |> 
  dplyr::mutate(
  annotation = glue::glue("tumor-{library_id}")
)

all_results_df <- dplyr::bind_rows(list(results_df, tumor_annotations_df)) |> 
  # rename columns
  dplyr::rename(
    singler_tumor_annotation = annotation,
    singler_tumor_ontology = ontology
  )

# export tsv file with barcodes and singler_tumor_annotation and singler_tumor_ontology columns 
readr::write_tsv(all_results_df, opt$output_annotations_file)
