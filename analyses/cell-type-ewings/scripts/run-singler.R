#!/usr/bin/env Rscript

# this script is used to run SingleR using a reference SCE + BlueprintEncode from celldex 
# SingleR is run three times using different references. 
# 1. Both the full BlueprintEncodeData reference from celldex and the ref SCE where cells are either annotated as tumor or normal 
# Ambiguous cells are removed from the reference 
# 2. Ref SCE where all tumor cells are labeled as tumor and all other cells keep the original annotation obtained from running SingleR as part of scpca-nf 
# 3. Ref SCE with tumor cells labeled as in #2 combined with the full BlueprintEncodeData reference 

project_root <- here::here()
renv::load(project_root)

library(optparse)
library(SingleCellExperiment)

option_list <- list(
  make_option(
    opt_str = c("--input_sce_file"),
    type = "character",
    default = NULL,
    help = "Path to RDS file containing a processed SingleCellExperiment object from scpca-nf to annotate cells in"
  ),
  make_option(
    opt_str = c("--ref_sce_file"),
    type = "character",
    default = NULL,
    help = "Path to RDS file containing a processed SingleCellExperiment object from scpca-nf to use as a reference"
  ),
  make_option(
    opt_str = c("--ref_annotations_file"),
    type = "character",
    default = NULL,
    help = "Path to file containing tumor cell annotations for reference "
  ),
  make_option(
    opt_str = c("--output_file"),
    type = "character",
    default = NULL,
    help = "Path to save SingleR annotations"
  ),
  make_option(
    opt_str = c("--scratch_dir"),
    type = "character",
    default = file.path(project_root, "scratch", "SingleR"),
    help = "path to scratch directory to store full output files from SingleR"
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

# make sure all input files exist 
stopifnot(
  "input_sce_file does not exist" = file.exists(opt$input_sce_file),
  "ref_sce_file does not exist" = file.exists(opt$ref_sce_file),
  "ref_annotations_file does not exist" = file.exists(opt$ref_annotations_file)
)

# set up multiprocessing params
if (opt$threads > 1) {
  bp_param <- BiocParallel::MulticoreParam(opt$threads)
} else {
  bp_param <- BiocParallel::SerialParam()
}

# read in sce file 
sce <- readr::read_rds(opt$input_sce_file)
library_id <- metadata(sce)$library_id

# read in ref sce and ref annotations for comparing between samples 
# ref is SCPCL000822
ref_sce <- readr::read_rds(opt$ref_sce_file)
ref_labels_df <- readr::read_tsv(opt$ref_annotations_file)


# output singler files 
fs::dir_create(opt$scratch_dir)
full_singler_output <- file.path(opt$scratch_dir, 
                                 glue::glue("{library_id}_singler-results.rds"))

# Prep references --------------------------------------------------------------

# grab HumanBlueprintEncode from celldex 
blueprint_ref <- celldex::BlueprintEncodeData(ensembl = TRUE)

# add in tumor classification for ref 
ref_coldata_df <- colData(ref_sce) |> 
  as.data.frame() |>
  dplyr::left_join(ref_labels_df, by = c("barcodes" = "cell_barcode")) |> 
  # create an updated annotation that replaces tumor cells with tumor from original SingleR annotation
  dplyr::mutate(singler_with_tumor_annotation = 
                  dplyr::if_else(
                    tumor_cell_classification == "Tumor", 
                    "Tumor", 
                    singler_celltype_annotation
                  ))

# add back to coldata 
colData(ref_sce) <- DataFrame(ref_coldata_df, row.names = ref_coldata_df$barcodes)

# remove ambiguous calls since we don't want those in the reference 
filtered_ref_sce <- ref_sce[, ref_sce$tumor_cell_classification != "Ambiguous"]

# Run SingleR ------------------------------------------------------------------

# run with tumor annotations + blueprint ref 
blueprint_tumor_results <- SingleR::SingleR(
  test = sce,
  ref = list(Blueprint = blueprint_ref,
             tumor_ref = filtered_ref_sce),
  labels = list(blueprint_ref$label.ont, filtered_ref_sce$tumor_cell_classification),
  BPPARAM = BiocParallel::MulticoreParam(4),
  restrict = rownames(sce)
)

# run with existing SingleR annotations, replacing SingleR annotations for tumor cells 
updated_singler_results <- SingleR::SingleR(
  test = sce, 
  ref = ref_sce,
  labels = ref_sce$singler_with_tumor_annotation,
  BPPARAM = BiocParallel::MulticoreParam(4),
  restrict = rownames(sce)
)

# combine updated SingleR annotations with blueprint reference 
updated_singler_blueprint_results <- SingleR::SingleR(
  test = sce, 
  ref = list(Blueprint = blueprint_ref,
             tumor_ref = ref_sce),
  labels = list(blueprint_ref$label.ont, ref_sce$singler_with_tumor_annotation),
  BPPARAM = BiocParallel::MulticoreParam(4),
  restrict = rownames(sce)
)

# Save output ------------------------------------------------------------------

ref_names <- c("blueprint_tumor", "updated_singler", "updated_singler_blueprint")

# create a list of all objects and save to rds file 
singler_results_list <- list(
  blueprint_tumor_results, 
  updated_singler_results, 
  updated_singler_blueprint_results
) |> 
  purrr::set_names(ref_names)


# get ontology labels 
cl_ont <- ontoProc::getOnto("cellOnto")
cl_df <- data.frame(
  annotation = cl_ont$name, # CL ID
  ontology = names(cl_ont$name) # human readable name 
)

# create a single df with results from all singler runs
annotations_df <- singler_results_list |> 
  purrr::imap(\(results, ref_name){

    
    # save data frame with _matching_ ontology ids and cell names
    df <- data.frame(
      barcodes = rownames(results),
      ontology = results$pruned.labels
    ) |>
      # replace ont labels with full human readable labels
      dplyr::left_join(cl_df, by = c("ontology")) |> 
      # if its not found in blueprint use the original annotation 
      dplyr::mutate(annotation = dplyr::if_else(is.na(annotation), ontology, annotation)) 
    
    colnames(df) <- c(
      "barcodes", 
      # add ref name to colnames for easier joining
      glue::glue("{ref_name}_ontology"),
      glue::glue("{ref_name}_classification") # this will make it easier to select columns with existing classifications
    )
    
    return(df)
    
  }) |> 
  purrr::reduce(dplyr::inner_join, by = "barcodes") |> 
  unique()

# save results 
readr::write_rds(singler_results_list, full_singler_output)
readr::write_tsv(annotations_df, opt$output_file)
