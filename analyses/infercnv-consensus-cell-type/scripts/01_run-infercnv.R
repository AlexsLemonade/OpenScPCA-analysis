#!/usr/bin/env Rscript
# This script runs inferCNV on an SCE object using either a pooled reference previously created, or an internal reference using cells in the input SCE

project_root <- here::here()

suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(optparse)
})

# source helper functions for preparing inferCNV annotations
source(
  project_root,
  "scripts",
  "utils.R"
)

option_list <- list(
  make_option(
    opt_str = "--sce_file",
    type = "character",
    help = "Path to the SCE file to run inferCNV on"
  ),
  make_option(
    opt_str = "--reference_type",
    type = "character",
    help = "Whether a pooled or internal normal reference should be used with inferCNV. The provided value should be one of 'pooled' or 'internal'."
  ),
  make_option(
    opt_str = "--pooled_reference_sce",
    type = "character",
    help = "When a pooled reference is used, path to the normal reference SCE. This is required if reference_type is 'pooled'."
  ),
  make_option(
    opt_str = "--internal_reference_group",
    type = "character",
    help = "When an internal reference is used, cell type grouping of the reference to use. This is required if reference_type is 'internal'."
  ),
  make_option(
    opt_str = "--celltype_tsv",
    type = "character",
    help = "Path to TSV file with consensus cell types for the library to use when building an internal reference. This is required if reference_type is 'internal'."
  ),
  make_option(
    opt_str = "--reference_celltype_tsv",
    type = "character",
    help = "TSV with cell type groups to use in the reference. This is required if reference_type is 'internal'."
  ),
  make_option(
    opt_str = c("--output_dir"),
    type = "character",
    default = "scratch",
    help = "Folder to save final infercnv results"
  ),
  make_option(
    opt_str = c("--hmm_model"),
    type = "character",
    default = "i6",
    help = "The HMM model to use with inferCNV, either 'i3' or 'i6'. Default is i6."
  ),
  make_option(
    opt_str = c("--skip_hmm"),
    action = "store_true",
    default = FALSE,
    help = "Use this flag to turn off fitting the HMM."
  ),
  # clustering options. context:
  # https://github.com/broadinstitute/infercnv/wiki/infercnv-tumor-subclusters#tumor-subclustering-by-leiden-clustering-preferred
  make_option(
    opt_str = c("--k_nn"),
    type = "numeric",
    default = 20,
    help = "Number of nearest neighbors to use during leiden clsutering in inferCNV; default is 20"
  ),
  make_option(
    opt_str = c("--leiden_method"),
    type = "character",
    default = "PCA",
    help = "Data on which clustering should be performed in inferCNV, either PCA or simple (expression data itself) for leiden clustering"
  ),
  make_option(
    opt_str = c("--leiden_function"),
    type = "character",
    default = "CPM",
    help = "Objective function to use during leiden clustering"
  ),
  make_option(
    opt_str = c("--leiden_resolution"),
    type = "character", # note this will get converted to numeric later if needed
    default = "auto",
    help = "Resolution value to use with leiden clustering in inferCNV. By default this is automatically chosen, but a numeric value can be provided.."
  ),
  make_option(
    opt_str = c("--gene_order_file"),
    type = "character",
    default = file.path(project_root, "references", "Homo_sapiens.GRCh38.104.gene_order.txt"),
    help = "Path to gene order file as tab delimited .txt file with no column headers.
      Columns are: Ensembl gene id, chr, start, stop."
  ),
  make_option(
    opt_str = "--annotation_file",
    type = "character",
    default = file.path(project_root, "scratch", "infercnv_annotation.txt"),
    help = "Path to annotation file for inferCNV input"
  ),
  make_option(
    opt_str = c("--scratch_dir"),
    type = "character",
    default = file.path(project_root, "scratch", "infercnv"),
    help = "path to scratch directory to store intermediate files from infercnv"
  ),
  make_option(
    opt_str = c("--testing"),
    action = "store_true",
    default = FALSE,
    help = "Flag to use if we are running with test data, to ensure some cells are in the reference when an internal reference is specified."
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
    default = 2025,
    help = "Random seed to set for reproducibility. Note that inferCNV is only reproducible on a given operating system."
  )
)

# Parse and check input options ----------------
opts <- parse_args(OptionParser(option_list = option_list))

stopifnot(
  "sce_file does not exist" = file.exists(opts$sce_file),
  "gene_order_file does not exist" = file.exists(opts$gene_order_file),
  "annotation_file not provided" = !is.null(opts$annotation_file),
  "output_dir was not specified" = !is.null(opts$output_dir),
  "reference_type was not specified" = !is.null(opts$reference_type),
  "hmm_model not properly specified" = opts$hmm_model %in% c("i3", "i6")
)


# check the reference information for both types of references, as different arguments are required for either type
# also, define the ref_name while here
stopifnot("reference_type must be one of 'pooled' or 'internal'." = opts$reference_type %in% c("pooled", "internal"))
if (opts$reference_type == "pooled") {
  stopifnot("pooled_reference_sce was not specified. This is required if reference_type is 'pooled'." = !is.null(opts$pooled_reference_sce))
  stopifnot("pooled_reference_sce does not exist." = file.exists(opts$pooled_reference_sce))

  ref_name <- stringr::str_remove(
    basename(opts$pooled_reference_sce),
    "^ref_"
  ) |>
    stringr::str_remove("\\.rds$")
  ref_name <- glue::glue("{ref_name}_pooled")
} else {
  stopifnot(
    "internal_reference_group was not specified. This is required if reference_type is 'internal'." = !is.null(opts$internal_reference_group),
    "celltype_tsv was not specified. This is required if reference_type is 'internal'." = !is.null(opts$celltype_tsv),
    "reference_celltype_tsv was not specified. This is required if reference_type is 'internal'." = !is.null(opts$reference_celltype_tsv)
  )
  stopifnot(
    "celltype_tsv does not exist." = file.exists(opts$celltype_tsv),
    "reference_celltype_tsv does not exist." = file.exists(opts$reference_celltype_tsv)
  )

  ref_name <- glue::glue("{opts$internal_reference_group}_internal")
}

# make the resolution parameter numeric if it's not auto
resolution <- opts$leiden_resolution
if (resolution != "auto") {
  resolution <- as.numeric(resolution)
  stopifnot(
    "The resolution parameter should either be a number or the string 'auto'." = !is.na(resolution)
  )
}

# Setup files and paths -----------------
set.seed(opts$seed)

# read SCE file
sce <- readRDS(opts$sce_file)

# get library id and construct scratch directory
library_id <- metadata(sce)$library_id

# organize scratch by library id
scratch_dir <- file.path(opts$scratch_dir, library_id)

# if scratch_dir exists, we need to remove it to prevent conflicts from previous inferCNV runs
if (dir.exists(scratch_dir)) {
  fs::dir_delete(scratch_dir)
}

# ensure directories we need exist
fs::dir_create(c(
  opts$output_dir,
  scratch_dir,
  dirname(opts$annotation_file)
))

# define output metadata file
cnv_metadata_file <- file.path(opts$output_dir, glue::glue("{library_id}_cnv-metadata.tsv"))
# define path to output cnv metadata file saved in scratch dir
scratch_metadata_file <- file.path(scratch_dir, "map_metadata_from_infercnv.txt")

# define output infercnv obj
scratch_obj_file <- file.path(scratch_dir, "run.final.infercnv_obj")
output_obj_file <- file.path(opts$output_dir, glue::glue("{library_id}_cnv-obj.rds"))

# png file to save
scratch_png <- file.path(scratch_dir, "infercnv.png")
output_png <- file.path(opts$output_dir, glue::glue("{library_id}_infercnv.png"))

# Prepare reference and annotations file  ---------------------------------------------

# Update SCE column names to always prefix with library for downstream bookkeeping
colnames(sce) <- glue::glue("{library_id}-{colnames(sce)}")

if (opts$reference_type == "pooled") {
  # Read in reference SCE file
  ref_sce <- readRDS(opts$pooled_reference_sce)

  # Remove any reference cells which also appear in this library
  duplicated_cells <- intersect(colnames(ref_sce), colnames(sce))
  sce <- sce[, !(colnames(sce) %in% duplicated_cells)]

  # Create input matrix by combining library SCE and reference SCE counts assays
  raw_counts_matrix <- cbind(
    counts(sce),
    counts(ref_sce)
  )

  # Check that we removed the duplicates:
  stopifnot("Duplicate cells present in input to inferCNV" = !all(duplicated(colnames(raw_counts_matrix))))

  # Export annotations file
  prepare_pooled_reference_annotations(
    colnames(raw_counts_matrix), # all cell ids in the matrix
    colnames(ref_sce), # reference names specifically
    opts$annotation_file # output file
  )
} else {
  # Create input matrix for inferCNV as the input counts matrix
  raw_counts_matrix <- counts(sce)

  prepare_internal_reference_annotations(
    opts$internal_reference_group, # cell type groups to include in reference
    opts$reference_celltype_tsv, # map between reference groups and consensus cell types
    opts$celltype_tsv, # consensus cell type annotations
    opts$annotation_file, # annotations file to export
    library_id, # SCE of interest library id
    opts$testing # logical if we're running with test data
  )
}

# Run infercnv ------------------------

# create the infercnv object
infercnv_obj <- infercnv::CreateInfercnvObject(
  raw_counts_matrix = raw_counts_matrix,
  annotations_file = opts$annotation_file,
  delim = "\t",
  gene_order_file = opts$gene_order_file,
  ref_group_name = "reference" # we use the label "reference" to designate reference cells
)

# Clean up some memory before running inferCNV
rm(sce, raw_counts_matrix)
if (opts$reference_type == "pooled") {
  rm(ref_sce)
}
gc()

# run infercnv
infercnv_obj <- infercnv::run(
  infercnv_obj,
  cutoff = 0.1, # use 1 for smart-seq, 0.1 for 10x-genomics
  out_dir = scratch_dir, # save all intermediate files to scratch dir
  denoise = TRUE,
  HMM = !opts$skip_hmm,
  HMM_type = opts$hmm_model, # specifies i3 or i6
  save_rds = FALSE, # don't save the intermediate rds files
  num_threads = opts$threads,
  ###### clustering settings
  leiden_function = opts$leiden_function,
  leiden_method = opts$leiden_method,
  leiden_resolution = resolution
)

# Save final results -----------------------------------------------------------

if (!opts$skip_hmm) {
  # create table with barcodes and CNVs for each chromosome
  infercnv::add_to_seurat(
    seurat_obj = NULL,
    infercnv_output_path = scratch_dir
  )

  # add reference information to metadata file and save to output directory
  # we have to read in with base R, since there are rownames
  read.table(scratch_metadata_file, header = TRUE, sep = "\t") |>
    dplyr::mutate(normal_reference = ref_name) |>
    # pull out row names into cell_id column
    tibble::rownames_to_column(var = "cell_id") |>
    readr::write_tsv(cnv_metadata_file)
}

# add reference information the options slot and save to output directory
dat <- readRDS(scratch_obj_file)
dat@options$normal_reference <- ref_name
readr::write_rds(dat, output_obj_file)

# copy png file to output directory
fs::file_copy(scratch_png, output_png, overwrite = TRUE)
