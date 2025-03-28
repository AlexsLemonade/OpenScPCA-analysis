#!/usr/bin/env Rscript
# This script runs inferCNV on an SCE object using the provided normal reference

project_root <- here::here()

suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(optparse)
})

option_list <- list(
  make_option(
    opt_str = "--sce_file",
    type = "character",
    help = "Path to the SCE file to run inferCNV on"
  ),
  make_option(
    opt_str = "--reference_file",
    type = "character",
    help = "Path to the normal reference SCE"
  ),
  make_option(
    opt_str = "--annotation_file",
    type = "character",
    default = file.path(project_root, "scratch", "infercnv_annotation.txt"),
    help = "Path to annotation file for inferCNV input"
  ),
  make_option(
    opt_str = c("--gene_order_file"),
    type = "character",
    default = file.path(project_root, "references", "Homo_sapiens.GRCh38.104.gene_order.txt"),
    help = "Path to gene order file as tab delimited .txt file with no column headers.
      Columns are: Ensembl gene id, chr, start, stop."
  ),
  # TODO: we may want some clustering options here later:
  # https://github.com/AlexsLemonade/OpenScPCA-analysis/issues/1089
  # https://github.com/broadinstitute/infercnv/wiki/infercnv-tumor-subclusters#tumor-subclustering-by-leiden-clustering-preferred
  make_option(
    opt_str = c("--output_dir"),
    type = "character",
    help = "Folder to save final infercnv results"
  ),
  make_option(
    opt_str = c("--scratch_dir"),
    type = "character",
    default = file.path(project_root, "scratch", "infercnv"),
    help = "path to scratch directory to store intermediate files from infercnv"
  ),
  make_option(
    opt_str = c("--hmm_model"),
    type = "character",
    default = "i3",
    help = "The HMM model to use with inferCNV, either 'i3' or 'i6'. Default is i3."
  ),
  make_option(
    opt_str = c("-t", "--threads"),
    type = "integer",
    default = 1,
    help = "Number of multiprocessing threads to use."
  ),
  make_option(
    opt_str = c("--seed"),
    type = "integer",
    default = 2025,
    help = "Random seed to set for reproducibility. Note that inferCNV is only reproducible on a given operating system."
  )
)


# Parse options ----------------
opts <- parse_args(OptionParser(option_list = option_list))


# check input directories
stopifnot(
  "sce_file does not exist" = file.exists(opts$sce_file),
  "reference_file does not exist" = file.exists(opts$reference_file),
  "gene_order_file does not exist" = file.exists(opts$gene_order_file),
  "annotation_file was not specified" = !is.null(opts$annotation_file),
  "output_dir was not specified" = !is.null(opts$output_dir),
  "scratch_dir was not specified" = !is.null(opts$scratch_dir),
  "hmm_model not properly specified" = opts$hmm_model %in% c("i3", "i6")
)

# Setup files and paths -----------------
set.seed(opts$seed)

# read SCE file so we can set up paths with its library id
sce <- readRDS(opts$sce_file)

# define reference name to set up output files and paths
ref_name <- stringr::str_remove(
  basename(opts$reference_file),
  "^ref-"
) |>
  stringr::str_remove("\\.rds$")

# get library id and construct scratch directory
library_id <- metadata(sce)$library_id

# create output and scratch directory if not already present for library id
fs::dir_create(opts$output_dir)

# organize scratch by library id
scratch_dir <- file.path(opts$scratch_dir, library_id)
fs::dir_create(scratch_dir)

# define output metadata file
cnv_metadata_file <- file.path(opts$output_dir, glue::glue("{library_id}_reference-{ref_name}_cnv-metadata.tsv"))
# define path to output cnv metadata file saved in scratch dir
scratch_metadata_file <- file.path(scratch_dir, "map_metadata_from_infercnv.txt")

# define output infercnv obj
scratch_obj_file <- file.path(scratch_dir, "run.final.infercnv_obj")
output_obj_file <- file.path(opts$output_dir, glue::glue("{library_id}_reference-{ref_name}_cnv-obj.rds"))

# png file to save
scratch_png <- file.path(scratch_dir, "infercnv.png")
output_png <- file.path(opts$output_dir, glue::glue("{library_id}_reference-{ref_name}_infercnv.png"))


# Prepare input data  ---------------------------------------------
ref_sce <- readRDS(opts$reference_file)

# Update SCE column names to always prefix with library for downstream bookkeeping
colnames(sce) <- glue::glue("{library_id}-{colnames(sce)}")

# Remove any reference cells which also appear in this library
duplicated_cells <- intersect(colnames(ref_sce), colnames(sce))
sce <- sce[, !(colnames(sce) %in% duplicated_cells)]

# Create input matrix by combining library SCE and reference SCE counts assays
raw_counts_matrix <- cbind(
  as.matrix(counts(sce)),
  as.matrix(counts(ref_sce))
)

# Check that we removed the duplicates:
stopifnot(
  "Duplicate cells present in input to inferCNV" = !all(duplicated(colnames(raw_counts_matrix)))
)

# Define reference group name for inferCNV
reference_group_name <- "reference"

# Create and export annotation table for inferCNV
data.frame(
  cell_id = colnames(sce),
  annotation = "unknown" # TODO: Should we make this a consensus annotation once that information is in the SCEs? Is this overkill?
) |>
  dplyr::bind_rows(
    data.frame(
      cell_id = colnames(ref_sce),
      annotation = reference_group_name
    )
  ) |>
  readr::write_tsv(opts$annotation_file, col_names = FALSE)


# Run infercnv ------------------------

# create the infercnv object
infercnv_obj <- infercnv::CreateInfercnvObject(
  raw_counts_matrix = raw_counts_matrix,
  annotations_file = opts$annotation_file,
  delim = "\t",
  gene_order_file = opts$gene_order_file,
  ref_group_names = reference_group_name
)


# run infercnv
infercnv_obj <- infercnv::run(
  infercnv_obj,
  cutoff = 0.1, # use 1 for smart-seq, 0.1 for 10x-genomics
  out_dir = scratch_dir, # save all intermediate files to scratch dir
  HMM = TRUE,
  HMM_type = opts$hmm_model, # specifies i3 or i6
  save_rds = FALSE, # don't save the intermediate rds files
  num_threads = opts$threads
)


# Save final results -----------------------------------------------------------

# create table with barcodes and CNVs for each chromosome
infercnv::add_to_seurat(
  seurat_obj = NULL,
  infercnv_output_path = scratch_dir
)

# format metadata file and save to output directory
fs::file_copy(scratch_metadata_file, cnv_metadata_file, overwrite = TRUE)

# copy final object to output directory
fs::file_copy(scratch_obj_file, output_obj_file, overwrite = TRUE)

# copy png file to output directory
fs::file_copy(scratch_png, output_png, overwrite = TRUE)
