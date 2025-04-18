#!/usr/bin/env Rscript
# This script runs inferCNV on an SCE object using the provided normal reference
#
# Usage:
#
# Rscript 01_run-infercnv.R \
#   --sce_file <input sce file> \
#   --reference_file <input normal reference file> \
#   --output_dir <output directory for infercnv results> \
#   --hmm_type <i3 or i6>
#
# For information on additional arguments, run:
#
# Rscript 01_run-infercnv.R --help
#

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
    opt_str = c("--output_dir"),
    type = "character",
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


# Parse options ----------------
opts <- parse_args(OptionParser(option_list = option_list))


# check input directories
stopifnot(
  "sce_file does not exist" = file.exists(opts$sce_file),
  "reference_file does not exist" = file.exists(opts$reference_file),
  "gene_order_file does not exist" = file.exists(opts$gene_order_file),
  "annotation_file not provided" = !is.null(opts$annotation_file),
  "output_dir was not specified" = !is.null(opts$output_dir),
  "hmm_model not properly specified" = opts$hmm_model %in% c("i3", "i6")
)


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

# Read in reference SCE file
ref_sce <- readRDS(opts$reference_file)

# define reference name to save in output
ref_name <- stringr::str_remove(
  basename(opts$reference_file),
  "^ref-"
) |>
  stringr::str_remove("\\.rds$")

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


# Prepare input data  ---------------------------------------------

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

# Create and export annotation table for inferCNV:
# "unknown" cells are uncharacterized, and "reference" cells are in the reference
data.frame(cell_id = colnames(raw_counts_matrix)) |>
  dplyr::mutate(annotations = dplyr::if_else(cell_id %in% colnames(sce), "unknown", reference_group_name)) |>
  readr::write_tsv(opts$annotation_file, col_names = FALSE)


# Run infercnv ------------------------

# create the infercnv object
infercnv_obj <- infercnv::CreateInfercnvObject(
  raw_counts_matrix = raw_counts_matrix,
  annotations_file = opts$annotation_file,
  delim = "\t",
  gene_order_file = opts$gene_order_file,
  ref_group_name = reference_group_name
)


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
