#!/usr/bin/env Rscript
# This script runs inferCNV on an SCE object where reference cells are taken
# from the object itself based on a given cell type category.
#
# Usage:
#
# Rscript 01b_run-infercnv-internal-reference.R \
#   --sce_file <input sce file> \
#   --celltype_file <file with with cell-type-ewings results for this library> \
#   --reference_celltype_group <group of celltypes to include in the reference> \
#   --output_dir <output directory for infercnv results>
#
# For information on additional arguments, run:
#
# Rscript 01b_run-infercnv-internal-reference.R --help
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
    opt_str = "--reference_celltype_group",
    type = "character",
    help = "Name of a group of cell types to include in the normal reference based on `cell-type-consensus` annotations"
  ),
  make_option(
    opt_str = "--celltype_file",
    type = "character",
    help = "Path to TSV file with results from the `cell-type-ewings` module"
  ),
  make_option(
    opt_str = c("--output_dir"),
    type = "character",
    help = "Folder to save final infercnv results"
  ),
  make_option(
    opt_str = c("--testing"),
    action = "store_true",
    default = FALSE,
    help = "Flag to use if we are running with test data, to ensure some cells are in the reference group"
  ),
  make_option(
    opt_str = "--immune_ref_url",
    type = "character",
    default = "https://raw.githubusercontent.com/AlexsLemonade/OpenScPCA-analysis/refs/heads/main/analyses/cell-type-consensus/references/consensus-immune-cell-types.tsv",
    help = "URL of the OpenScPCA consensus cell type immune cells"
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


# check inputs
stopifnot(
  "sce_file does not exist" = file.exists(opts$sce_file),
  "celltype_file does not exist" = file.exists(opts$celltype_file),
  "reference_celltype_group not provided" = !is.null(opts$reference_celltype_group),
  "gene_order_file does not exist" = file.exists(opts$gene_order_file),
  "annotation_file not provided" = !is.null(opts$annotation_file),
  "output_dir was not specified" = !is.null(opts$output_dir),
  "hmm_model not properly specified" = opts$hmm_model %in% c("i3", "i6")
)

# Check the reference group
celltype_groups <- c("endo", "immune", "endo-immune")
stopifnot(
  "The reference celltype group should be one of endo, immune, or endo-immune." =
    opts$reference_celltype_group %in% celltype_groups
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

# Read in cell type TSV
celltype_df <- readr::read_tsv(opts$celltype_file)

# Define cell groups
immune_celltypes <- readr::read_tsv(opts$immune_ref_url) |>
  dplyr::pull(consensus_annotation)
endo_celltypes <- c(
  "endothelial cell",
  "blood vessel endothelial cell",
  "microvascular endothelial cell"
)
if (opts$reference_celltype_group == "immune") {
  reference_celltypes <- immune_celltypes
} else if (opts$reference_celltype_group == "endo") {
  reference_celltypes <- endo_celltypes
} else {
  reference_celltypes <- c(immune_celltypes, endo_celltypes)
}

# get library id and construct scratch directory
library_id <- metadata(sce)$library_id

# define a reference name
ref_name <- glue::glue("{opts$reference_celltype_group}_not-pooled")

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

# Prepare annotation file  ---------------------------------------------

# First, update SCE column names (and hence cell ids) to always prefix with library
# While we don't need this bookkeeping here, this will help keep data format consistent
# when running the downstream code
colnames(sce) <- glue::glue("{library_id}-{colnames(sce)}")

# Determine set of reference cells and save to annotation file
# If we're testing, assign 20% of cells to the reference. Otherwise, use cell types categories appropriately
if (opts$testing) {
  ncells <- floor(nrow(celltype_df) * 0.2)
  annotation_df <- celltype_df |>
    dplyr::mutate(
      row_index = dplyr::row_number(),
      annotation = ifelse(
        row_index <= ncells, "reference", "unknown"
      )
    )
} else {
  annotation_df <- celltype_df |>
    # add indicator for cell types intended for the reference
    dplyr::mutate(annotation = ifelse(
      consensus_annotation %in% reference_celltypes & !stringr::str_detect(ewing_annotation, "tumor"),
      "reference",
      "unknown"
    ))
}

# Export annotation TSV in expected format
annotation_df |>
  dplyr::mutate(barcodes = glue::glue("{library_id}-{barcodes}")) |>
  dplyr::select(barcodes, annotation) |>
  readr::write_tsv(opts$annotation_file, col_names = FALSE)

# Run infercnv ------------------------

# create the infercnv object
infercnv_obj <- infercnv::CreateInfercnvObject(
  raw_counts_matrix = counts(sce),
  annotations_file = opts$annotation_file,
  delim = "\t",
  gene_order_file = opts$gene_order_file,
  ref_group_name = "reference"
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
