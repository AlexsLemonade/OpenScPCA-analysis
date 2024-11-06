#!/usr/bin/env Rscript

# Run `infercnv` for one sample with or without a healthy reference, where reference cells must pass the provided annotation score threshold
# infercnv
#
# USAGE:
# Rscript 06_infercnv.R \
#   --sample_id SCPCS000179
#   --reference one of none, immune, endothelium, both, pull

# OUTPUT :
# For every condition, the `infercnv` object and final `infercnv` heatmap are saved in the corresponding {sample_id}/06_infercnv results subfolder
# When running the HMM- CNV prediction, the proportion of CNV per chromosome is saved as a metadata for each cells.
# The updated `Seurat` object is additionally saved under the {sample_id} results subfolder for easy use downstream.


library(optparse)
library(Seurat)
library(infercnv)

# Parse arguments --------------------------------------------------------------
# set up arguments
option_list <- list(
  make_option(
    opt_str = c("-s", "--sample_id"),
    type = "character",
    default = "SCPCS000179",
    help = "The sample_id of the sample to be used for inference of genomic copy number using infercnv "
  ),
  make_option(
    opt_str = c("-r", "--reference"),
    type = "character",
    default = "both",
    help = "Reference cells to use as normal cells, either none, immune, endothelium, both or pull. The 'pull` option will pull in normal cells from other samples to include in the reference."
  ),
  make_option(
    opt_str = c("-m", "--HMM"),
    type = "character",
    default = "i3",
    help = "If running an additional HMM model to call CNV, either no, or i3 or i6"
  ),
  make_option(
    opt_str = c("-t", "--threshold"),
    type = "numeric",
    default = 0.85,
    help = "Threshold prediction score from label transfer to consider a normal cell in the reference"
  ),
  make_option(
    opt_str = c("--seed"),
    type = "character",
    default = "12345",
    help = "Random seed to set"
  )
)

opts <- parse_args(OptionParser(option_list = option_list))

set.seed(opts$seed)
# paths to data ----------------------------------------------------------------

# The base path for the OpenScPCA repository, found by its (hidden) .git directory
repository_base <- rprojroot::find_root(rprojroot::is_git_root)

# The path to this module
module_base <- file.path(repository_base, "analyses", "cell-type-wilms-tumor-06")

# Path to the result directory
result_dir <- file.path(module_base, "results", opts$sample_id)

# path to the pull of normal cells (inter-patient)
srat_normal_file <- file.path(module_base, "results", "references", "06b_normal-cell-reference.rds")


# path to output infercnv object
output_dir <- file.path(result_dir, "06_infercnv", glue::glue("reference-", opts$reference, "_HMM-", opts$HMM))
output_rds <- file.path(output_dir, glue::glue("06_infercnv_", opts$sample_id, "_reference-", opts$reference, "_HMM-", opts$HMM, ".rds"))
# path to heatmap png
png_file <- glue::glue("infercnv.png")
scratch_png <- file.path(output_dir, png_file)
output_png <- file.path(output_dir, glue::glue("06_infercnv_", opts$sample_id, "_reference-", opts$reference, "_HMM-", opts$HMM, "_heatmap.png"))
# path to updated seurat object
output_srat <- file.path(result_dir, glue::glue("06_infercnv_", "HMM-", opts$HMM, "_", opts$sample_id, "_reference-", opts$reference, ".rds"))

# Define functions -------------------------------------------------------------
# read_infercnv_mat will read outputs saved automatically by of infercnv in file_path
read_infercnv_mat <- function(file_path) {
  obs_table <- readr::read_delim(
    file = file_path, delim = " ", quote = '"', skip = 1,
    col_names = FALSE, progress = FALSE, show_col_types = FALSE
  )
  mat <- t(as.matrix(obs_table[, -1]))
  cell_names <- readr::read_delim(
    file = file_path, delim = " ", quote = '"', n_max = 1,
    col_names = FALSE, progress = FALSE, show_col_types = FALSE
  )
  cell_names <- as.character(cell_names[1, ])
  rownames(mat) <- cell_names
  colnames(mat) <- as.character(dplyr::pull(obs_table, 1))
  return(mat)
}


# Read in data -----------------------------------------------------------------
srat <- readRDS(
  file.path(result_dir, paste0("02b-fetal_kidney_label-transfer_", opts$sample_id, ".Rds"))
)

stopifnot("Incorrect reference provided" = opts$reference %in% c("none", "immune", "endothelium", "both", "pull"))

normal_label <- "normal" # label normal cells to use as normal

if (opts$reference %in% c("both", "endothelium", "immune")) {
  if (opts$reference == "both") {
    normal_cells <- c("endothelium", "immune")
  } else {
    normal_cells <- opts$reference
  }

  # Determine which cells pass the threshold and rename them to "normal"
  to_rename <- srat@meta.data$fetal_kidney_predicted.compartment
  names(to_rename) <- rownames(srat@meta.data)
  to_rename[to_rename %in% normal_cells &
    srat@meta.data$fetal_kidney_predicted.compartment.score > opts$threshold] <- normal_label
  srat@meta.data$fetal_kidney_predicted.compartment <- to_rename
  normal_cells <- normal_label # now they are called "normal"

  # the total count of normal cells should be be greater than 3
  total_normal <- sum(srat@meta.data$fetal_kidney_predicted.compartment == normal_label)
  stopifnot("There must be at least 3 normal cells to use a reference." = total_normal >= 3)
} else if (opts$reference == "none") {
  normal_cells <- NULL
} else { # pull
  srat_normal <- readRDS(srat_normal_file)
  # we merge the spike-in cells into the `Seurat` object
  srat <- merge(srat, srat_normal)
  srat <- JoinLayers(srat) # else GetAssayData won't work
  # we rename the `fetal_kidney_predicted.compartment` for these cells as "spike"
  to_rename <- srat@meta.data$fetal_kidney_predicted.compartment
  names(to_rename) <- rownames(srat@meta.data)
  to_rename[grepl("spike", names(to_rename))] <- normal_label
  srat@meta.data$fetal_kidney_predicted.compartment <- to_rename
  normal_cells <- normal_label
}

# Extract raw counts -----------------------------------------------------------
counts <- GetAssayData(object = srat, assay = "RNA", layer = "counts")

# Create a dataframe of annotation ---------------------------------------------
annot_df <- data.frame(condition = as.character(srat$fetal_kidney_predicted.compartment))
rownames(annot_df) <- colnames(counts)

HMM_logical <- TRUE
HMM_type <- opts$HMM

if (opts$HMM == "no") {
  HMM_logical <- FALSE
  HMM_type <- NULL
}

# create output directory if it does not exist
dir.create(output_dir, recursive = TRUE)

# retrieve the gene order file created in `06a_build-geneposition.R`
gene_order_file <- file.path(module_base, "results", "references", "gencode_v19_gene_pos.txt")

# Run infercnv ------------------------------------------------------------------
# create inferCNV object and run method
options(future.globals.maxSize = 89128960000000)
options(scipen = 100) # recommended for infercnv analysis_mode = "subclusters"


infercnv_obj <- infercnv::CreateInfercnvObject(
  raw_counts_matrix = as.matrix(counts),
  annotations_file = annot_df,
  ref_group_names = normal_cells,
  gene_order_file = gene_order_file,
  # ensure all cells are included
  min_max_counts_per_cell = c(-Inf, +Inf)
)

infercnv_obj <- infercnv::run(
  infercnv_obj,
  cutoff = 0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
  out_dir = output_dir,
  analysis_mode = "subclusters",
  cluster_by_groups = T,
  denoise = TRUE,
  HMM = HMM_logical,
  HMM_type = HMM_type,
  save_rds = TRUE,
  save_final_rds = TRUE
)

if (HMM_logical) {
  # Add `infercnv` data to the `Seurat` object
  srat <- infercnv::add_to_seurat(
    infercnv_output_path = output_dir,
    seurat_obj = srat,
    top_n = 10
  )

  # save `Seurat` object
  saveRDS(srat, output_srat)
}



# save some infercnv outputs
saveRDS(infercnv_obj, output_rds)
fs::file_copy(scratch_png, output_png, overwrite = TRUE)

# remove outputs from infercnv that we don't want to keep
files.in.dir <- list.files(output_dir, full.names = T)
files.to.keep <- c(output_png, output_rds)
files.to.remove <- list(files.in.dir[!(files.in.dir %in% files.to.keep)])
do.call(unlink, files.to.remove)
