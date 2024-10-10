#!/usr/bin/env Rscript

# Run `infercnv` for one sample with or without a healthy reference
# infercnv
#
# USAGE:
# Rscript 06_infercnv.R \
#   --sample_id SCPCS000205
#   --reference one of none, immune, endothelium, both


library(optparse)
library(Seurat)
library(infercnv)

# Parse arguments --------------------------------------------------------------
# set up arguments
option_list <- list(
  make_option(
    opt_str = c("-s", "--sample_id"),
    type = "character",
    default = "SCPCS000205",
    help = "The sample_id of the sample to be used for inference of genomic copy number using infercnv "
  ),
  make_option(
    opt_str = c("-r", "--reference"),
    type = "character",
    default = "both",
    help = "Reference cells to use as normal cells, either none, immune, endothelium or both"
  )
)

opts <- parse_args(OptionParser(option_list = option_list))

# paths to data ----------------------------------------------------------------

# The base path for the OpenScPCA repository, found by its (hidden) .git directory
repository_base <- rprojroot::find_root(rprojroot::is_git_root)

# The path to this module
module_base <- file.path(repository_base, "analyses", "cell-type-wilms-tumor-06")

# Path to the result directory
result_dir <- file.path(module_base, "results", opts$sample_id)

# path to output infercnv object
output_dir <- file.path(result_dir,  "06_infercnv", glue::glue("reference-",opts$reference ))
output_rds <- file.path(output_dir, glue::glue("06_infercnv_",opts$sample_id,"_reference-", opts$reference, ".rds"))
# path to heatmap png
png_file <- glue::glue("infercnv.png")
scratch_png <- file.path(output_dir, png_file)
output_png <- file.path(output_dir,  glue::glue("06_infercnv_",opts$sample_id,"_reference-", opts$reference,  "_heatmap.png"))


# Define functions -------------------------------------------------------------
# read_infercnv_mat will read outputs saved automatically by of infercnv in file_path
read_infercnv_mat <- function(file_path) {
  obs_table <- readr::read_delim(file = file_path, delim = ' ', quote = '"', skip = 1,
                                 col_names = FALSE, progress = FALSE, show_col_types = FALSE)
  mat <- t(as.matrix(obs_table[, -1]))
  cell_names <- readr::read_delim(file = file_path, delim = ' ', quote = '"', n_max = 1,
                                  col_names = FALSE, progress = FALSE, show_col_types = FALSE)
  cell_names <- as.character(cell_names[1, ])
  rownames(mat) <- cell_names
  colnames(mat) <- as.character(dplyr::pull(obs_table, 1))
  return(mat)
}


# Read in data -----------------------------------------------------------------
srat <- readRDS(
  file.path(result_dir,  paste0("02b-fetal_kidney_label-transfer_",  opts$sample_id, ".Rds"))
)

# Extract raw counts -----------------------------------------------------------
counts <- GetAssayData(object = srat, assay = "RNA", layer = "counts")

# Create a dataframe of annotation ---------------------------------------------
annot_df <- data.frame(condition = as.character(srat$fetal_kidney_predicted.compartment))
rownames(annot_df) <- colnames(counts)

stopifnot("Incorrect reference provided" = opts$reference %in% c("none", "immune", "endothelium", "both"))

if(opts$reference == "both"){
  normal_cells <- c("endothelium", "immune")
} else if(opts$reference == "none"){
  normal_cells <- NULL
} else{
  normal_cells <- opts$reference
}

# create output directory if it does not exist
dir.create(output_dir, recursive = TRUE)

# retrieve the gene order file created in `06a_build-geneposition.R`
gene_order_file <- file.path(module_base, "results", "references", "gencode_v19_gene_pos.txt")

# Run infercnv ------------------------------------------------------------------
# create inferCNV object and run method

infercnv_obj <- infercnv::CreateInfercnvObject(
  raw_counts_matrix = as.matrix(counts),
  annotations_file = annot_df,
  ref_group_names = normal_cells,
  gene_order_file = gene_order_file)

infercnv_obj <- infercnv::run(
  infercnv_obj,
  cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
  out_dir=output_dir, 
  cluster_by_groups=T, 
  denoise=FALSE,
  HMM=FALSE,
  save_rds = FALSE,
  save_final_rds = FALSE
)


saveRDS(infercnv_obj, output_rds)
fs::file_copy(scratch_png, output_png, overwrite = TRUE)

# remove outputs from infercnv that we don't want to keep
files.in.dir <- list.files(output_dir, full.names = T)
files.to.keep <- c(output_png, output_rds)
files.to.remove <- list(files.in.dir[!(files.in.dir %in% files.to.keep)])
do.call(unlink, files.to.remove)
