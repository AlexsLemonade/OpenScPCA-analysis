#!/usr/bin/env Rscript

# this script is used to run inferCNV on a processed SCE file 

project_root <- here::here()
renv::load(project_root)

library(optparse)
library(SingleCellExperiment)

option_list <- list(
  make_option(
    opt_str = c("--sce_file"),
    type = "character",
    help = "Path to RDS file containing a processed SingleCellExperiment object from scpca-nf"
  ),
  make_option(
    opt_str = c("--annotations_file"),
    type = "character",
    help = "Path to save annotations file as tab delimited .txt file with no column headers. 
      First column is cell name and second column is annotations (normal or unknown)."
  ),
  make_option(
    opt_str = c("--gene_order_file"),
    type = "character",
    default = file.path(project_root, "references", "infercnv_refs", "Homo_sapiens.GRCh38.104.gene_order.txt"),
    help = "Path to gene order file as tab delimited .txt file with no column headers.
      Columns are: Ensembl gene id, chr, start, stop."
  ),
  make_option(
    opt_str = c("--reference_cell_file"),
    type = "character",
    default = NULL,
    help = "Optional file that contains a table of cell barcodes. 
      Any cells with Normal in the reference_cell_class column will be used as a reference for InferCNV"
  ),
  make_option(
    opt_str = c("--results_dir"),
    type = "character",
    default = NULL,
    help = "Folder within `--infercnv_results_prefix/library_id` to save infercnv results"
  ),
  make_option(
    opt_str = c("--infercnv_results_prefix"),
    type = "character",
    default = file.path(project_root, "results", "infercnv"),
    help = "path to folder where final infercnv results live"
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
    default = 1,
    help = "Number of multiprocessing threads to use."
  )
)

# Parse options
opt <- parse_args(OptionParser(option_list = option_list))


# Set up -----------------------------------------------------------------------

# make sure path to sce file exists
stopifnot("--sce_file does not exist." = file.exists(opt$sce_file))

# read in sce file
sce <- readr::read_rds(opt$sce_file)

# get library id and construct output directory
library_id <- metadata(sce)$library_id

# create output and scratch directory if not already present for library id
output_dir <- file.path(opt$infercnv_results_prefix, library_id, opt$results_dir)
fs::dir_create(output_dir)

# scratch directory organized the same way that final results are 
scratch_dir <- file.path(opt$scratch_dir, library_id, opt$results_dir)
fs::dir_create(scratch_dir)

# define output metadata file 
cnv_metadata_file <- file.path(output_dir, glue::glue("{library_id}_cnv-metadata.tsv"))

# define output infercnv obj
scratch_obj_file <- file.path(scratch_dir, "run.final.infercnv_obj")
output_obj_file <- file.path(output_dir, glue::glue("{library_id}_cnv-obj.rds"))

# png file to save 
scratch_png <- file.path(scratch_dir, "infercnv.png")
output_png <- file.path(output_dir, glue::glue("{library_id}_infercnv.png"))

# Define normal cells ----------------------------------------------------------

if(!is.null(opt$reference_cell_file)){
  
  # make sure normal cells file exists 
  stopifnot("reference file does not exist" = file.exists(opt$reference_cell_file))
  normal_cells <- readr::read_tsv(opt$reference_cell_file) |> 
    dplyr::filter(reference_cell_class == "Normal")
  
  # check that all normal cells are in colnames of sce 
  stopifnot("All barcodes in the normal_cells file are not found in the sce object" = 
              all(normal_cells %in% colnames(sce)))
  
  # define reference 
  ref_name <- "reference"
  
} else {
  # otherwise set normal cells to default 
  normal_cells <- ""
  # define reference 
  ref_name <- NULL
}

# pull out normal cells and create annotations
annotation_df <- data.frame(barcodes = colData(sce)$barcodes) |> 
  dplyr::mutate(annotations = dplyr::if_else(
    barcodes %in% normal_cells,
    "reference",
    "unknown"
  )) 

readr::write_tsv(annotation_df, opt$annotations_file, col_names = FALSE)

# Run inferCNV -----------------------------------------------------------------

# create the infercnv object
infercnv_obj <- infercnv::CreateInfercnvObject(raw_counts_matrix=counts(sce),
                                               annotations_file=opt$annotations_file,
                                               delim="\t",
                                               gene_order_file=opt$gene_order_file,
                                               ref_group_names=ref_name)

# run infercnv and report by cell rather than group
infercnv_obj <- infercnv::run(infercnv_obj,
                              cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                              out_dir=scratch_dir, # save all intermediate files to scratch dir
                              denoise=T,
                              HMM=T ,
                              save_rds = FALSE, # don't save the intermediate rds files
                              num_threads = opt$threads
)

# Save final results -----------------------------------------------------------

# create table with barcodes and CNVs for each chromosome
infercnv::add_to_seurat(seurat_obj = NULL,
                                      infercnv_output_path = scratch_dir)

# define path to output cnv metadata file saved in scratch dir
scratch_metadata_file <- file.path(scratch_dir, "map_metadata_from_infercnv.txt")

# copy metadata file to output directory 
fs::file_copy(scratch_metadata_file, cnv_metadata_file, overwrite = TRUE)

# copy final object to output directory 
fs::file_copy(scratch_obj_file, output_obj_file)

# copy png file to output directory 
fs::file_copy(scratch_png, output_png, overwrite = TRUE)

