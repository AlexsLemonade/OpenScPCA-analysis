#!/usr/bin/env Rscript

project_root <- here::here()
renv::load(project_root)

library(optparse)
library(copykat)
library(SingleCellExperiment)

option_list <- list(
  make_option(
    opt_str = c("--sce_file"),
    type = "character",
    help = "Path to RDS file containing a processed SingleCellExperiment object from scpca-nf"
  ),
  make_option(
    opt_str = c("--singler_normal_cells"),
    type = "character",
    default = NULL,
    help = "Comma separated list of normal cell types to pull from SingleR annotations to use as reference.
      If both singler_normal_cells and cellassign_normal_cells are used, 
      the intersection of cells from both will be used as the normal reference."
  ), 
  make_option(
    opt_str = c("--cellassign_normal_cells"),
    type = "character",
    default = NULL,
    help = "Comma separated list of normal cell types to pull from CellAssign annotations to use as reference.
      If both singler_normal_cells and cellassign_normal_cells are used, 
      the intersection of cells from both will be used as the normal reference."
  ), 
  make_option(
    opt_str = c("--output_dir"),
    type = "character",
    default = file.path(project_root, "results", "copykat"),
    help = "path to folder to save all output files"
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
if(!file.exists(opt$sce_file)){
  stop("--sce_file does not exist.")
}

# read in sce file 
sce <- readr::read_rds(opt$sce_file)

# get library id and construct output directory 
library_id <- metadata(sce)$library_id

# create output directory if not already present for library id 
output_dir <- file.path(opt$output_dir, library_id)
fs::dir_create(output_dir)

# define paths to output rds files 
copykat_output <- file.path(output_dir, glue::glue("{library_id}_no_normal.rds"))
copykat_normal_output <- file.path(output_dir, glue::glue("{library_id}_with_normal.rds"))

coldata_df <- colData(sce) |> 
  as.data.frame() 

# get list of singleR cell types 
if(!is.null(opt$singler_normal_cells)){
  singler_normal_cells <- stringr::str_split(opt$singler_normal_cells, pattern = ",") |> 
    unlist()
  if(!all(singler_normal_cells %in% unique(sce$singler_celltype_annotation))){
    stop("--singler_normal_cells must be present in the singler_celltype_annotation column of the SCE object.")
  }
  
  singler_cells <- coldata_df |> 
    dplyr::filter(singler_celltype_annotation %in% singler_normal_cells) |> 
    dplyr::pull(barcodes)
    
} else {
  singler_cells <- c()
}

# get list of CelllAssign cell types 
if(!is.null(opt$cellassign_normal_cells)){
  cellassign_normal_cells <- stringr::str_split(opt$cellassign_normal_cells, pattern = ",") |> 
    unlist()
  if(!all(cellassign_normal_cells %in% unique(sce$cellassign_celltype_annotation))){
    stop("--cellassign_normal_cells must be present in the cellassign_celltype_annotation column of the SCE object.")
  }
  
  cellassign_cells <- coldata_df |> 
    dplyr::filter(cellassign_celltype_annotation %in% cellassign_normal_cells) |> 
    dplyr::pull(barcodes)
} else {
  cellassign_cells <- c()
}

# create vector of cells that should be used as normal reference 
all_normal_cells <- intersect(singler_cells, cellassign_cells)

# Run copyKat 
copykat_result <- copykat(
  rawmat = as.matrix(counts(sce)),
  id.type = "E",
  sam.name = glue::glue("{library_id}_no_normal"),
  plot.genes = FALSE,
  output.seg = FALSE,
  n.cores = opt$threads
)

# run copykat but with normal cells as reference 
copykat_normal_result <- copykat(
  rawmat = as.matrix(counts(sce)),
  id.type = "E", 
  sam.name = glue::glue("{library_id}_with_normal"),
  norm.cell.names = all_normal_cells,
  plot.genes = FALSE,
  output.seg = FALSE,
  n.cores = opt$threads
)

# save full results files 
readr::write_rds(copykat_result, copykat_output)
readr::write_rds(copykat_normal_result, copykat_normal_output)

