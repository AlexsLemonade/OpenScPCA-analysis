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
    default = file.path(project_root, "results", "infercnv"),
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

# Define normal cells ----------------------------------------------------------

# make sure normal cells are defined 
if(is.null(opt$singler_normal_cells) && is.null(opt$cellassign_normal_cells)){
  stop("Must provide at least one normal cell type to use with either --singler_normal_cells or --cellassign_normal_cells")
}

coldata_df <- colData(sce) |> 
  as.data.frame()

#  get list of singleR cell types
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

# pull out normal cells and create annotations
annotation_df <- coldata_df |> 
  dplyr::mutate(annotations = dplyr::if_else(
    barcodes %in% all_normal_cells,
    "reference",
    "unknown"
  )) |> 
  dplyr::select(barcodes, annotations)

readr::write_tsv(annotation_df, opt$annotations_file, col_names = FALSE)

# Run inferCNV -----------------------------------------------------------------

# create the infercnv object
infercnv_obj <- infercnv::CreateInfercnvObject(raw_counts_matrix=counts(sce),
                                               annotations_file=opt$annotations_file,
                                               delim="\t",
                                               gene_order_file=opt$gene_order_file,
                                               ref_group_names=c("reference"))

# run infercnv and report by cell rather than group 
infercnv_obj <- infercnv::run(infercnv_obj,
                              cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                              out_dir=output_dir, 
                              denoise=T,
                              HMM=T ,
                              HMM_report_by = "cell",
                              save_rds = FALSE # don't save the intermediate rds files 
)
