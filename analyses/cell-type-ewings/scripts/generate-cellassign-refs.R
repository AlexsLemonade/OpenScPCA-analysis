#!/usr/bin/env Rscript

# this script is used to create binary gene x cell marker genes matrices
# for use as references with CellAssign

project_root <- here::here()
renv::load(project_root)

library(optparse)

option_list <- list(
  make_option(
    opt_str = c("--local_ref_dir"),
    type = "character",
    default = file.path(project_root, "references"),
    help = "Path to where reference files live"
  ),
  make_option(
    opt_str = "--all_markers",
    type = "character",
    default = file.path(project_root, "references", "visser-all-marker-genes.tsv"),
    help = "Path to marker gene file with markers for tumor cells and all normal cells from Visser et al"
  ),
  make_option(
    opt_str = c("--panglao_file"),
    type = "character",
    default = file.path(project_root, "references", "PanglaoDB_markers_2020-03-27.tsv")
  ),
  make_option(
    opt_str = c("--gtf_file"),
    type = "character",
    default = "s3://scpca-references/homo_sapiens/ensembl-104/annotation/Homo_sapiens.GRCh38.104.gtf.gz",
    help = "URI to gtf file."
  ),
  make_option(
    opt_str = c("--output_dir"),
    type = "character",
    default = file.path(project_root, "references", "cellassign_refs"),
    help = "Directory to save all reference matrix files"
  )
)

opt <- parse_args(OptionParser(option_list = option_list))

# Functions --------------------------------------------------------------------

# function for creating a binary matrix with cell types as columns and marker genes as rows 
build_binary_df <- function(marker_genes_df) {
  
  # check column names exist 
  stopifnot("cell_type and ensembl_gene_id must be columns in marker genes file" = 
              c("cell_type", "ensembl_gene_id") %in% colnames(marker_genes_df))
  
  # build binary df 
  binary_df <- marker_genes_df |>
    dplyr::select(cell_type, ensembl_gene_id) |> 
    unique() |> 
    tidyr::pivot_wider(
      id_cols = ensembl_gene_id,
      names_from = cell_type,
      values_from = cell_type,
      values_fn = length,
      values_fill = 0
    ) |>
    # add a column with no marker genes
    # cell assign will assign cells to "other" when no other cell types are appropriate
    dplyr::mutate(other = 0) 
  
  return(binary_df)
}

# Set up -----------------------------------------------------------------------

# make sure the marker gene files exist
stopifnot("Provided marker gene file does not exist" = 
            file.exists(opt$all_markers),
          "PanglaoDB reference file does not exist" = 
            file.exists(opt$panglao_file))

# read in marker gene files panglao file 
all_markers_df <- readr::read_tsv(opt$all_markers) |> 
  # account for genes being from multiple sources
  dplyr::select(cell_type, ensembl_gene_id, gene_symbol) |> 
  dplyr::distinct()

tumor_markers_df <- all_markers_df |> 
  dplyr::filter(cell_type == "tumor")

panglao_markers_df <- readr::read_tsv(opt$panglao_file) |> 
  dplyr::filter(`cell type` %in% c("Endothelial cells", "Fibroblasts")) |> 
  dplyr::select("gene_symbol" = `official gene symbol`, "cell_type" = `cell type`)

# Grab GTF file ----------------------------------------------------------------

# sync gtf file to local directory 
gtf_filename <- stringr::word(opt$gtf_file, -1, sep = "/")

local_gtf_file <- file.path(opt$local_ref_dir, gtf_filename)
sync_call <- glue::glue('aws s3 cp {opt$gtf_file} {local_gtf_file}')
system(sync_call)

# Prep panglao df --------------------------------------------------------------

# read in gtf file (genes only for speed)
gtf <- rtracklayer::import(local_gtf_file, feature.type = "gene")

# create a data frame with ensembl id and gene symbol
gene_id_map <- gtf |>
  as.data.frame() |>
  dplyr::select(
    "ensembl_gene_id" = "gene_id",
    "gene_symbol" = "gene_name"
  ) |>
  dplyr::filter(gene_symbol != "NA") |>
  dplyr::distinct() |> 
  ## in case there are any duplicate gene_ids (there shouldn't be!)
  dplyr::group_by(ensembl_gene_id) |>
  dplyr::summarize(gene_symbol = paste(gene_symbol, collapse = ";"))

# add in gene id 
panglao_markers_df <- panglao_markers_df |> 
  dplyr::left_join(gene_id_map) |> 
  # remove symbols that don't have a match 
  dplyr::filter(!is.na(ensembl_gene_id))
  # add in the tumor markers
  dplyr::bind_rows(tumor_markers_df)

# Build and save ref tables ----------------------------------------------------

# build a list of references
marker_gene_ref_list <- list(
  "tumor-marker" = tumor_markers_df,
  "visser-all-marker" = all_markers_df,
  "panglao-endo-fibro" = panglao_markers_df
)

# turn into ref dfs and save to output files 
marker_gene_ref_list |> 
  purrr::map(build_binary_df) |> 
  purrr::iwalk(\(df, name){
    filename <- glue::glue("{name}_cellassign.tsv")
    readr::write_tsv(df, file.path(opt$output_dir, filename))
  })
  
