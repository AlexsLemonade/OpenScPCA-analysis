#!/usr/bin/env Rscript

# This script was adapted from:
# https://github.com/AlexsLemonade/OpenScPCA-analysis/blob/b1b8dea054d30f081279be9277501bd189a3c4d0/analyses/cell-type-ewings/scripts/cnv-workflow/00-make-gene-order-file.R

project_root <- here::here()
library(optparse)

option_list <- list(
  make_option(
    opt_str = c("--gtf_file"),
    type = "character",
    default = "s3://scpca-references/homo_sapiens/ensembl-104/annotation/Homo_sapiens.GRCh38.104.gtf.gz",
    help = "URI to gtf file to create gene order file."
  ),
  make_option(
    opt_str = c("--local_ref_dir"),
    type = "character",
    default = file.path(project_root, "references"),
    help = "Directory to use for saving GTF file and gene order file."
  ),
  make_option(
    opt_str = c("--scratch_dir"),
    type = "character",
    default = file.path(project_root, "scratch"),
    help = "Path to store copied GTF file"
  )
)

# Parse options
opt <- parse_args(OptionParser(option_list = option_list))

# create ref directory if doesn't already exist
fs::dir_create(opt$local_ref_dir)

# sync gtf file to local directory
gtf_filename <- basename(opt$gtf_file)

local_gtf_file <- file.path(opt$scratch_dir, gtf_filename)
if (!file.exists(local_gtf_file)) {
  sync_call <- glue::glue("aws s3 cp {opt$gtf_file} {local_gtf_file} --no-sign-request")
  system(sync_call)
}

# define gene order file name using gtf file
gtf_basename <- stringr::str_remove(gtf_filename, ".gtf.gz")
gene_order_filename <- glue::glue("{gtf_basename}.gene_order.txt")
gene_order_file <- file.path(opt$local_ref_dir, gene_order_filename)

# read in gtf file
gtf <- rtracklayer::import(local_gtf_file, feature.type = "gene")

# format gene order file
gtf_df <- gtf |>
  as.data.frame() |>
  dplyr::select(gene_id, seqnames, start, end) |>
  dplyr::mutate(seqnames = glue::glue("chr{seqnames}")) |>
  # only keep chr1 - 22 and chrX and chrY
  dplyr::filter(grepl("^chr([1-9]|1[0-9]|2[0-2]|X|Y)$", seqnames))

# export gene order file
readr::write_tsv(gtf_df, gene_order_file, col_names = FALSE)
