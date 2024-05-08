#!/usr/bin/env Rscript

project_root <- here::here()
renv::load(project_root)

library(optparse)

option_list <- list(
  make_option(
    opt_str = c("--gtf_file"),
    type = "character",
    default = file.path(project_root, "references", "Homo_sapiens.GRCh38.104.gtf.gz"),
    help = "Path to gtf file to create gene order file."
  ),
  make_option(
    opt_str = c("--gene_order_file"),
    type = "character",
    default = file.path(project_root, "references", "Homo_sapiens.GRCh38.104.gene_order.txt"),
    help = "Path to save gene order file as tab delimited .txt file with no column headers.
      Columns are: Ensembl gene id, chr, start, stop."
  )
)

# Parse options
opt <- parse_args(OptionParser(option_list = option_list))

# check that input exists 
if(!file.exists(opt$gtf_file)){
  stop("--gtf_file does not exist.")
}

# create gene order file 
gtf <- rtracklayer::import(opt$gtf_file, feature.type = "gene")

gtf_df <- gtf |> 
  as.data.frame() |> 
  dplyr::select(gene_id, seqnames, start, end) |> 
  dplyr::mutate(seqnames = glue::glue("chr{seqnames}"))

readr::write_tsv(gtf_df, opt$gene_order_file, col_names = FALSE)
