#!/usr/bin/env Rscript
#
# Export SCE reference versions

# sets limit to 48 GB, needed for reading in the qs2 files
mem.maxVSize(48000)


library(optparse)

option_list <- list(
  make_option(
    opt_str = c("--input_ref_file"),
    type = "character",
    default = "",
    help = "Path to Seurat version of an OsteoCar object in qs2 format"
  ),
  make_option(
    opt_str = c("--output_sce_file"),
    type = "character",
    help = "Path to RDS file with an SCE version of the reference."
  ),
  make_option(
    opt_str = c("--output_anndata_file"),
    type = "character",
    help = "Path to output H5AD file to hold an AnnData version of the OsteoCAR object. If not provided, the AnnData object will not be saved."
  ),
  make_option(
    opt_str = c("--seed"),
    type = "integer",
    default = 2025,
    help = "Random seed used for reproducibility of selecting HVGs"
  )
)

# Parse options and check arguments
opts <- parse_args(OptionParser(option_list = option_list))

# Check that input and output arguments are valid
stopifnot("input_ref_file does not exist" = file.exists(opts$input_ref_file))


# load the bigger libraries after passing checks
suppressPackageStartupMessages({
  library(Seurat)
  library(SingleCellExperiment)
})

set.seed(opts$seed)

# read input file
osteocar_seurat <- qs2::qs_read(opts$input_ref_file)

# convert to SCE, using `as` to avoid CI error
# the clean up up for memory right away
osteocar_sce <- as.SingleCellExperiment(osteocar_seurat)
rm(osteocar_seurat)
gc()

# remove unneeded items from object
reducedDims(osteocar_sce) <- NULL

if (!is.null(opts$output_sce_file)) {
  readr::write_rds(
    osteocar_sce,
    opts$output_sce_file,
    compress = "gz"
  )
}


if (!is.null(opts$output_anndata_file)) {
  # Pare the object down to save space:
  # - subset to top 2000 HVG (batch-aware)
  # - remove logcounts
  gene_var <- scrapper::modelGeneVariances(
    logcounts(osteocar_sce),
    block = osteocar_sce$sample_name
  )
  hv_index <- scrapper::chooseHighlyVariableGenes(
    gene_var$statistics$residuals,
    top = 2000
  )
  hv_genes <- rownames(osteocar_sce)[hv_index]

  osteocar_sce <- osteocar_sce[hv_genes, ]
  logcounts(osteocar_sce) <- NULL

  # export the AnnData object, using counts as X
  anndataR::write_h5ad(
    osteocar_sce,
    opts$output_anndata_file,
    x_mapping = "counts",
    layers_mapping = FALSE,
    compression = "gzip",
    mode = "w"
  )

}
