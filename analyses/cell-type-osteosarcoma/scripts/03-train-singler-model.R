#!/usr/bin/env Rscript
#
# Train a SingleR model from the aggregated OsteoCAR reference
#
# The reference is created by the script '02-create-reference.R' ,
# which also updates outdated gene symbols to match ScPCA. The reference is already
# aggregated with SingleR::aggregateReference(), so no further aggregation is done here.
# An ScPCA SCE object is used to define the genes the model expects in query data:
# its rownames are converted to the gene symbols used in the OsteoCAR reference
# with rOpenScPCA::sce_to_symbols(reference = "10x2020"),
# so query SCEs must be converted the same way before classification.

suppressWarnings({
  suppressPackageStartupMessages({
    library(optparse)
    library(SingleCellExperiment)
  })
})

option_list <- list(
  make_option(
    opt_str = c("--ref_file"),
    type = "character",
    default = "",
    help = "Path to RDS file with the harmonized, aggregated OsteoCAR reference"
  ),
  make_option(
    opt_str = c("--sce_file"),
    type = "character",
    default = "",
    help = "Path to an ScPCA SCE object for determining genes to restrict to"
  ),
  make_option(
    opt_str = c("--singler_model_file"),
    type = "character",
    default = "",
    help = "Path to RDS file to save trained SingleR model"
  ),
  make_option(
    opt_str = c("--filter_genes"),
    action = "store_true",
    default = FALSE,
    help = "Whether to remove mitochondrial and ribosomal genes from the reference before training the model. Default: FALSE"
  ),
  make_option(
    opt_str = c("--threads"),
    type = "integer",
    default = 4,
    help = "Number of threads for SingleR to use"
  ),
  make_option(
    opt_str = c("--seed"),
    type = "integer",
    default = 2025,
    help = "Random seed, used for the clustering in SingleR::aggregateReference()"
  )
)

# Parse options and check arguments
opts <- parse_args(OptionParser(option_list = option_list))

stopifnot(
  "ref_file does not exist" = file.exists(opts$ref_file),
  "sce_file does not exist" = file.exists(opts$sce_file),
  "singler_model_file is required" = opts$singler_model_file != ""
)

# Set random seed for reproducibility
set.seed(opts$seed)

# Set up parallelization for SingleR
if (opts$threads > 1) {
  bp_param <- BiocParallel::MulticoreParam(opts$threads)
} else {
  bp_param <- BiocParallel::SerialParam()
}

# Read the osteocar aggregated reference
osteocar_ref <- readr::read_rds(opts$ref_file)

# Make sure the label column and logcounts assay exist in the reference
stopifnot(
  "label column is missing from the reference" = "label" %in% colnames(colData(osteocar_ref)),
  "logcounts assay is missing from the reference" = "logcounts" %in% assayNames(osteocar_ref)
)


# Gene symbols in the ScPCA SCE ------------------------------------------------

# ScPCA SCE rownames are Ensembl IDs, and the reference uses gene symbols.
# The OsteoCAR reference was built with the 10x Genomics reference
# refdata-cellranger-arc-GRCh38-2020-A-2.0.0 (Ensembl v98), while ScPCA uses
# Ensembl v104, so many gene symbols differ between the two (e.g., YARS vs. YARS1).
# To match the reference, convert the SCE rownames to the gene symbols used in the
# 10x 2020 reference with rOpenScPCA::sce_to_symbols(reference = "10x2020").
# Genes without a symbol in that reference keep their Ensembl ID.
# Query SCEs must be converted the same way before classification,
# so their rownames match `test_genes`.
test_genes <- readr::read_rds(opts$sce_file) |>
  rOpenScPCA::sce_to_symbols(reference = "10x2020") |>
  rownames()


# Define restrict vector for model training ------------------------------------

restrict_genes <- intersect(test_genes, rownames(osteocar_ref))

# if we are filtering genes, remove those from the restrict_genes list
# this way they will not be included in the model
if (opts$filter_genes) {
  # Remove mitochondrial and ribosomal genes from the restrict_genes list
  remove_genes <- c(
    grep("^MT-", restrict_genes, value = TRUE),
    grep("^RP[SL]", restrict_genes, value = TRUE)
  )
  restrict_genes <- setdiff(restrict_genes, remove_genes)
}

# Print the number of genes used for training and check if it is sufficient
message("Number of genes used for training: ", length(restrict_genes))
# a very small overlap usually means the gene identifiers don't match
stopifnot(
  "Too few genes are shared between the reference and the SCE; check gene identifiers" =
    length(restrict_genes) > 1000
)


# Train the model --------------------------------------------------------------

osteocar_trained <- SingleR::trainSingleR(
  ref = osteocar_ref,
  labels = osteocar_ref$label,
  # the aggregated reference is still fairly sparse, so wilcox is appropriate
  de.method = "wilcox",
  # the genes the query SCEs will have, after converting them to symbols with
  # rOpenScPCA::sce_to_symbols(); without this, SingleR expects the query to have
  # exactly the same genes, in the same order, as the reference
  test.genes = test_genes,
  # restrict marker selection to genes shared by the reference and the query
  restrict = restrict_genes,
  # the reference is already aggregated, make sure we set it to FALSE,
  # otherwise SingleR will aggregate it again
  aggr.ref = FALSE,
  BPPARAM = bp_param
)

# Print the labels in the trained model
message("Labels in the trained model: ", paste(sort(osteocar_trained$labels$unique),
  collapse = ", "
))

# Write the trained model to the specified file
readr::write_rds(osteocar_trained, opts$singler_model_file)
