#!/usr/bin/env Rscript
#
# Create a harmonized, aggregated OsteoCAR reference for SingleR
#
# This script takes the primary and metastasis OsteoCAR references (SCE objects)
# and creates a single harmonized reference with one level of annotation:
# - Tumor and Mesenchymal labels come from Ann_Level1
# - All other (normal) cell type labels come from Ann_Level2
# Labels are harmonized across both references, so e.g. tumor cells in the
# primary and metastasis references receive the same label.


# sets limit to 48 GB, needed for reading in the large SCE files
mem.maxVSize(48000)

library(optparse)

option_list <- list(
  make_option(
    opt_str = c("--prim_sce_file"),
    type = "character",
    help = "Path to RDS file with the SCE version of the primary OsteoCAR reference"
  ),
  make_option(
    opt_str = c("--mets_sce_file"),
    type = "character",
    help = "Path to RDS file with the SCE version of the metastasis OsteoCAR reference"
  ),
  make_option(
    opt_str = c("--output_ref_file"),
    type = "character",
    help = "Path to output RDS file to hold the harmonized, aggregated reference"
  ),
  make_option(
    opt_str = c("--ntop"),
    type = "integer",
    default = 2000,
    help = "Number of top genes to use for PCA in SingleR::aggregateReference()"
  ),
  make_option(
    opt_str = c("--threads"),
    type = "integer",
    default = 4,
    help = "Number of threads to use"
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
  "prim_sce_file does not exist" = !is.null(opts$prim_sce_file) && file.exists(opts$prim_sce_file),
  "mets_sce_file does not exist" = !is.null(opts$mets_sce_file) && file.exists(opts$mets_sce_file),
  "output_ref_file is required" = !is.null(opts$output_ref_file)
)

# load the bigger libraries after passing checks
suppressPackageStartupMessages({
  library(SingleCellExperiment)
})

set.seed(opts$seed)

if (opts$threads > 1) {
  bp_param <- BiocParallel::MulticoreParam(opts$threads)
} else {
  bp_param <- BiocParallel::SerialParam()
}

# Ann_Level1 labels to use directly; all other cells use their Ann_Level2 label
level1_labels <- c("Tumor", "Mesenchymal")


# Helper functions -------------------------------------------------------------

# Read in an SCE reference, keep only what is needed, and add `combined_label`
prepare_ref <- function(sce_file, level1_labels) {
  sce <- readr::read_rds(sce_file)

  stopifnot(
    "logcounts assay is missing from the reference" = "logcounts" %in% assayNames(sce),
    "Ann_Level1 and Ann_Level2 must be present in colData" =
      all(c("Ann_Level1", "Ann_Level2") %in% colnames(colData(sce)))
  )

  # only logcounts are needed for aggregation, so drop other assays to save memory
  assays(sce) <- assays(sce)["logcounts"]

  ann_level1 <- as.character(sce$Ann_Level1)
  ann_level2 <- as.character(sce$Ann_Level2)
  sce$combined_label <- ifelse(
    ann_level1 %in% level1_labels,
    ann_level1,
    ann_level2
  )

  # remove any cells without a label if any
  sce[, !is.na(sce$combined_label)]
}

# Aggregate a reference on `combined_label` and give it unique column names
aggregate_ref <- function(sce, ntop, bp_param) {
  ref_aggr <- SingleR::aggregateReference(
    ref = sce,
    labels = sce$combined_label,
    ntop = ntop,
    BPPARAM = bp_param
  )
}


# Read and prepare references --------------------------------------------------

prim_sce <- prepare_ref(opts$prim_sce_file, level1_labels)
mets_sce <- prepare_ref(opts$mets_sce_file, level1_labels)

# Report any labels found in only one reference.
# These may be real biological differences, but they could also be naming
# inconsistencies that need to be fixed so labels match across references.
prim_labels <- unique(prim_sce$combined_label)
mets_labels <- unique(mets_sce$combined_label)
prim_only <- setdiff(prim_labels, mets_labels)
mets_only <- setdiff(mets_labels, prim_labels)

if (length(prim_only) > 0) {
  message("Labels only in the primary reference: ", paste(prim_only, collapse = ", "))
}
if (length(mets_only) > 0) {
  message("Labels only in the metastasis reference: ", paste(mets_only, collapse = ", "))
}


# Subset to shared genes -------------------------------------------------------

shared_genes <- intersect(rownames(prim_sce), rownames(mets_sce))
message("Number of shared genes: ", length(shared_genes))

prim_sce <- prim_sce[shared_genes, ]
mets_sce <- mets_sce[shared_genes, ]
gc()


# Aggregate each reference separately ------------------------------------------

prim_aggr <- aggregate_ref(prim_sce, opts$ntop, bp_param)
rm(prim_sce)
gc()

mets_aggr <- aggregate_ref(mets_sce, opts$ntop, bp_param)
rm(mets_sce)
gc()


# Combine the aggregated references --------------------------------------------

stopifnot(
  "Gene order differs between aggregated references" =
    identical(rownames(prim_aggr), rownames(mets_aggr))
)

harmonized_ref <- cbind(prim_aggr, mets_aggr)

# Name columns using the harmonized label plus an index, so names are
# unique without distinguishing primary vs metastasis (e.g. T_cell_1, T_cell_2)
colnames(harmonized_ref) <- paste(
  harmonized_ref$label,
  seq_len(ncol(harmonized_ref)),
  sep = "_"
)

# Export the harmonized reference; use `harmonized_ref$label` as labels in SingleR
readr::write_rds(
  harmonized_ref,
  opts$output_ref_file,
  compress = "gz"
)
