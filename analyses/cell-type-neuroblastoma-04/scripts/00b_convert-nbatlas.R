#!/usr/bin/env Rscript
#
# This script exports an SCE and AnnData version of a given NBAtlas Seurat object:
# - We retain only the raw counts, normalized counts, and cell metadata in the converted objects
# - Only cells present in the provided `cell_id_file` are included in the final objects
# - Gene symbols are converted to ensembl gene ids to the extent possible; as a consequence the
#   converted objects will have fewer genes compared to the original NBAtlas object

library(optparse)

option_list <- list(
  make_option(
    opt_str = c("--nbatlas_file"),
    type = "character",
    default = "",
    help = "Path to Seurat version of an NBAtlas object"
  ),
  make_option(
    opt_str = c("--tumor_metadata_file"),
    type = "character",
    default = "s",
    help = "Path to RDS file with data frame containing NBAtlas tumor metadata."
  ),
  make_option(
    opt_str = c("--cell_id_file"),
    type = "character",
    default = "",
    help = "Path to text file with cell ids present in the full atlas dated `20241203`. Any cell ids not present in this list will be excluded from the converted reference.."
  ),
  make_option(
    opt_str = c("--sce_file"),
    type = "character",
    help = "Path to output RDS file to hold an SCE version of the NBAtlas object"
  ),
  make_option(
    opt_str = c("--anndata_file"),
    type = "character",
    help = "Path to output H5AD file to hold an AnnData version of the NBAtlas object"
  ),
  make_option(
    opt_str = c("--testing"),
    action = "store_true",
    default = FALSE,
    help = "Use this flag when running in CI to output a smaller version of NBAtlas for testing"
  ),
  make_option(
    opt_str = c("--seed"),
    type = "integer",
    default = 2025,
    help = "Random seed used to subset NBAtlas when --testing is specified"
  )
)

# Parse options and check arguments
opts <- parse_args(OptionParser(option_list = option_list))

stopifnot(
  "nbatlas_file does not exist" = file.exists(opts$nbatlas_file),
  "tumor_metadata_file does not exist" = file.exists(opts$tumor_metadata_file),
  "cell_id_file does not exist" = file.exists(opts$cell_id_file),
  "sce_file was not provided" = !is.null(opts$sce_file),
  "anndata_file was not provided" = !is.null(opts$anndata_file)
)

# load the bigger libraries after passing checks
suppressPackageStartupMessages({
  library(Seurat)
  library(SingleCellExperiment)
  library(zellkonverter)
})
set.seed(opts$seed)

# read input files and determine relevant cell ids --------------
nbatlas_seurat <- readRDS(opts$nbatlas_file)
tumor_cells <- readRDS(opts$tumor_metadata_file) |>
  rownames()
all_cell_ids <- readr::read_lines(opts$cell_id_file)

# keep only cells that are present in in `all_ids` 0 -------------
nbatlas_seurat <- subset(nbatlas_seurat, cells = all_cell_ids)

# convert Seurat to SCE object directly, to save space in the final object ---------
nbatlas_sce <- SingleCellExperiment(
  assays = list(
    counts = nbatlas_seurat[["RNA"]]$counts,
    logcounts = nbatlas_seurat[["RNA"]]$data
  )
)

# convert rownames to Ensembl as much as possible -----------
# see this issue for context https://github.com/AlexsLemonade/OpenScPCA-analysis/issues/1180

# variable definitions for convenience
map_df <- rOpenScPCA::scpca_gene_reference # data frame mapping ScPCA ensembl ids to various versions of gene symbols
nbatlas_sym <- rownames(nbatlas_sce)

# define initial new vector as ensembl ids from matching scpca symbols as default, toggling in v10x2020 symbols only when the scpca symbol is undefined
scpca_ensembl <- map_df$gene_ids[match(nbatlas_sym, map_df$gene_symbol_scpca)]
v10x2020_ensembl <- map_df$gene_ids[match(nbatlas_sym, map_df$gene_symbol_10x2020)]
combined_ensembl <- ifelse(
  is.na(scpca_ensembl),
  v10x2020_ensembl,
  scpca_ensembl
)

# Now account for the symbols that are _only_ in `gene_symbol_10x2024` of which there are 4 (see linked issue)
genes_10x2024 <- nbatlas_sym[nbatlas_sym %in% map_df$gene_symbol_10x2024 &
  !(nbatlas_sym %in% c(map_df$gene_symbol_scpca, map_df$gene_symbol_10x2020))]
ensembl_10x2024 <- map_df$gene_ids[match(nbatlas_sym, genes_10x2024)]

combined_ensembl <- ifelse(
  is.na(ensembl_10x2024),
  combined_ensembl,
  ensembl_10x2024
) # ---> sum(!is.na(combined_ensembl)) ===> 31498 genes


# change over the names to Ensembl with some checks
keep_indices <- !is.na(combined_ensembl)
ensembl_rownames <- combined_ensembl[keep_indices]
stopifnot("New ensembl rownames are not actually ensembl ids" = all(stringr::str_detect(ensembl_rownames, "^ENSG\\d+$")))

nbatlas_sce <- nbatlas_sce[keep_indices, ]
rownames(nbatlas_sce) <- ensembl_rownames
stopifnot("Failed to convert gene symbols to ensembl" = nrow(nbatlas_sce) == sum(keep_indices))


# add in colData ------
# this includes updating cell labels with `neuroendocrine-tumor` label for the "tumor zoom" cells
colData(nbatlas_sce) <- nbatlas_seurat@meta.data |>
  dplyr::mutate(
    cell_id = rownames(colData(nbatlas_sce)),
    NBAtlas_label = ifelse(
      cell_id %in% tumor_cells,
      "Neuroendocrine-tumor",
      Cell_type
    )
  ) |>
  DataFrame(row.names = rownames(colData(nbatlas_sce)))

# remove Seurat file to save space -------
rm(nbatlas_seurat)
gc()


# if testing, subset the SCE to fewer cells: keep 5% of each label --------
if (opts$testing) {
  keep_cells <- colData(nbatlas_sce) |>
    as.data.frame() |>
    dplyr::group_by(NBAtlas_label) |>
    dplyr::sample_frac(0.05) |>
    dplyr::pull(cell_id)

  nbatlas_sce <- nbatlas_sce[, keep_cells]
}

# export reformatted NBAtlas objects: SCE and AnnData ----------
readr::write_rds(
  nbatlas_sce,
  opts$sce_file,
  compress = "gz"
)

zellkonverter::writeH5AD(
  nbatlas_sce,
  opts$anndata_file,
  X_name = "counts",
  compression = "gzip"
)
