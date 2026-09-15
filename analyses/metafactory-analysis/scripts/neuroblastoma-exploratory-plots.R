#!/usr/bin/env Rscript

# This script creates two figures used in the poster/presentation for Pediatric AACR on metafactory
# The merged object for SCPCP000004 is read in and subset to Neuroendocrine cells
# These cells are then integrated using Harmony
# Metaprogram scores for 4 programs are plot on a UMAP

# The metaprogram object is used to extract the ORA results table to make a heatmap
# This script uses the output from running Neuroblastoma samples through `ews-nf` 
# specifying only to include Neuroendocrine cells (possible tumor cells)


suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(ggplot2)
  library(patchwork)
})

# Setup ------------------------------------------------------------------------
# define input files 
# path to data directory with merged results for SCPCP000004
repository_base <- rprojroot::find_root(rprojroot::is_git_root)
data_dir <-file.path(repository_base, "data/2026-03-24/results/merge-sce")
merged_sce <- file.path(data_dir, "SCPCP000004", "SCPCP000004_merged.rds")

# path to metaprogram scores and metrics object containing ora results
module_dir <- file.path(repository_base, "analyses/metafactory-analysis")
ews_nf_out_dir <- file.path(module_dir, "results/SCPCP000004/ews-nf-output")
mp_scores_file <- file.path(ews_nf_out_dir, "k-10_metaprogram_scores.tsv.gz")
mp_obj_file <- file.path(ews_nf_out_dir, "k-10_metaprograms-metrics.rds")

# define output files
plots_dir <- here::here(module_dir, "plots/SCPCP000004")
fs::dir_create(plots_dir)
umap_file <- file.path(plots_dir, "metaprogram_scores_umap.png")
heatmap_file <- file.path(plots_dir, "ora-heatmap.pdf")

# set a seed
set.seed(2026)

# read in the merged sce
sce <- readr::read_rds(merged_sce)

# subset to the NE cells only since that's all that was used for metaprogram generation
cells_to_keep <- sce$openscpca_celltype_annotation == "Neuroendocrine"
sce <- sce[, cells_to_keep]

# read in the metaprogram scores 
mp_long <- readr::read_tsv(mp_scores_file) |>
  # add cell id column to use for joining
  dplyr::mutate(
    library_id = stringr::word(unique_id, 1, sep = "-"),
    cell_id = glue::glue("{library_id}-{barcodes}")
  ) |>
  # get rid of extra columns
  dplyr::select(-c("barcodes", "unique_id", "library_id"))

# read in mp object file
mp_metrics_obj <- readr::read_rds(mp_obj_file)


# Integration with Harmony -----------------------------------------------------
# grab the hvgs from the original object
hvg <- metadata(sce)$merged_highly_variable_gene

# run pca using the hvgs
sce <- scater::runPCA(sce, subset_row = hvg)

# run harmony with library id and participant id as covariates 
sce <- harmony::RunHarmony(
  sce,
  group.by.vars = c("library_id", "participant_id")
)

# add UMAP from Harmony embeddings 
sce <- scater::runUMAP(sce, dimred = "HARMONY", name = "UMAP_harmony")

# Add metaprograms to cell metadata for plotting -------------------------------
# pivot the scores to have one column for each metaprogram
mp_scores_df <- mp_long |>
  tidyr::pivot_wider(
    names_from = "metaprogram",
    values_from = "mp_score"
  )

# extract UMAP and join in metaprogram scores 
coldata_df <- scuttle::makePerCellDF(
  sce,
  use.dimred = "UMAP_harmony",
  features = c("cell_id")
) |>
  dplyr::left_join(
    mp_scores_df,
    by = c("cell_id")
  )

# Make the faceted plot --------------------------------------------------------
# list of metaprograms and associated names
metaprograms_to_plot <- c(
  "Translation" = "MP02", 
  "MYC" = "MP04", 
  "Differentiated" = "MP05", 
  "Cycling" = "MP06"
)

# palettes to use for each metaprogram
mp_palettes <- c(
  "MP02" = "Blues",
  "MP04" = "Greens",
  "MP05" = "Oranges",
  "MP06" = "Purples"
)

# make a faceted plot for the specified metaprograms
metaprogram_plot <- metaprograms_to_plot |>
  purrr::imap(\(mp, mp_name){

    ggplot(coldata_df,
           aes(x = UMAP_harmony.1, y = UMAP_harmony.2, color = .data[[mp]])) +
      geom_point(size = 0.01, alpha = 0.5) +
      theme_classic() +
      scale_color_distiller(palette = mp_palettes[[mp]], direction = 1) +
      labs(
        x = "UMAP1",
        y = "UMAP2",
        title = mp_name
      ) +
      theme(
        text = element_text(size = 12), 
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        aspect.ratio = 1
      )

  }) |>
  wrap_plots()

# save to a png file
ggsave(umap_file, metaprogram_plot, height = 7, width = 7)

# Heatmap of ORA results -------------------------------------------------------
# grab the ora results df
ora_results_df <- mp_metrics_obj$ora_results_df

# convert ora results into a heatmap of pvalues
heatmap_mtx <- ora_results_df |>
  dplyr::mutate(logp = -log10(p.adjust)) |>
  # only keep the top 4 gene sets
  dplyr::slice_max(logp, n = 4, by = metaprogram) |>
  dplyr::select(metaprogram, ID, logp) |>
  tidyr::pivot_wider(
    names_from = metaprogram,
    values_from = logp,
    values_fill = 0  # set missing gene sets to 0
  ) |>
  # remove _ and wrap gene set names
  dplyr::mutate(
    ID = stringr::str_replace_all(ID, "_", " ") |>
      stringr::str_to_lower() |>
      stringr::str_wrap(width = 50)
  ) |>
  tibble::column_to_rownames("ID") |>
  as.matrix()

ora_ht <- ComplexHeatmap::Heatmap(
  heatmap_mtx,
  name = "-log10(pvalue)\n",
  row_names_side = "left",
  row_dend_side = "right",
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  col = circlize::colorRamp2(c(0, 20), colors = c("gray95", "darkslateblue")),
  rect_gp = grid::gpar(col = "white", lwd = 2), # add some white lines around each box
  row_names_gp = grid::gpar(fontsize = 12, lineheight = 0.8),
  column_names_gp = grid::gpar(fontsize = 16),
  border = TRUE,
  heatmap_legend_param = list(
    title_gp = grid::gpar(fontsize = 16),
    labels_gp = grid::gpar(fontsize = 16),
    at = c(0, 5, 10, 15, 20), # custom ticks and labels to indicate it goes higher than 40
    labels = c("0", "5", "10", "15", "≥20")
  )
)

pdf(heatmap_file, width = 12, height = 12)

ora_ht |>
  # don't cut off super long gene set names
  ComplexHeatmap::draw(padding = unit(c(2, 10, 2, 2), "cm"))

dev.off()


