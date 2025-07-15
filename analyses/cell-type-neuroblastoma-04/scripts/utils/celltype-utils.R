# This file contains the following functions:
# - harmonize_celltypes()
# - faceted_umap()
# - generate_dotplot()

#' Function to update NBAtlas labels
#' This function matches NBAtlas labels to consensus labels where appropriate to facilitate plots
#'
#' @param df Data frame with columns to recode
#' @param label_column Name of column with current labels
#' @param recoded_column Name of new column to create with harmonized labels
#'
#' @returns Data frame with new column as specified containing recoded labels
harmonize_celltypes <- function(df, label_column, recoded_column) {

  # groups of NBAtlas cell type labels to collapse to a broader label
  t_cells <- c("CD8+ T cell", "CD4+ T cell", "Treg", "NKT cell")
  dendritic_cells <- c("cDC2/DC3", "Migratory cDC", "cDC1")
  nk_cells <- c("Circulating NK cell", "TOX2+/KIT+ NK cell", "Resident NK cell")
  monocyte_cells <- c("Classical monocyte", "Patrolling monocyte")


  df |>
    dplyr::mutate(
      {{ recoded_column }} := dplyr::case_when(
        # nbatlas ~ consensus _broad_ annotation label
        {{ label_column }} == "Fibroblast" ~ "fibroblast",
        {{ label_column }} == "Macrophage" ~ "macrophage",
        {{ label_column }} == "Neutrophil" ~ "myeloid", # this is the only NBAtlas category that maps to our myeloid group
        {{ label_column }} == "Endothelial" ~ "endothelial cell",
        {{ label_column }} == "Plasma" ~ "plasma cell",
        {{ label_column }} %in% t_cells ~ "T cell",
        {{ label_column }} %in% dendritic_cells ~ "dendritic cell",
        {{ label_column }} %in% nk_cells ~ "natural killer cell",
        {{ label_column }} %in% monocyte_cells ~ "monocyte",
        # Update these labels for clarity, including making the stromal labels more distinguishable
        {{ label_column }} == "stromal cell" ~ "stromal cell (consensus)",
        {{ label_column }} == "Stromal other" ~ "Stromal other (NBAtlas)",
        {{ label_column }} == "NE" ~ "Neuroendocrine",
        .default = {{ label_column }}
      )
    )
}



#' Create a faceted UMAP panel where each panel has only one cell type colored
#' Adapted from scpca-nf:
#' https://github.com/AlexsLemonade/scpca-nf/blob/4bb82aa635b572a62f2028dbec587fcfd2155e26/templates/qc_report/celltypes_qc.rmd#L134-L221
#'
#' @param umap_df Data frame with UMAP1 and UMAP2 columns
#' @param annotation_column Column containing broad cell type annotation
#' @param celltype_colors Named vector of colors to use for each broader validation group
#' @param facet_type Whether to use facet_wrap or facet_grid
#' @param annotation_type_column Additional column to use if facet_type is "grid"
#'
#' @return ggplot object containing a faceted UMAP where each cell type is a facet.
#'   In each panel, the cell type of interest is colored red and all other cells are grey.
faceted_umap <- function(umap_df,
                         annotation_column,
                         celltype_colors,
                         facet_type = c("wrap", "grid"),
                         annotation_type_column = NULL) {
  facet_type <- match.arg(facet_type)

  # color by the annotation column but only color one cell type at a time
  faceted_umap <- ggplot(
    umap_df,
    aes(x = UMAP1, y = UMAP2, color = {{ annotation_column }})
  ) +
    # set points for all "other" points
    geom_point(
      data = dplyr::select(
        umap_df, -{{ annotation_column }}
      ),
      color = "gray90",
      alpha = 0.5,
      size = 0.5
    ) +
    scale_color_manual(values = celltype_colors) +
    # set points for desired cell type
    geom_point(size = 0.5, alpha = 0.5) +
    coord_equal() +
    theme_bw() +
    theme(
      legend.position = "none",
      axis.ticks = element_blank(),
      axis.text = element_blank()
    )

  # add faceting as specified
  if (facet_type == "wrap") {
    faceted_umap <- faceted_umap +
      facet_wrap(
        vars({{ annotation_column }}),
        ncol = 4
      )
  } else {
    faceted_umap <- faceted_umap +
      facet_grid(
        rows = vars({{ annotation_column }}),
        cols = vars({{ annotation_type_column }})
      )
  }
  return(faceted_umap)
}





#' Prepare and create a marker gene expression dotplot
#' Adapted from corresponding code in this ews-nf report:
#' https://github.com/AlexsLemonade/ews-nf/blob/bd38f2bd5ae581ea7dcbd98c5ef717afb9015766/templates/summary-report.Rmd
#'
#' @param merged_sce Merged SCE object to plot with a `cell_id` column in the colData
#' @param markers_df Data frame of marker genes for validation.
#' This function expects columns `marker_gene_label` (cell types), `ensembl_gene_id`, and `gene_symbol`
#' @param singler_df Data frame of singler annotations to plot.
#' This function expects columns called `label_recoded` (cell types) and `cell_id`.
#' The `cell_id` values should match the `cell_id` colData in the merged_sce
#' @param total_cells_df Data frame of cell counts and plot order.
#' This function expects columns called `label_recoded` (cell types), `y_label` (cell types with (total cells) as factor for plot order), and `total_cells`
#' @param expressed_genes Vector of genes that are expressed in the merged_sce
#' @param bar_order Vector for the annotation bar order
#' @param celltype_palette Named vector of colors for the annotation bar
#' @param min_cells Only include genes present in at least this many cells
#'
#' @returns Dotplot object
generate_dotplot <- function(
    merged_sce,
    markers_df,
    singler_df,
    total_cells_df,
    expressed_genes,
    bar_order,
    celltype_palette,
    min_cells = 0) {
  all_markers <- markers_df |>
    dplyr::pull(ensembl_gene_id) |>
    unique()

  # consider only markers that are expressed
  expressed_markers <- intersect(all_markers, expressed_genes)

  # get logcounts from merged_sce for expressed genes
  gene_exp_df <- scuttle::makePerCellDF(
    merged_sce,
    features = expressed_markers,
    assay.type = "logcounts",
    use.coldata = "cell_id",
    use.dimred = FALSE
  ) |>
    tidyr::pivot_longer(starts_with("ENSG"), names_to = "ensembl_gene_id", values_to = "logcounts") |>
    dplyr::mutate(detected = logcounts > 0)

  # Join with cell type results and marker gene info
  all_info_df <- singler_df |>
    dplyr::left_join(gene_exp_df, by = "cell_id") |>
    # account for the same gene being present in multiple cell types
    dplyr::left_join(markers_df, by = "ensembl_gene_id", relationship = "many-to-many")

  # define y- and x- axis orders

  # y axis - singler
  y_label_order <- rev(total_cells_df$y_label)

  # x axis - marker genes. this uses the passed in bar_order
  markers_df <- markers_df |>
    dplyr::filter(ensembl_gene_id %in% expressed_markers) |>
    dplyr::mutate(marker_gene_label = factor(marker_gene_label, levels = bar_order)) |>
    dplyr::arrange(marker_gene_label)

  marker_gene_order <- markers_df |>
    dplyr::pull(gene_symbol) |>
    unique()

  group_stats_df <- all_info_df |>
    # remove genes that aren't present in final annotations
    dplyr::filter(gene_symbol %in% marker_gene_order) |>
    # for each assigned cell type/marker gene combo get total detected and mean expression
    dplyr::group_by(label_recoded, marker_gene_label, ensembl_gene_id) |>
    dplyr::summarize(
      detected_count = sum(detected),
      mean_exp = mean(logcounts)
    ) |>
    dplyr::ungroup() |>
    # add total cells
    dplyr::left_join(total_cells_df, by = "label_recoded") |>
    # get total percent expressed
    dplyr::mutate(percent_exp = (detected_count / total_cells) * 100) |>
    # add in validation group for marker genes
    # this includes all possible marker genes and all possible validation group assignments
    dplyr::left_join(
      markers_df,
      by = c("ensembl_gene_id", "marker_gene_label"),
      relationship = "many-to-many"
    )

  dotplot_df <- group_stats_df |>
    dplyr::filter(
      # remove lowly expressed marker genes
      mean_exp > 0,
      percent_exp > 10,
      # only show cells that have above the minimum total number
      total_cells > min_cells
    ) |>
    # set orders of gene symbol and validation groups
    dplyr::mutate(
      y_label = factor(y_label, y_label_order),
      gene_symbol = factor(gene_symbol, levels = marker_gene_order),
      marker_gene_label = factor(marker_gene_label, levels = bar_order)
    )


  ### Make the dotplot!
  dotplot <- ggplot(dotplot_df) +
    aes(
      y = y_label,
      x = gene_symbol,
      color = mean_exp,
      size = percent_exp
    ) +
    geom_point() +
    scale_color_viridis_c(option = "magma") +
    facet_grid(cols = vars(marker_gene_label), scales = "free", space = "free") +
    theme_classic() +
    theme(
      strip.background = element_rect(fill = "transparent", color = NA),
      strip.placement = "outside",
      strip.text.x = element_blank(),
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
      axis.ticks.x = element_blank(),
      text = element_text(size = 14),
      panel.spacing = unit(0.5, "lines") # adjust spacing and match with annotation bar
    ) +
    labs(
      x = "Validation marker genes",
      y = "SingleR label",
      color = "Mean gene expression",
      size = "Percent cells expressed"
    )


  # add annotation bar aligning marker genes with validation group
  color_bar <- ggplot(dotplot_df, aes(x = gene_symbol, y = 1, fill = marker_gene_label)) +
    geom_tile() +
    facet_grid(cols = vars(marker_gene_label), scales = "free", space = "free") +
    scale_fill_manual(values = celltype_palette, breaks = levels(dotplot_df$marker_gene_label)) +
    ggmap::theme_nothing() +
    theme(
      strip.background = element_rect(fill = "transparent", color = NA),
      strip.text.x = element_text(angle = 90, hjust = 0, vjust = 0.5, size = 12),
      strip.placement = "outside",
      legend.position = "none",
      panel.spacing = unit(0.5, "lines"),
      strip.clip = "off"
    ) +
    labs(fill = "")

  combined_plot <- color_bar / dotplot +
    patchwork::plot_layout(ncol = 1, heights = c(0.1, 4))

  return(combined_plot)
}
