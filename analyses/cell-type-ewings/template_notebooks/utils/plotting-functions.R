# These functions are used in `celltype-exploration.Rmd` 
# They are used for creating the summary plots in the report 


# source in helper plotting functions 
# plot_density() and plot_gene_heatmap()
repository_base <- rprojroot::find_root(rprojroot::is_git_root)
module_base <- file.path(repository_base, "analyses", "cell-type-ewings") 
validation_helper_functions <- file.path(module_base, "scripts", "utils", "tumor-validation-helpers.R")
source(validation_helper_functions)


#' Plot expression for a set of genes or gene sets on the UMAP
#' UMAP is faceted by gene or gene set and color corresponds to expression value (mean, sum, AUC, etc.)
#'
#' @param df Data frame with UMAP embeddings and gene set expression value
#' @param color_column Column to use for coloring the points in the UMAP. 
#'   This should correspond to an expression value, like counts, mean, sum, or AUC
#' @param facet_column Column to use for faceting the UMAP
#'   This should correspond to the column with a gene or gene set name
#'
#' @return UMAP faceted by facet_column and colored by color_column
#'
expression_umap <- function(
    df, 
    color_column,
    facet_column) {

  ggplot(df, aes(x = UMAP1, y = UMAP2, color = {{color_column}})) +
    geom_point(size = 0.1, alpha = 0.5) +
    scale_color_viridis_c(option = "turbo") +
    facet_wrap(vars({{facet_column}})) +
    # make sure there's a box around every facet 
    theme(
      aspect.ratio = 1,
      strip.background = element_rect(fill = "transparent", linewidth = 0.5),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
    ) +
    # remove axis numbers and background grid
    scale_x_continuous(labels = NULL, breaks = NULL) +
    scale_y_continuous(labels = NULL, breaks = NULL)
  
}


#' Density plot showing expression across clusters
#'
#' @param df Data frame with clusters for each cell and expresion values to be plotted 
#'   Must contain the `cluster` column 
#' @param expression_columns Vector of columns present in the data frame that contain some
#'   expression value to show on the x-axis of the density plot, such as AUC values or mean expression
#' @param x_label label to use for labeling all x-axes 
#'
#' @return Faceted density plot with clusters as rows and each column from expression_columns
#'   as a single panel 

cluster_density_plot <- function(
    df,
    expression_columns,
    annotation_column = "cluster",
    x_label){
  
  # create density plot for each column and combine into one figure 
  expression_columns |>
    purrr::map(\(column){
      plot_density(
        df,
        column,
        annotation_column
      ) +
        labs(x = x_label) +
        theme(text = element_text(size = 8))
    }) |>
    patchwork::wrap_plots(ncol = 2)
  
}

#' Heatmap with genes or gene sets as rows and cells as columns
#' Values in the heatmap correspond to an expression value or AUC
#' Two annotations will be included, one for cell type and one for clusters
#'
#' @param df Data frame containing exp_columns, cell_type_column, and cluster_column
#' @param exp_columns Vector of column names that contain values to be shown in the heatmap
#' @param cell_type_column Column indicating cell types for each cell in the df
#' @param cluster_column Column indicating cluster assignments for each cell in the df
#' @param legend_title Title to use for heatmap legend 
#'
#' @return Heatmap where each exp_column is a row and each cell in the df is a column

annotated_exp_heatmap <- function(
    df,
    exp_columns,
    cell_type_column,
    cluster_column,
    legend_title
){
  
  
  # get annotation colors for cell type and cluster 
  cell_types <- unique(df[[cell_type_column]])
  num_cell_types <- length(cell_types)
  cell_type_colors <- palette.colors(palette = "Dark2") |>
    head(n = num_cell_types) |>
    purrr::set_names(cell_types)
  
  clusters <- unique(df[[cluster_column]])
  num_clusters <- length(clusters)
  # use a larger palette for clusters to account for any cases with a large number of clusters 
  cluster_colors <- palette.colors(palette = "alphabet") |> 
    head(n = num_clusters) |>
    purrr::set_names(clusters)
  
  # create annotation for heatmap
  annotation <- ComplexHeatmap::columnAnnotation(
    cell_type = df[[cell_type_column]],
    cluster = df[[cluster_column]], 
    col = list(
      cell_type = cell_type_colors, 
      cluster = cluster_colors
    )
  )
  
  # build matrix for heatmap cells x gene set/ expression columns 
  heatmap_mtx <- df |>
    dplyr::select(barcodes, all_of(exp_columns)) |>
    tibble::column_to_rownames("barcodes") |>
    as.matrix() |>
    t()
  # remove any extra strings from the expression column names 
  rownames(heatmap_mtx) <- stringr::str_remove(rownames(heatmap_mtx), "auc_|_mean")
  
  # plot heatmap of marker genes
  plot_gene_heatmap(heatmap_mtx,
                    row_title = "",
                    legend_title = legend_title,
                    annotation = annotation,
                    cluster_columns = FALSE
  )
}
