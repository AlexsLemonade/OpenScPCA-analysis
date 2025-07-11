# This file contains functions to support making Jaccard similarity heatmaps
# Code is adapted from:
# https://github.com/AlexsLemonade/scpca-nf/blob/4bb82aa635b572a62f2028dbec587fcfd2155e26/templates/qc_report/celltypes_supplemental_report.rmd#L50

#' Function to calculate Jaccard similarity on two vectors
#'
#' @param vec1 First vector
#' @param vec2 Second vector
#'
#' @return Jaccard similarity between the vectors
calc_jaccard <- function(vec1, vec2) {
  length(intersect(vec1, vec2)) / length(union(vec1, vec2))
}


# Wrapper function to calculate jaccard similarity matrix for two categorical variables
#'
#' @param celltype_df The celltype_df data frame which must contain these columns:
#'   `colname1`, `colname2`, and `barcodes`
#' @param colname1 Column name, as a string, of first categorical variable of interest
#' @param colname2 Column name, as a string, of second categorical variable of interest
#'
#' @return Jaccard similarity matrix for the two columns. `colname1` values will
#'   be row names and `colname2` values will be column names in the final matrix
make_jaccard_matrix <- function(celltype_df, colname1, colname2) {
  # make lists of barcodes for each category, named by the category
  id1_list <- split(celltype_df$barcodes, celltype_df[[colname1]])
  id2_list <- split(celltype_df$barcodes, celltype_df[[colname2]])

  # create the grid of comparisons
  cross_df <- tidyr::expand_grid(id1 = names(id1_list), id2 = names(id2_list))

  # calculate a single Jaccard index for each combination using split lists & ids
  jaccard_scores <- cross_df |>
    purrr::pmap_dbl(\(id1, id2){
      calc_jaccard(id1_list[[id1]], id2_list[[id2]])
    })

  # add scores to the comparison grid and convert to matrix
  jaccard_matrix <- cross_df |>
    dplyr::mutate(jaccard = jaccard_scores) |>
    # convert to matrix
    tidyr::pivot_wider(
      names_from = "id2",
      values_from = "jaccard"
    ) |>
    tibble::column_to_rownames(var = "id1") |>
    as.matrix()

  return(jaccard_matrix)
}

#' Create a ComplexHeatmap from a matrix
#'
#' @param mat Matrix to create heatmap from.
#' @param row_title Label for row title.
#' @param column_title Label for column title.
#' @param labels_font_size Font size to use for rows and column labels.
#' @param keep_legend_name The name to use in the legend
#' @param col_fun Color function for the heatmap palette. Default is `heatmap_col_fun`.
#' @param ... Additional arguments to pass to `ComplexHeatmap::Heatmap()`
#'
#' @return A ComplexHeatmap object
create_single_heatmap <- function(
    mat,
    row_title,
    column_title,
    keep_legend_name,
    col_fun = heatmap_col_fun,
    ...) {
  heat <- ComplexHeatmap::Heatmap(
    t(mat), # transpose because matrix rows are in common & we want a vertical arrangement
    col = col_fun,
    border = TRUE, # each heatmap gets its own outline
    ## Row parameters
    cluster_rows = FALSE,
    row_title = row_title, # each heatmap gets its own title
    row_title_gp = grid::gpar(fontsize = 18),
    row_title_side = "right",
    row_names_side = "left",
    row_names_gp = grid::gpar(fontsize = 16),
    ## Column parameters
    cluster_columns = FALSE,
    column_title = column_title,
    column_title_gp = grid::gpar(fontsize = 18),
    column_names_side = "bottom",
    column_names_gp = grid::gpar(fontsize = 16),
    ### passed in args
    ...,
    ## Legend parameters
    heatmap_legend_param = list(
      title = "Jaccard index",
      direction = "horizontal",
      legend_width = unit(2, "in"),
      legend_gp = grid::gpar(fontsize = 12) # doesn't seem to work :/
    ),
    # only keep legends that match `keep_legend_name`
    show_heatmap_legend = row_title == keep_legend_name,
  )

  return(heat)
}





# This function calculates Jaccard similarity between two vectors and plots a heatmap
# Code is based on this scpca-nf code:
# https://github.com/AlexsLemonade/scpca-nf/blob/4bb82aa635b572a62f2028dbec587fcfd2155e26/templates/qc_report/celltypes_supplemental_report.rmd#L50
#
# df: data frame with columns of interest
# annotation_col1: first column to compare; arg should be a string
# annotation_col2: second column to compare; arg should be a string
# label1: label to use in heatmap for annotation_col1 axis
# label2: label to use in heatmap for annotation_col2 axis
make_jaccard_heatmap <- function(
    df,
    annotation_col1,
    annotation_col2,
    label1,
    label2) {
  # Calculate jaccard matrix
  jaccard_matrix <- make_jaccard_matrix(
    df,
    annotation_col1,
    annotation_col2
  )

  create_single_heatmap(
    jaccard_matrix,
    row_title = label2,
    column_title = label1,
    keep_legend_name = label2,
    col_fun = heatmap_col_fun,
    # additional arguments
    column_names_rot = 90
  ) |>
    ComplexHeatmap::draw(heatmap_legend_side = "bottom")
}
