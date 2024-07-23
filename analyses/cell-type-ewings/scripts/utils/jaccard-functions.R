# These functions are taken directly from the supplemental cell type report produced in scpca-nf
# https://github.com/AlexsLemonade/scpca-nf/blob/main/templates/qc_report/celltypes_supplemental_report.rmd

#' Function to calculate Jaccard similarity on two vectors
#'
#' @param vec1 First vector
#' @param vec2 Second vector
#'
#' @return Jaccard similarity between the vectors
jaccard <- function(vec1, vec2) {
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
      jaccard(id1_list[[id1]], id2_list[[id2]])
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

# function that turns jaccard matrices into a list of heatmaps
make_heatmap_list <- function(jaccard_matrices, column_title, legend_match) {
  # Set heatmap padding option
  heatmap_padding <- 0.2
  ComplexHeatmap::ht_opt(TITLE_PADDING = grid::unit(heatmap_padding, "in"))
  # list of heatmaps looking at SingleR/ CellAssign vs tumor/normal
  heatmap <- jaccard_matrices |>
    purrr::imap(
      \(celltype_mat, celltype_method) {
        ComplexHeatmap::Heatmap(
          t(celltype_mat), # transpose because matrix rows are in common & we want a vertical arrangement
          col = circlize::colorRamp2(c(0, 1), colors = c("white", "darkslateblue")),
          border = TRUE,
          ## Row parameters
          cluster_rows = TRUE,
          row_title = celltype_method,
          row_title_gp = grid::gpar(fontsize = 12),
          row_title_side = "left",
          row_names_side = "left",
          row_dend_side = "right",
          row_names_gp = grid::gpar(fontsize = 10),
          ## Column parameters
          cluster_columns = TRUE,
          column_title = column_title,
          column_title_gp = grid::gpar(fontsize = 12),
          column_names_side = "bottom",
          column_names_gp = grid::gpar(fontsize = 10),
          column_names_rot = 90,
          ## Legend parameters
          heatmap_legend_param = list(
            title = "Jaccard index",
            direction = "vertical",
            legend_width = unit(1.5, "in")
          ),
          show_heatmap_legend = celltype_method == legend_match,
        )
      }
    ) |>
    # concatenate vertically into HeatmapList object
    purrr::reduce(ComplexHeatmap::`%v%`)

  return(heatmap)
}

# function to plot jaccard matrices to compare one set of annotations to other methods
# `annotation_column` will be the columns and will be compared to each method listed in `methods_to_compare`
plot_jaccard <- function(classification_df,
                         annotation_column, # single column of annotations to compare to all methods
                         methods_to_compare, # list of methods to compare to annotations
                         column_title, # title for columns
                         legend_match # what legend to keep, should be a name in `methods_to_compare`
) {
  # create jaccard matrices
  jaccard_matrices <- methods_to_compare |>
    purrr::map(\(name) {
      make_jaccard_matrix(
        classification_df,
        annotation_column,
        name
      )
    })

  heatmap <- make_heatmap_list(jaccard_matrices, column_title, legend_match)

  return(heatmap)
}
