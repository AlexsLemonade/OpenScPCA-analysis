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
