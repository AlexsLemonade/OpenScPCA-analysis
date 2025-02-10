# function to read in project data frames with all cells in a project
# output is a summarized table with total cells per sample, total cells per annotation, and number of cell types 
summarize_celltypes <- function(file, id){
  
  # read in data
  df <- readr::read_tsv(file) 
  
  # get total cell count and number of assigned cell types per library
  df <- df |> 
    dplyr::group_by(library_id) |> 
    dplyr::mutate(
      total_cells_per_library = dplyr::n(),
      num_celltypes = length(unique(consensus_annotation))
    ) |>
    dplyr::ungroup()
  
  summary_df <- df |> 
    dplyr::group_by(library_id, sample_type, consensus_annotation, consensus_ontology) |> 
    dplyr::summarize(total_cells_per_annotation = length(consensus_annotation)) |>
    dplyr::left_join(total_cells_df, by = "library_id") |> 
    dplyr::mutate(
      # add percentage 
      percent_cells_annotation = round((total_cells_per_annotation / total_cells_per_library) * 100, 2)
    ) |> 
    dplyr::ungroup()
  
  return(summary_df)
  
}
