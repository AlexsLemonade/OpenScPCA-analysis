# These functions are used in `celltype-exploration.Rmd` 
# They are used for creating the summary plots in the report 


# create faceted UMAP coloring each cell by an expression value 
expression_umap <- function(
    df, 
    color_column) {
  
  
  ggplot(df, aes(x = UMAP1, y = UMAP2, color = {{color_column}})) +
    geom_point(size = 0.1, alpha = 0.5) +
    scale_color_viridis_c() +
    facet_wrap(vars(geneset)) +
    theme(
      aspect.ratio = 1,
      strip.background = element_rect(fill = "transparent", linewidth = 0.5),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
    )
  
}

# faceted density plots 
cluster_density_plot <- function(
    df,
    expression_columns,
    label){
  
  expression_columns |>
    purrr::map(\(column){
      plot_density(
        df,
        column,
        "cluster"
      ) +
        labs(y = "Cluster", 
             x = label) +
        theme(text = element_text(size = 8))
    }) |>
    patchwork::wrap_plots(ncol = 2)
  
}

# heatmap 
annotated_exp_heatmap <- function(
    df,
    exp_columns,
    annotation_column,
    cluster_column,
    legend_title 
){
  
  
  cell_types <- unique(df[[annotation_column]])
  num_cell_types <- length(cell_types)
  cell_type_colors <- palette.colors(palette = "Dark2") |>
    head(n = num_cell_types) |>
    purrr::set_names(cell_types)
  
  clusters <- unique(df[[cluster_column]])
  num_clusters <- length(clusters)
  cluster_colors <- palette.colors(palette = "alphabet") |> 
    head(n = num_clusters) |>
    purrr::set_names(clusters)
  
  # create annotation for heatmap
  annotation <- ComplexHeatmap::columnAnnotation(
    cell_type = df[[annotation_column]],
    cluster = df[[cluster_column]], 
    col = list(
      cell_type = cell_type_colors, 
      cluster = cluster_colors
    )
  )
  
  # build matrix for heatmap cells x gene set sum or mean
  heatmap_mtx <- all_info_df |>
    dplyr::select(barcodes, all_of(auc_columns)) |>
    tibble::column_to_rownames("barcodes") |>
    as.matrix() |>
    t()
  rownames(heatmap_mtx) <- stringr::str_remove(rownames(heatmap_mtx), "auc_")
  
  # plot heatmap of marker genes
  plot_gene_heatmap(heatmap_mtx,
                    row_title = "",
                    legend_title = legend_title,
                    annotation = annotation,
                    cluster_columns = TRUE
  )
}
