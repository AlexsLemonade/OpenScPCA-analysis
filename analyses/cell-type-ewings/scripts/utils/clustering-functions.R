# These functions are used to perform clustering, calculate cluster stats,
# and produce plots related to clusters and cluster stats

# source in jaccard functions
module_base <- rprojroot::find_root(rprojroot::is_renv_project)
jaccard_functions <- file.path(module_base, "scripts", "utils", "jaccard-functions.R")
source(jaccard_functions)

# source in helper functions: plot_density()
validation_functions <- file.path(module_base, "scripts", "utils", "tumor-validation-helpers.R")
source(validation_functions)

# cluster statistics functions -------------------------------------------------


# get silhouette width and cluster purity for each cluster
# calculates values across all parameters used to determine clustering
# all_cluster_results must have cluster_params column
get_cluster_stats <- function(sce,
                              all_cluster_results) {
  pcs <- reducedDim(sce, "PCA")

  # split clustering results by param used
  split_clusters <- all_cluster_results |>
    split(all_cluster_results$cluster_params)

  # for each nn_param get cluster width and purity
  all_stats_df <- split_clusters |>
    purrr::map(\(df){
      
      # make sure there are multiple clusters present, otherwise width can't be computed
      cluster_num <- length(unique(df$cluster))
      
      if(cluster_num == 1){
        sil_df <- data.frame(
          cell_id = df$cell_id,
          cluster = df$cluster,
          other = NA,
          width = NA
        )
      } else {
        sil_df <- bluster::approxSilhouette(pcs, df$cluster) |>
          as.data.frame() |>
          tibble::rownames_to_column("cell_id") 
      }
      
      # purity can all be from one cluster 
      purity_df <- bluster::neighborPurity(pcs, df$cluster) |>
        as.data.frame() |>
        tibble::rownames_to_column("cell_id")

      # join into one data frame to return
      stats_df <- sil_df |>
        dplyr::left_join(purity_df, by = "cell_id")

      return(stats_df)
    }) |>
    dplyr::bind_rows(.id = "cluster_params") |> 
    dplyr::left_join(all_cluster_results, by = c("cell_id", "cluster_params"))

  return(all_stats_df)
}

# get cluster stability for each unique combination of params used for clustering
# must have `cluster_params` column
get_cluster_stability <- function(sce,
                                  all_cluster_results,
                                  threads) {
  pcs <- reducedDim(sce, "PCA")
  
  # split clustering results by param used
  cluster_df_list <- all_cluster_results |>
    split(all_cluster_results$cluster_params)
  
  # for each parameter, get ari values
  cluster_stability_df <- cluster_df_list |>
    purrr::map(\(df){
      
      # make sure we set objective function to available options
      objective_function <- dplyr::if_else(!is.na(unique(df$objective_function)),
                                           unique(df$objective_function),
                                           "CPM")
                                           
      
      # run stability 
      rOpenScPCA::calculate_stability(sce,
                                      cluster_df = df,
                                      algorithm = unique(df$algorithm),
                                      nn = unique(df$nn),
                                      resolution = unique(df$resolution),
                                      objective_function = objective_function,
                                      threads = threads)
      
    }) |>
    dplyr::bind_rows(.id = "cluster_params")
  
  return(cluster_stability_df)
}

# Plotting ---------------------------------------------------------------------

# plot individual stats for clusters, either purity or width
plot_cluster_stats <- function(all_stats_df,
                               stat_column,
                               plot_title) {
  
  nn_range <- unique(all_stats_df$nn)
  
  ggplot(all_stats_df, aes(x = nn, y = {{ stat_column }})) +
    # ggforce::geom_sina(size = .2) +
    ggbeeswarm::geom_quasirandom(method = "smiley", size = 0.1) +
    facet_wrap(vars(resolution),
               labeller = labeller(resolution = ~ glue::glue("{.}-res"))) +
    stat_summary(
      aes(group = nn),
      color = "red",
      # median and quartiles for point range
      fun = "median",
      fun.min = function(x) {
        quantile(x, 0.25)
      },
      fun.max = function(x) {
        quantile(x, 0.75)
      }
    ) +
    labs(
      title = plot_title
    ) +
    scale_x_continuous(breaks = nn_range)
}

# plot cluster stability 
plot_cluster_stability <- function(stat_df,
                                   plot_title){
  
  nn_range <- unique(all_stats_df$nn)
  
  ggplot(stability_df, aes(x = nn, y = ari)) +
    geom_jitter(width = 0.1) +
    facet_wrap(vars(resolution),
               labeller = labeller(resolution = ~ glue::glue("{.}-res"))) +
    labs(title = "Cluster stability") +
    stat_summary(
      aes(group = nn),
      color = "red",
      # median and quartiles for point range
      fun = "median",
      fun.min = function(x) {
        quantile(x, 0.25)
      },
      fun.max = function(x) {
        quantile(x, 0.75)
      }
    ) +
    labs(
      title = plot_title
    ) +
    scale_x_continuous(breaks = nn_range)
  
}


# heatmap comparing cluster assignments to SingleR labels
cluster_celltype_heatmap <- function(cluster_classification_df) {
  # get a jaccard mtx for each cluster param
  jaccard_df_list <- cluster_classification_df |>
    split(cluster_classification_df$nn_param) |>
    purrr::map(\(df) {
      make_jaccard_matrix(
        df,
        "cluster",
        "singler_lumped"
      )
    })

  # turn into heatmap list
  make_heatmap_list(jaccard_df_list, column_title = "clusters", legend_match = "k_5", cluster_rows = FALSE)
}


# Density plot looking at marker gene expression across all clusters
# each panel is from a different marker gene list
# each row shows marker gene expression for that cluster
plot_marker_genes <- function(cluster_exp_df,
                              k_value) {
  # pick clustering to use and select those columns
  final_clusters_df <- cluster_exp_df |>
    dplyr::filter(nn_param == k_value)

  # grab columns that contain marker gene sums
  marker_gene_columns <- colnames(final_clusters_df)[which(endsWith(colnames(final_clusters_df), "_sum"))]

  # create individual density plots and combine into one
  marker_gene_columns |>
    purrr::map(\(column){
      plot_density(
        final_clusters_df,
        column,
        "cluster"
      ) +
        labs(y = "Cluster")
    }) |>
    patchwork::wrap_plots(ncol = 2) + patchwork::plot_annotation(glue::glue("{k_value}-clusters"))
}
