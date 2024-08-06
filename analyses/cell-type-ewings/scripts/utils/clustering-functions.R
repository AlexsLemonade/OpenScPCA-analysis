# These functions are used to perform clustering, calculate cluster stats,
# and produce plots related to clusters and cluster stats

# source in jaccard functions
module_base <- rprojroot::find_root(rprojroot::is_renv_project)
jaccard_functions <- file.path(module_base, "scripts", "utils", "jaccard-functions.R")
source(jaccard_functions)

# source in helper functions: plot_density()
validation_functions <- file.path(module_base, "scripts", "utils", "tumor-validation-helpers.R")
source(validation_functions)

# Perform clustering -----------------------------------------------------------

# get louvain, jaccard clusters for a specified value of k (nearest neighbors)
get_clusters <- function(pcs, k) {
  clusters <- bluster::clusterRows(
    pcs,
    bluster::NNGraphParam(
      k = k,
      type = "jaccard",
      cluster.fun = "louvain"
    )
  )

  return(clusters)
}

# define a function to perform clustersweep and get clusters across multiple values of k (5,40,5)
cluster_sweep <- function(sce) {
  # first perform clustering across parameters
  cluster_results <- bluster::clusterSweep(reducedDim(sce, "PCA"),
    bluster::NNGraphParam(),
    k = as.integer(seq(5, 40, 5)),
    cluster.fun = "louvain",
    type = "jaccard"
  )

  # turn results into a data frame
  cluster_df <- cluster_results$clusters |>
    as.data.frame() |>
    # add barcode column
    dplyr::mutate(barcodes = colnames(sce)) |>
    # combine all cluster results into one column
    tidyr::pivot_longer(
      cols = ends_with("jaccard"),
      names_to = "params",
      values_to = "cluster"
    ) |>
    # separate out parameters, nn, function, and type into their own columns
    dplyr::mutate(
      nn_param = stringr::word(params, 1, sep = "_") |>
        stringr::str_replace("k.", "k_"),
      cluster_fun = stringr::word(params, 2, sep = "_") |>
        stringr::str_remove("cluster.fun."),
      cluster_type = stringr::word(params, -1, sep = "_") |>
        stringr::str_remove("type.")
    ) |>
    # remove combined params column
    dplyr::select(-params)

  return(cluster_df)
}

# cluster statistics functions -------------------------------------------------


# get silhouette width and cluster purity for each cluster
# calculates values across all nn_param options used to determine clustering
# all_cluster_results must have nn_param column
get_cluster_stats <- function(sce,
                              all_cluster_results) {
  pcs <- reducedDim(sce, "PCA")

  # split clustering results by param used
  split_clusters <- all_cluster_results |>
    split(all_cluster_results$nn_param)

  # for each nn_param get cluster width and purity
  all_stats_df <- split_clusters |>
    purrr::map(\(df){
      sil_df <- bluster::approxSilhouette(pcs, df$cluster) |>
        as.data.frame() |>
        tibble::rownames_to_column("barcodes")

      purity_df <- bluster::neighborPurity(pcs, df$cluster) |>
        as.data.frame() |>
        tibble::rownames_to_column("barcodes")

      # join into one data frame to return
      stats_df <- sil_df |>
        dplyr::left_join(purity_df, by = "barcodes")

      return(stats_df)
    }) |>
    dplyr::bind_rows(.id = "nn_param")

  return(all_stats_df)
}

# calculate cluster stability for a single set of clusters using ari
# bootstrap and get ari for clusters compared to sampled clusters
# re-clusters and gets ari across 20 iterations
get_ari <- function(pcs,
                    clusters,
                    k) {
  ari <- c()
  for (iter in 1:20) {
    # sample cells with replacement
    sample_cells <- sample(nrow(pcs), nrow(pcs), replace = TRUE)
    resampled_pca <- pcs[sample_cells, , drop = FALSE]

    # perform clustering on sampled cells
    resampled_clusters <- get_clusters(resampled_pca, k)

    # calculate ARI between new clustering and original clustering
    ari[iter] <- pdfCluster::adj.rand.index(resampled_clusters, clusters[sample_cells])
  }

  ari_df <- data.frame(
    ari = ari,
    k_value = k
  )
}

# get cluster stability for each nn_param cluster results are available for
get_cluster_stability <- function(sce,
                                  all_cluster_results) {
  pcs <- reducedDim(sce, "PCA")

  # split clustering results by param used
  cluster_df_list <- all_cluster_results |>
    split(all_cluster_results$nn_param)

  # for each parameter, get ari values
  cluster_stability_df <- cluster_df_list |>
    purrr::imap(\(df, k_value){
      # make sure k is numeric and remove extra k_
      k <- stringr::str_remove(k_value, "k_") |>
        as.numeric()

      get_ari(pcs, df$cluster, k)
    }) |>
    dplyr::bind_rows()

  return(cluster_stability_df)
}

# Plotting ---------------------------------------------------------------------

# plot individual stats for clusters, either purity or width
plot_cluster_stats <- function(all_stats_df,
                               stat_column) {
  ggplot(all_stats_df, aes(x = nn_param, y = {{ stat_column }})) +
    # ggforce::geom_sina(size = .2) +
    ggbeeswarm::geom_quasirandom(method = "smiley", size = 0.1) +
    stat_summary(
      aes(group = nn_param),
      color = "red",
      # median and quartiles for point range
      fun = "median",
      fun.min = function(x) {
        quantile(x, 0.25)
      },
      fun.max = function(x) {
        quantile(x, 0.75)
      }
    )
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
