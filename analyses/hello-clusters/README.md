# hello-clusters

## Description

This module provides examples of how to use clustering functionality in the [`rOpenScPCA` package](https://github.com/AlexsLemonade/rOpenScPCA/).

When clustering scRNA-seq data, in particular when those clusters are used in downstream analyses, it is important to evaluate the quality of the clusters.
The `rOpenScPCA` package provides several functions that leverage the [`bluster` package](https://bioconductor.org/packages/release/bioc/html/bluster.html) to facilitate performing and evaluating graph-based clustering in a reproducible manner.

### Perform clustering

The function `rOpenScPCA::calculate_clusters()` can be used to perform graph-based clustering.
By default, this function uses the Louvain algorithm with Jaccard weighting.


### Evaluate clustering

`rOpenScPCA` contains several functions to calculate quality metrics for a particular clustering result:

- `rOpenScPCA::calculate_silhouette()`
  - This function calculates the [silhouette width](https://bioconductor.org/books/3.19/OSCA.advanced/clustering-redux.html#silhouette-width)
  - Clusters with higher average silhouette widths indicate that clusters are highly-separated from other clusters
- `rOpenScPCA::calculate_purity()`
  - This function calculates [neighborhood purity](https://bioconductor.org/books/3.19/OSCA.advanced/clustering-redux.html#cluster-purity)
  - Higher purity values indicate that most neighboring cells in a cluster are assigned to the same cluster, therefore indicating good cluster separation
- `rOpenScPCA::calculate_stability()`
  - This function calculates [cluster stability](https://bioconductor.org/books/3.19/OSCA.advanced/clustering-redux.html#cluster-bootstrapping) with Adjusted Rand Index (ARI) and a bootstrapping procedure
  - Higher stability values indicate that clusters are more robust

### Identify optimal clustering parameters

It is often helpful to explore and evaluate results from different clustering algorithms and/or parameters to choose an optimal clustering scheme.
The function `rOpenScPCA::sweep_clusters()` allows you to generate clustering results from a provided set of algorithms and/or parameters, whose quality can then be assessed to select a set of clusters to proceed with.


## Installing rOpenScPCA

The `rOpenScPCA` package is available in the [`AlexsLemonade/rOpenScPCA` repository](https://github.com/AlexsLemonade/rOpenScPCA/).

If you use this package in your analysis module, you should install and track it with `renv` as follows:

```
# First, install rOpenScPCA
renv::install("AlexsLemonade/rOpenScPCA")

# Second, run snapshot to add the package to renv.lock
renv::snapshot()
```

To update the package to its most recent version, you can use the following:

```
renv::update("rOpenScPCA")
```

## Example notebooks

1. The `01_perform-evaluate-clustering.Rmd` notebook shows examples of:
    - Performing clustering with `rOpenScPCA::calculate_clusters()`
    - Evaluating clustering with `rOpenScPCA::calculate_silhouette()`, `rOpenScPCA::calculate_purity()`, and `rOpenScPCA::calculate_stability()`
It also contains explanations for how to interpret cluster quality metrics.

2. The `02_compare-clustering-parameters.Rmd` notebook shows examples of:
    - Performing clustering across a set of parameterizations with `rOpenScPCA::sweep_clusters()`
    - Comparing and visualizing multiple sets of clustering results
