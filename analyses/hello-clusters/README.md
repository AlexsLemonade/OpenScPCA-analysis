# hello-clusters

## Description

This module provides examples of how to use clustering functionality in the `rOpenScPCA` package. 

When clustering scRNA-seq data, in particular when those clusters are used in downstream analyses, it is important to evaluate the quality of the clusters.
The `rOpenScPCA` package provides several functions that leverage the [`bluster` package](https://bioconductor.org/packages/release/bioc/html/bluster.html) to facilitate performing and evaluating graph-based clustering in a reproducible manner.

### Perform clustering

The function `calculate_clusters()` can be used to perform graph-based clustering.
By default, this function uses the Louvain algorithm with Jaccard weighting.
### Evaluate clustering

`rOpenScPCA` contains several functions to calculate quality metrics for a particular clustering result:

- `calculate_silhouette()`
  - This function calculates the [silhouette width](https://bioconductor.org/books/3.19/OSCA.advanced/clustering-redux.html#silhouette-width)
  - Clusters with higher average silhouette widths indicate that clusters are highly-separated from other clusters
- `calculate_purity()`
  - This function calculates [neighborhood purity](https://bioconductor.org/books/3.19/OSCA.advanced/clustering-redux.html#cluster-purity)
  - Higher purity values indicate that most neighboring cells in a cluster are assigned to the same cluster, therefore indicating good cluster separation
- `calculate_stability()`
  - This function calculates [cluster stability](https://bioconductor.org/books/3.19/OSCA.advanced/clustering-redux.html#cluster-bootstrapping) with Adjusted Rand Index (ARI) and a bootstrapping procedure
  - Higher stability values indicate that clusters are more robust
  
### Identify optimal clustering parameters

It is often helpful to explore and evaluate results from different clustering algorithms and/or parameters to choose an optimal clustering scheme.
The function `sweep_clusters()` allows you to generate clustering results from a provided set of algorithms and/or parameters, whose quality can then be assessed to select a set of clusters to proceed with.
  


## Installing rOpenScPCA

The `rOpenScPCA` package is disseminated in the `OpenScPCA-analysis` repository in the `packages` directory. 
If you use this package in your analysis module, you should install and track it with `renv` as follows:

```
# First, install rOpenScPCA
renv::install("AlexsLemonade/OpenScPCA-analysis:packages/rOpenScPCA")



# Second, run snapshot to add the package to renv.lock
renv::snapshot()
```

## Example notebooks

_Content forthcoming._
