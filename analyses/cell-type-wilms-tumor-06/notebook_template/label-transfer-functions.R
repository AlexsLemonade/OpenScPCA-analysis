# This script contains functions used for label transfer in:
# - "02a_label-transfer_fetal_full_reference_Cao.Rmd"
# - "02b_label-transfer_fetal_full_reference_Stewart.Rmd"
#
# All code was adapted from the RunAzimuth function (https://github.com/satijalab/azimuth/blob/243ee5db80fcbffa3452c944254a325a3da2ef9e/R/azimuth.R), retaining only the specific components of
# this function used for label transfer. For example, additional aspects of this function project
# the query data onto the reference UMAP and perform associated calculations, but this is beyond the
# scope of label transfer for this project, so that code was not imported here.



#' Prepare query Seurat object for label transfer
#'
#' This function prepares Seurat objects by:
#' - Converting gene names to gene symbols and subsetting to only shared features with the reference
#' - Ensuring nCount_RNA and nFeature_RNA are present in the Seurat object
#' - Ensuring the percentage of mitochondrial genes is present in the Seurat object
#'
#' This code was adapted from the RunAzimuth function:
#' https://github.com/satijalab/azimuth/blob/243ee5db80fcbffa3452c944254a325a3da2ef9e/R/azimuth.R
#'
#' @param query The Seurat object which will undergo label transfer
#' @param reference_rownames The rownames (aka, features) in the reference object
#' @param homolog_file Path to the homologs.rds file obtained from Seurat
#'
#' @return Seurat object prepared for label transfer
prepare_query <- function(query, reference_rownames, homolog_file = homologs_file) {
  # Convert the query (sample) row names from ensembl IDs to gene names to match what
  # the Azimuth reference uses
  # Source: https://github.com/satijalab/azimuth/blob/243ee5db80fcbffa3452c944254a325a3da2ef9e/R/azimuth.R#L99-L104
  query <- ConvertGeneNames(
    object = query,
    reference.names = rownames(x = reference),
    homolog.table = "https://seurat.nygenome.org/azimuth/references/homologs.rds"
  )
  query <- Azimuth::ConvertGeneNames(
    object = query,
    reference.names = reference_rownames,
    homolog.table = homolog_file
  )

  # Calculate nCount_RNA and nFeature_RNA if the query does not
  # contain them already
  # Source: https://github.com/satijalab/azimuth/blob/243ee5db80fcbffa3452c944254a325a3da2ef9e/R/azimuth.R#L106-L120
  if (!all(c("nCount_RNA", "nFeature_RNA") %in% c(colnames(x = query[[]])))) {
    calcn <- as.data.frame(x = Seurat:::CalcN(object = query[["RNA"]]))
    colnames(x = calcn) <- paste(
      colnames(x = calcn),
      "RNA",
      sep = "_"
    )
    query <- AddMetaData(
      object = query,
      metadata = calcn
    )
    rm(calcn)
  }

  # Calculate percent mitochondrial genes if the query contains genes
  # matching the regular expression "^MT-"
  # Source: https://github.com/satijalab/azimuth/blob/243ee5db80fcbffa3452c944254a325a3da2ef9e/R/azimuth.R#L122-L131
  if (any(grepl(pattern = "^MT-", x = rownames(x = query)))) {
    query_symbols <- PercentageFeatureSet(
      object = query,
      pattern = "^MT-",
      col.name = "percent.mt"
    )
  }

  return(query)
}

#' Perform label transfer using an Azimuth-adapted approach
#'
#' This function adapts code from the RunAzimuth function to perform label transfer:
#' https://github.com/satijalab/azimuth/blob/243ee5db80fcbffa3452c944254a325a3da2ef9e/R/azimuth.R
#'
#' Note that several params were documented below based on their descriptions here:
#' https://github.com/satijalab/seurat/blob/1549dcb3075eaeac01c925c4b4bb73c73450fc50/R/integration.R
#'
#' @param query The Seurat object which will undergo label transfer
#' @param reference The reference dataset
#' @param reference_dims Dimensions calculated from the reference dataset
#' @param refdata object used by Azimuth
#' @param reference_dims Number of dimensions to use in the anchor weighting procedure
#' @param k.weight Number of neighbors to consider when weighting anchors. This should be
#' <=15 when running on OpenScPCA test data
#' @param n.trees More trees gives higher precision when using annoy approximate
#' nearest neighbor search
#' @param mapping.score.k Compute and store nearest k query neighbors in the
#' AnchorSet object that is returned.
#' @param ksmooth Number of cells to average over when computing transition
#' probabilities
#' @param verbose Display messages/progress
#'
#' @return Seurat object containing metadata from label transfer
transfer_labels <- function(
    query,
    reference,
    reference_dims,
    refdata,
    k.weight = 10,
    n.trees = 20,
    mapping.score.k = 80,
    ksmooth = 80,
    verbose = FALSE) {
  # Find anchors between query and reference
  # Source: https://github.com/satijalab/azimuth/blob/243ee5db80fcbffa3452c944254a325a3da2ef9e/R/azimuth.R#L132-L147
  anchors <- FindTransferAnchors(
    reference = reference,
    query = query,
    k.filter = NA,
    reference.neighbors = "refdr.annoy.neighbors",
    reference.assay = "refAssay",
    query.assay = "RNA",
    reference.reduction = "refDR",
    normalization.method = "SCT",
    features = rownames(Loadings(reference[["refDR"]])),
    dims = 1:reference_dims,
    n.trees = n.trees,
    mapping.score.k = mapping.score.k,
    verbose = verbose
  )

  # Perform label transfer
  # Source: https://github.com/satijalab/azimuth/blob/243ee5db80fcbffa3452c944254a325a3da2ef9e/R/azimuth.R#L164-L175
  query <- TransferData(
    reference = reference,
    query = query,
    query.assay = "RNA",
    dims = 1:reference_dims,
    anchorset = anchors,
    refdata = refdata,
    n.trees = n.trees,
    store.weights = TRUE,
    k.weight = k.weight,
    verbose = verbose
  )


  return(query)
}
