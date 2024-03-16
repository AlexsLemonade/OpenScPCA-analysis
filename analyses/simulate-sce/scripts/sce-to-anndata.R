#!/usr/bin/env Rscript

# This script takes a SingleCellExperiment stored in a .rds file and converts the main experiment
# (usually RNA) to an AnnData object saved as an hdf5 file

# The AnnData object being exported by this script is formatted to fit CZI schema: 3.0.0

# adapted from https://github.com/AlexsLemonade/scpca-nf/blob/8bef82d853d19e5aeddd75401aa54cf8bfbced13/bin/sce_to_anndata.R


library(optparse)

# set up arguments
option_list <- list(
  make_option(
    opt_str = c("-d", "--dir"),
    type = "character",
    help = "Path where rds files are located. Will search recursively for files ending in `_processed.rds`, `_filtered.rds`, or `_unfiltered.rds.`",
  )
)

opts <- parse_args(OptionParser(option_list = option_list))

# Load libraries ---------------------------------------------------------------

# wait to load these until after the command line arguments are parsed
# Some take a while to load, and we'd rather have the script be responsive

# There are some annoying warning here, let's suppress them
suppressPackageStartupMessages({
  suppressWarnings(library(SingleCellExperiment))
})


# Functions  -------------------------------------------------------------------

# this function updates merged object formatting for anndata export
outfile_name <- function(sce_filename, feature) {
  # paste X to any present reduced dim names
  suffix <- glue::glue("_{feature}.hdf5")
  feature_filename <- stringr::str_replace(sce_filename, ".rds$", suffix)
  return(feature_filename)
}

# Merged object function  ------------------------------------------------------

# this function updates merged object formatting for anndata export
format_merged_sce <- function(sce) {
  # paste X to any present reduced dim names
  reducedDimNames(sce) <- glue::glue("X_{reducedDimNames(sce)}")
  return(sce)
}



# CZI compliance function ------------------------------------------------------

# this function applies any necessary reformatting or changes needed to make
# sure that the sce that is getting converted to AnnData is compliant with
# CZI 3.0.0 requirements: https://github.com/chanzuckerberg/single-cell-curation/blob/b641130fe53b8163e50c39af09ee3fcaa14c5ea7/schema/3.0.0/schema.md
format_czi <- function(sce) {
  # add schema version
  metadata(sce)$schema_version <- "3.0.0"

  # add library_id as an sce colData column
  # need this column to join in the sample metadata with the colData
  if (!("library_id" %in% colnames(colData(sce)))) {
    sce$library_id <- metadata(sce)$library_id
  }

  # only move sample metadata if not a multiplexed library
  if (!("cellhash" %in% altExpNames(sce))) {
    # add sample metadata to colData sce
    sce <- scpcaTools::metadata_to_coldata(
      sce,
      join_columns = "library_id"
    )
  }

  # modify colData to be AnnData and CZI compliant
  coldata_df <- colData(sce) |>
    as.data.frame() |>
    dplyr::mutate(
      # create columns for assay and suspension ontology terms
      assay_ontology_term_id = metadata(sce)$assay_ontology_term_id,
      suspension_type = metadata(sce)$seq_unit,
      # add is_primary_data column; only needed for anndata objects
      is_primary_data = FALSE
    )

  # add colData back to sce object
  colData(sce) <- DataFrame(
    coldata_df,
    row.names = rownames(colData(sce))
  )

  # remove sample metadata from sce metadata, otherwise conflicts with converting object
  metadata(sce) <- metadata(sce)[names(metadata(sce)) != "sample_metadata"]

  # modify rowData
  # we don't do any gene filtering between normalized and raw counts matrix
  # so everything gets set to false
  rowData(sce)$feature_is_filtered <- FALSE

  # paste X to any present reduced dim names
  reducedDimNames(sce) <- glue::glue("X_{reducedDimNames(sce)}")

  return(sce)
}

# SCE conversion function ------------------------------------------------------
convert_sce_file <- function(sce_file) {
  # read in sce
  sce <- readr::read_rds(sce_file)

  # grab sample metadata
  # we need this if we have any feature data that we need to add it o
  sample_metadata <- metadata(sce)$sample_metadata

  # make main sce czi compliant for single objects, or format merged objects
  sce <- format_czi(sce)

  # export sce as anndata object
  # this function will also remove any R-specific object types from the SCE metadata
  #   before converting to AnnData
  scpcaTools::sce_to_anndata(
    sce,
    anndata_file = outfile_name(sce_file, "rna"),
    compression = "gzip"
  ) |> suppressMessages() # suppress notes about metadata conversion

  # convert altExps to AnnData
  altExpNames(sce) |>
    purrr::walk(\(feature_name){
      convert_altexp(
        sce,
        feature_name,
        outfile_name(sce_file, feature_name),
        opt$compress_output
      )
    })
}

# AltExp to AnnData -----------------------------------------------------------

# if feature data exists, grab it and export to AnnData
convert_altexp <- function(sce, feature_name, output_feature_h5, compress_output) {
  # make sure the feature data is present
  if (feature_name == "cellhash") {
    warning("Conversion of altExp data from multiplexed data is not supported.
             The altExp will not be converted.")
    return()
  }
  # create output file name


  # extract altExp
  alt_sce <- altExp(sce, feature_name)

  # only convert altExp with > 1 rows
  if (nrow(alt_sce) > 1) {
    # add sample metadata from main sce to alt sce metadata
    metadata(alt_sce)$sample_metadata <- sample_metadata

    # make sce czi compliant
    alt_sce <- format_czi(alt_sce)

    # export altExp sce as anndata object
    scpcaTools::sce_to_anndata(
      alt_sce,
      anndata_file = output_feature_h5
    )
  } else {
    # warn that the altExp cannot be converted
    warning(
      glue::glue("
        Only 1 row found in altExp named: {feature_name}.
        This altExp will not be converted to an AnnData object.
      ")
    )
  }
}

# Main script body -------------------------------------------------------------

# get file list
sce_files <- list.files(
  opts$dir,
  pattern = "_(processed|filtered|unfiltered).rds$",
  recursive = TRUE,
  full.names = TRUE
)

# convert each file
purrr::walk(sce_files, convert_sce_file)
