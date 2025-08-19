#!/usr/bin/env python3

# Script to annotate processed ScPCA objects with SCimilarity 
# Follows this tutorial: https://genentech.github.io/scimilarity/notebooks/cell_annotation_tutorial.html

import argparse
import sys
from pathlib import Path

import anndata
import pandas as pd
from scimilarity import CellAnnotation
from scimilarity.utils import align_dataset, lognorm_counts

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Annotate ScPCA samples using SCimilarity",
    )
    parser.add_argument(
        "--model_dir",
        type=Path,
        required=True,
        help="Path to the directory containing the SCimilarity foundational model",
    )
    parser.add_argument(
        "--processed_h5ad_file",
        type=Path,
        required=True,
        help="Path to the processed AnnData object stored as an h5ad file",
    )
    parser.add_argument(
        "--predictions_tsv",
        type=Path,
        required=True,
        help="Path to the output TSV file with cell type annotations",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=2025,
        help="Random seed to ensure reproducibility",
    )
    arg = parser.parse_args()

    ################################################
    ########### Input argument checks ##############
    ################################################
    arg_error = False

    # Check that input files exist
    if not arg.model_dir.is_dir():
        print(
            f"The provided reference SCimilarity model could not be found at: {arg.model_dir}.",
            file=sys.stderr,
        )
        arg_error = True
    if not arg.processed_h5ad_file.is_file():
        print(
            f"The provided input H5AD file could not be found at: {arg.processed_h5ad_file}.",
            file=sys.stderr,
        )
        arg_error = True

     # Exit if error(s)
    if arg_error:
        sys.exit(1)

    # Set seed for reproducibility
    scvi.settings.seed = arg.seed  # inherited by numpy and torch

    ################################################
    ################ Prep input data ###############
    ################################################

    # Read in model 
    scimlarity_model = CellAnnotation(model_path = arg.model_dir)

    # Read and check that query object has expected columns
    processed_anndata = anndata.read_h5ad(arg.processed_h5ad_file)

    # Ensure anndata is using a sparse matrix for faster processing
    processed_anndata.X = csr_matrix(processed_anndata.X)

    # Preprocess the data
    # Align the query dataset to the reference model
    processed_anndata = align_dataset(processed_anndata, scimlarity_model.gene_order)
    # Log-normalize the counts
    processed_anndata = lognorm_counts(processed_anndata)

    # compute embeddings
    processed_anndata.obsm["X_scimilarity"] = scimlarity_model.get_embeddings(processed_anndata.X)

    ################################################
    ############### Run Scimilarity ###############
    ################################################

    # Predict cell types
    predictions, nn_idxs, nn_dists, nn_stats = scimlarity_model.get_predictions_knn(
        processed_anndata.obsm["X_scimilarity"], weighting=True
    ) 

    ################################################
    ################ Export annotations ############
    ################################################

    # prepare the predictions with posterior probabilities for export
    predictions_df = predictions.values
    min_dist_df = nn_dists.min(axis=1)
    predictions_df = predictions_df.join(min_dist_df)

    # export TSV
    predictions_df.to_csv(arg.predictions_tsv, sep="\t", index=False)


if __name__ == "__main__":
    main()
