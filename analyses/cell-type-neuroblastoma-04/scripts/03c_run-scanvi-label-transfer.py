#!/usr/bin/env python3

# Script to perform label transfer using scANVI/scArches using the scANVI NBAtlas model on an ScPCA query object prepared with ./03b_prepare-scanvi-query.R
# This scANVI/scArches tutorial was used to help structure this script: https://docs.scvi-tools.org/en/1.3.2/tutorials/notebooks/multimodal/scarches_scvi_tools.html

import argparse
import sys
from pathlib import Path

import anndata
import pandas as pd
from scipy.sparse import csr_matrix
import scvi
import torch

# Define constants
BATCH_KEY = "Sample"
COVARIATE_KEYS = [
    "Assay",
    "Platform",
]
CELL_ID_KEY = "cell_id"
SCANVI_LATENT_KEY = "X_scANVI"
SCANVI_PREDICTIONS_KEY = "scanvi_predicton"


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Annotate SCPCP000004 using scANVI/scArches label transfer with the NBAtlas reference.",
    )
    parser.add_argument(
        "--query_file",
        type=Path,
        required=True,
        help="Path to the input AnnData file which has been prepared with scripts/03a_prepare-scanvi-query.R",
    )
    parser.add_argument(
        "--reference_scanvi_model_dir",
        type=Path,
        required=True,
        help="Path to the load the scANVI/scArches model trained on NBAtlas",
    )
    parser.add_argument(
        "--query_scanvi_model_dir",
        type=Path,
        required=True,
        help="Path to directory where the scANVI/scArches model trained with integrated query data will be saved."
        " This directory will be created at export."
    )
    parser.add_argument(
        "--predictions_tsv",
        type=Path,
        required=True,
        help="Path to the save TSV file of query scANVI/scArches model results including predictions and the scANVI latent representation",
    )
    parser.add_argument(
        "--testing",
        action="store_true",
        default=False,
        help="Flag to use if running on test data and/or in CI",
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
    if not arg.query_file.is_file():
        print(
            f"The provided input query file could not be found at: {arg.query_file}.",
            file=sys.stderr,
        )
        arg_error = True
    if not arg.reference_scanvi_model_dir.is_dir():
        print(
            f"The provided reference scANVI model could not be found at: {arg.reference_scanvi_model_dir}.",
            file=sys.stderr,
        )
        arg_error = True

    # Set up testing and accelerator settings
    if arg.testing:
        # limit max_epochs for faster runtime and ensure CPU
        common_train_kwargs = {"max_epochs": 5, "accelerator": "cpu"}
    else:
        accelerator = "gpu" if torch.cuda.is_available() else "cpu"
        # don't use max_epochs; let scvi pick the heuristic
        common_train_kwargs = {"accelerator": accelerator}

    # Read and check that query object has expected columns
    query = anndata.read_h5ad(arg.query_file)
    expected_columns = [BATCH_KEY, CELL_ID_KEY] + COVARIATE_KEYS

    if not set(expected_columns).issubset(query.obs.columns):
        print(
            f"The query AnnData object is missing one or more expected columns: {set(expected_columns).difference(query.obs.columns)}.",
            file=sys.stderr,
        )
        arg_error = True

    # Exit if error(s)
    if arg_error:
        sys.exit(1)

    # Set seed for reproducibility
    # torch commands are ignored if GPU not present
    scvi.settings.seed = arg.seed  # inherited by numpy and torch
    torch.cuda.manual_seed(arg.seed)
    torch.cuda.manual_seed_all(arg.seed)

    # Load the trained scANVI model for label transfer
    scanvi_model = scvi.model.SCANVI.load(arg.reference_scanvi_model_dir)

    # Ensure anndata is using a sparse matrix for faster processing
    query.X = csr_matrix(query.X)

    ################################################
    # Incorporate query data into the scANVI model #
    ################################################

    # Prepare query data for training
    scvi.model.SCANVI.prepare_query_anndata(query, scanvi_model)
    scanvi_query = scvi.model.SCANVI.load_query_data(query, scanvi_model)

    # Train model and get latent dimensions, cell type predictions
    scanvi_query.train(
        **common_train_kwargs,
        # additional scArches parameters
        plan_kwargs={"weight_decay": 0.0},
        check_val_every_n_epoch=5,
    )
    query.obsm[SCANVI_LATENT_KEY] = scanvi_query.get_latent_representation()
    query.obs[SCANVI_PREDICTIONS_KEY] = scanvi_query.predict()

    ################################################
    ################ Export objects ################
    ################################################

    # Export the query-trained scANVI model
    scanvi_query.save(arg.query_scanvi_model_dir, save_anndata=True, overwrite=True)

    # Extract data frame of latent representation and predictions
    latent_df = pd.DataFrame(query.obsm[SCANVI_LATENT_KEY])
    latent_df = latent_df.rename(columns=lambda x: f"{SCANVI_LATENT_KEY}_{x}")
    latent_df.index = query.obs.index # set index for joining
    combined_df = latent_df.join(
        query.obs[ expected_columns + [SCANVI_PREDICTIONS_KEY] ]
    )
    combined_df.to_csv(arg.predictions_tsv, sep="\t", index=False)

if __name__ == "__main__":
    main()
