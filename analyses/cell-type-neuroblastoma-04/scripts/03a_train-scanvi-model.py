#!/usr/bin/env python3

# Script to train a scANVI/scArches model from the NBAtlas reference
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


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Train a scANVI/scArches model from the NBAtlas reference.",
    )
    parser.add_argument(
        "--reference_file",
        type=Path,
        required=True,
        help="Path to the input NBAtlas reference file",
    )
    parser.add_argument(
        "--reference_celltype_column",
        type=str,
        default="Cell_type_wImmuneZoomAnnot",
        help="Column in the reference AnnData object that contains cell type annotations."
        " Default is 'Cell_type_wImmuneZoomAnnot', unless --testing is specified in which case `Cell_type` is the default",
    )
    parser.add_argument(
        "--reference_scanvi_model_dir",
        type=Path,
        required=True,
        help="Path to directory where the scANVI model trained on the reference object will be saved."
        " This directory will be created at export."
    )
    parser.add_argument(
        "--scanvi_latent_tsv",
        type=Path,
        help="Optionally, path to save a standalone TSV of the scANVI model latent representation on the reference object."
        " If not provided, no TSV is exported.",
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

    # Check that the input file exists
    if not arg.reference_file.is_file():
        print(
            f"The provided input reference file could not be found at: {arg.reference_file}.",
            file=sys.stderr,
        )
        arg_error = True

    # Set up testing and accelerator settings
    if arg.testing:
        # limit max_epochs for faster runtime and ensure CPU
        common_train_kwargs = {"max_epochs": 5, "accelerator": "cpu"}
        cell_type_column = "Cell_type"
    else:
        # don't use max_epochs; let scvi pick the heuristic
        common_train_kwargs = {
            "accelerator": "gpu" if torch.cuda.is_available() else "cpu"
        }
        cell_type_column = arg.reference_celltype_column

    # Read and check that the reference object has expected columns
    reference = anndata.read_h5ad(arg.reference_file)
    expected_columns = [BATCH_KEY, CELL_ID_KEY, cell_type_column] + COVARIATE_KEYS

    if not set(expected_columns).issubset(reference.obs.columns):
        print(
            f"The reference AnnData object is missing one or more expected columns: {set(expected_columns).difference(reference.obs.columns)}.",
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

    # Ensure anndata is using a sparse matrix for faster processing
    reference.X = csr_matrix(reference.X)

    ################################################
    ######## SCVI reference model training #########
    ################################################

    scvi.model.SCVI.setup_anndata(
        reference,
        batch_key=BATCH_KEY,
        categorical_covariate_keys=COVARIATE_KEYS,  # control for cell/nucleus and 10x2/3
    )

    scvi_model = scvi.model.SCVI(
        reference,
        # scArches parameters
        # from: https://docs.scvi-tools.org/en/1.3.2/tutorials/notebooks/multimodal/scarches_scvi_tools.html#train-reference
        use_layer_norm="both",
        use_batch_norm="none",
        encode_covariates=True,  # essential for scArches
        dropout_rate=0.2,
        n_layers=2,
    )

    # Train SCVI model
    scvi_model.train(**common_train_kwargs)

    ################################################
    ####### scANVI reference model training ########
    ################################################

    scanvi_model = scvi.model.SCANVI.from_scvi_model(
        scvi_model,
        unlabeled_category="Unknown", # will be used to set up query next; labels will start as `Unknown`
        labels_key=cell_type_column,
    )
    scanvi_model.train(**common_train_kwargs)
    reference.obsm[SCANVI_LATENT_KEY] = scanvi_model.get_latent_representation()

    ################################################
    ################ Export objects ################
    ################################################

    # Export the NBAtlas-trained scANVI model
    scanvi_model.save(arg.reference_scanvi_model_dir, overwrite=True, save_anndata=True)

    # Export standalone TSV of NBAtlas scANVI latent representation & labels if specified
    if arg.scanvi_latent_tsv:
        latent_df = pd.DataFrame(reference.obsm[SCANVI_LATENT_KEY])
        latent_df = latent_df.rename(columns=lambda x: SCANVI_LATENT_KEY + "_" + str(x))
        latent_df.index = reference.obs.index # set index for joining
        combined_df = latent_df.join(reference.obs[expected_columns]) # only save relevant columns here
        combined_df.to_csv(arg.scanvi_latent_tsv, sep="\t", index=False)

if __name__ == "__main__":
    main()
