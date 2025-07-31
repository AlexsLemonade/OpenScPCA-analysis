#!/usr/bin/env python3

# Script to train a scANVI/scArches model from the NBAtlas reference
# This scANVI/scArches tutorial was used to help structure this script: https://docs.scvi-tools.org/en/1.3.2/tutorials/notebooks/multimodal/scarches_scvi_tools.html

import os
import sys
from pathlib import Path
import argparse
import anndata
import scvi
import torch
from scipy.sparse import csr_matrix

# Define variables to use in the objects
BATCH_KEY = "Sample"
COVARIATE_KEYS = [
    "Assay",
    "Platform",
]  # TODO: make this a script option to compare with/without?
SCVI_LATENT_KEY = "X_scVI"
SCANVI_LATENT_KEY = "X_scANVI"
SCANVI_LABELS_KEY = "labels_scanvi"
UNLABELED_VALUE = "Unknown"
SCANVI_PREDICTIONS_KEY = "predictions_scanvi"


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
        "--reference_scvi_model_dir",
        type=Path,
        required=True,
        help="Path to the save the intermediate SCVI model trained on the reference object",
    )
    parser.add_argument(
        "--reference_scanvi_model_dir",
        type=Path,
        required=True,
        help="Path to the save the scANVI model trained on the reference object",
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
    if not os.path.isfile(arg.reference_file):
        print(
            f"The provided input reference file could not be found at: {arg.reference_file}.",
            file=sys.stderr,
        )
        arg_error = True

    # Set up testing and accelerator settings
    # Note we do this here so `cell_type_label` can be used to check the reference object next
    if arg.testing:
        cell_type_label = "Cell_type"
        # limit max_epochs for faster runtime and ensure CPU
        common_train_kwargs = {"max_epochs": 5, "accelerator": "cpu"}
    else:
        cell_type_label = arg.reference_celltype_column
        accelerator = "gpu" if torch.cuda.is_available() else "cpu"
        # don't use max_epochs; let scvi pick the heuristic
        common_train_kwargs = {"accelerator": accelerator}

    # Define lists of expected columns in the reference and query objects
    expected_columns = [BATCH_KEY, cell_type_label] + COVARIATE_KEYS

    # Read the reference and query objects and check that expected columns are present
    reference = anndata.read_h5ad(arg.reference_file)
    all_present_reference = all(
        col in reference.obs.columns for col in expected_columns
    )
    if not all_present_reference:
        print(
            f"The reference AnnData object is missing one or more expected columns: {expected_columns}.",
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

    ################################################
    ######## SCVI reference model training #########
    ################################################

    # Ensure anndata objects are using sparse matries for faster scvi processing
    reference.X = csr_matrix(reference.X)

    scvi.model.SCVI.setup_anndata(
        reference,
        batch_key=BATCH_KEY,
        categorical_covariate_keys=COVARIATE_KEYS,  # control for cell/nucleus and 10x2/3
    )

    scvi_model = scvi.model.SCVI(
        reference,
        # scArches parameters
        use_layer_norm="both",
        use_batch_norm="none",
        encode_covariates=True,  # essential for scArches
        dropout_rate=0.2,
        n_layers=2,
    )

    scvi_model.train(**common_train_kwargs)
    reference.obsm[SCVI_LATENT_KEY] = scvi_model.get_latent_representation()

    ################################################
    ####### scANVI reference model training ########
    ################################################

    reference.obs[SCANVI_LABELS_KEY] = reference.obs[cell_type_label].values

    scanvi_model = scvi.model.SCANVI.from_scvi_model(
        scvi_model,
        unlabeled_category=UNLABELED_VALUE,
        labels_key=SCANVI_LABELS_KEY,
    )
    scanvi_model.train(**common_train_kwargs)
    reference.obsm[SCANVI_LATENT_KEY] = scanvi_model.get_latent_representation()

    ################################################
    ################ Export objects ################
    ################################################

    # Export the NBAtlas-trained SCVI model
    scvi_model.save(arg.reference_scvi_model_dir, save_anndata=True, overwrite=True)

    # Export the NBAtlas-trained scANVI model
    scanvi_model.save(arg.reference_scanvi_model_dir, save_anndata=True, overwrite=True)

if __name__ == "__main__":
    main()
