#!/usr/bin/env python3

# Script to perform cell type annotation on SCPCP000004 using scANVI/scArches label transfer with the NBAtlas reference
# This scANVI/scArches tutorial was used to help structure this script: https://docs.scvi-tools.org/en/1.3.2/tutorials/notebooks/multimodal/scarches_scvi_tools.html

import argparse
import os
from pathlib import Path
import re
import subprocess
import sys
import anndata
import scvi
import torch
from scipy.sparse import csr_matrix


# Define variables to use in the objects
BATCH_KEY = "Sample"
COVARIATE_KEYS = ["Assay", "Platform"]
SCANVI_LABELS_KEY = "labels_scanvi"
UNLABELED_VALUE = "Unknown"
SCANVI_PREDICTIONS_KEY = "predictions_scanvi"
SCANVI_LATENT_KEY = "X_scANVI"


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Annotate SCPCP000004 using scANVI/scArches label transfer with the NBAtlas reference.",
    )
    parser.add_argument(
        "--reference_file",
        type=Path,
        required=True,
        help="Path to the input NBAtlas reference file",
    )
    parser.add_argument(
        "--query_file",
        type=Path,
        required=True,
        help="Path to the input prepared ScPCA merged AnnData file",
    )
    parser.add_argument(
        "--reference_celltype_column",
        type=str,
        required=True,
        default = "Cell_type_wImmuneZoomAnnot",
        help="Column in the reference AnnData object that contains cell type annotations. Default is 'Cell_type_wImmuneZoomAnnot'.",
    )
    parser.add_argument(
        "--reference_scvi_model_file",
        type=Path,
        required=True,
        help="Path to the save the SCVI model trained on the reference object",
    )
    parser.add_argument(
        "--reference_scanvi_model_file",
        type=Path,
        required=True,
        help="Path to the save the scANVI model trained on the reference object",
    )
    parser.add_argument(
        "--full_scanvi_model_file",
        type=Path,
        required=True,
        help="Path to the save the full scANVI/scArches model trained with integrated query data",
    )
    parser.add_argument(
        "--exclude_libraries",
        type=str,
        required=False,
        help="Optionally, a comma-separated list of ScPCA libraries to exclude from analysis",
    )
    parser.add_argument(
        "--output_tsv_file",
        type=Path,
        required=True,
        help="Path to the save a TSV with scANVI/scArches cell type annotations",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=2025,
        help="Random seed to ensure reproducibility",
    )
    parser.add_argument(
        "--accelerator",
        type=str,
        default="cpu",
        help="Use 'gpu' for GPU acceleration or 'cpu' for CPU only. Default is 'cpu'.",
    )
    args = parser.parse_args()

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
    if not os.path.isfile(arg.query_file):
        print(
            f"The provided input query file could not be found at: {arg.query_file}.",
            file=sys.stderr,
        )
        arg_error = True

    # Check that the accelerator is valid and that GPU is available is specified
    if arg.accelerator not in ["gpu", "cpu"]:
        print(
            f"The provided accelerator '{arg.accelerator}' is not valid. Use 'gpu' for GPU acceleration or 'cpu' for CPU only.",
            file=sys.stderr,
        )
        arg_error = True
    if arg.accelerator == "gpu" and not torch.cuda.is_available():
        print(
            "The specified accelerator is 'gpu', but no GPU is available. Please use 'cpu' instead or ensure a GPU is available.",
            file=sys.stderr,
        )
        arg_error = True


    # Define lists of expected columns in the reference and query objects
    expected_covariate_columns = [BATCH_KEY] + COVARIATE_KEYS
    reference_expected_columns = expected_covariate_columns + [arg.reference_celltype_column]

    # Read the reference and query objects and check that expected columns are present
    reference = anndata.read_h5ad(arg.reference_file)

    all_present_reference = all(col in reference.obs.columns for col in reference_expected_columns)
    if not all_present_reference:
        print(
            f"The reference AnnData object is missing one or more expected columns: {reference_expected_columns}.",
            file=sys.stderr,
        )
        arg_error = True

    query = anndata.read_h5ad(arg.query_file)
    all_present_query = all(col in query.obs.columns for col in expected_covariate_columns)
    if not all_present_query:
        print(
            f"The query AnnData object is missing one or more expected columns: {expected_covariate_columns}.",
            file=sys.stderr,
        )
        arg_error = True

    # Prepare and check library ids to exclude
    exclude_libraries = [library_id.strip() for library_id in arg.exclude_libraries.split(',')]
    if len(exclude_libraries) > 0:
        if not all(exclude_libraries in query.obs["library_id"]):
            print(
                f"Not all libraries in `exclude_libraries` are present in the query object. Please check ids: {arg.exclude_libraries}.",
                file=sys.stderr,
            )
            arg_error = True

        if arg_error:
            sys.exit(1)

    # Set seed for reproducibility
    # torch commands are ignored if GPU not present
    scvi.settings.seed = arg.seed # inherited by numpy and torch
    torch.cuda.manual_seed(arg.seed)
    torch.cuda.manual_seed_all(arg.seed)

    ################################################
    ######## SCVI reference model training #########
    ################################################

    # Ensure anndata objects are using sparse matries for faster scvi processing
    reference.X = csr_matrix(reference.X)
    query.X = csr_matrix(query.X)

    # remove libraries from the query if specified
    if len(exclude_libraries) > 0:
        query_mask = ~query.obs["library_id"].isin(exclude_libraries)
        query = query[query_mask, :]


    scvi.model.SCVI.setup_anndata(
        reference,
        batch_key = BATCH_KEY,
        categorical_covariate_keys = COVARIATE_KEYS # control for cell/nucleus and 10x2/3
    )

    scvi_model = scvi.model.SCVI(
        reference,
        # scArches parameters
        use_layer_norm="both",
        use_batch_norm="none",
        encode_covariates=True,
        dropout_rate=0.2,
        n_layers=2,
    )

    scvi_model.train(
        accelerator = arg.accelerator,
        early_stopping = True,
        early_stopping_patience = 3, # stop after 3 epochs without improvement; the max_epochs for this data is 22
    )

    # Export the trained SCVI model
    scvi_model.save(arg.reference_scvi_model_file, save_anndata = True, overwrite=True)


    ################################################
    ####### scANVI reference model training ########
    ################################################

    reference.obs[SCANVI_LABELS_KEY] = reference.obs[arg.reference_celltype_column].values

    scanvi_model = scvi.model.SCANVI.from_scvi_model(
        scvi_model,
        unlabeled_category = UNLABELED_VALUE,
        labels_key = SCANVI_LABELS_KEY,
    )
    scanvi_model.train(
        accelerator = arg.accelerator,
        early_stopping = True,
        early_stopping_patience = 3
    )
    # Export the trained scANVI model
    scanvi_model.save(arg.reference_scanvi_model_file, save_anndata = True, overwrite=True)


    ################################################
    # Incorporate query data into the scANVI model #
    ################################################

    # Incorporate query data into the scANVI model
    scvi.model.SCANVI.prepare_query_anndata(query, scanvi_model)
    scanvi_query = scvi.model.SCANVI.load_query_data(query, scanvi_model)

    # https://docs.scvi-tools.org/en/1.3.2/tutorials/notebooks/multimodal/scarches_scvi_tools.html#id2
    scanvi_query.train(
        accelerator = arg.accelerator,
        # scArches parameters
        plan_kwargs={"weight_decay": 0.0},
        check_val_every_n_epoch=10,
    )

    query.obsm[SCANVI_LATENT_KEY] = scanvi_query.get_latent_representation()
    query.obs[SCANVI_PREDICTIONS_KEY] = scanvi_query.predict()

    # Save the full scANVI model with integrated query data
    scanvi_query.save(arg.full_scanvi_model_file, save_anndata = True, overwrite=True)

    # Save TSV of the predictions
    scanvi_df = query.obs[['cell_id', SCANVI_PREDICTIONS_KEY]]
    scanvi_df.to_csv(arg.output_tsv_file, sep='\t', index=True)


if __name__ == "__main__":
    main()
