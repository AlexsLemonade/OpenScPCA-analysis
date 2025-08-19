#!/usr/bin/env python3

# Script to annotate processed ScPCA objects with SCimilarity 
# Follows this tutorial: https://genentech.github.io/scimilarity/notebooks/cell_annotation_tutorial.html

import argparse
import sys
from pathlib import Path

import anndata
import pandas
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
        default="test.tsv",
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
    
    # Set the seed
    random.seed(arg.seed)

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

    ################################################
    ################ Prep input data ###############
    ################################################

    # Read in model 
    scimilarity_model = CellAnnotation(model_path = model_dir)

    # Read and make sure object formatting is correct
    processed_anndata = anndata.read_h5ad(processed_h5ad_file)
    
    # counts should be stored as a layer and not in X
    processed_anndata.layers['counts'] = processed_anndata.X
    # rownames should correspond to gene symbols, not IDs
    processed_anndata.var_names = processed_anndata.var["gene_symbol"].astype(str)
    processed_anndata.var_names_make_unique()

    # Preprocess the data
    # Align the query dataset to the reference model
    processed_anndata = align_dataset(processed_anndata, scimlarity_model.gene_order)
    # Log-normalize the counts
    processed_anndata = lognorm_counts(processed_anndata)

    ################################################
    ############### Run Scimilarity ###############
    ################################################
    
    # compute embeddings
    processed_anndata.obsm["X_scimilarity"] = scimlarity_model.get_embeddings(processed_anndata.X)

    # Predict cell types
    predictions, nn_idxs, nn_dists, nn_stats = scimlarity_model.get_predictions_knn(
        processed_anndata.obsm["X_scimilarity"], weighting=True
    ) 

    ################################################
    ################ Export annotations ############
    ################################################

    # prepare the predictions with min distance for export 
    predictions_df = pandas.DataFrame({
      "barcode": processed_anndata.obs_names.to_list(), 
      "scimilarity_celltype_annotation": predictions.values, 
      "min_dist": nn_stats["min_dist"]
    })

    # export TSV
    predictions_df.to_csv(arg.predictions_tsv, sep="\t", index=False)


if __name__ == "__main__":
    main()
