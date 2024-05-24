#!/usr/bin/env python3

# This script runs scrublet on a given AnnData file and exports results as a TSV.


import argparse
import anndata
import scrublet
import pandas
from pathlib import Path
import sys

# Function to run scrublet
def run_scrublet(adata: anndata.AnnData, random_seed: int) -> pandas.DataFrame:
    """
        Run scrublet on a counts matrix from an AnnData object, following https://github.com/swolock/scrublet?tab=readme-ov-file#quick-start
        Returns a dataframe with the following columns:
            - `barcodes` is the droplet's unique barcode
            - `scrublet_score` is the probability that each droplet is a doublet
            - `scrublet_prediction` is the predicted droplet status ("doublet" or "singlet"), based on an automatically-set threshold
    """

    scrub = scrublet.Scrublet(adata.X, random_state = random_seed)
    # predicted_doublets is boolean array, where True are predicted doublets and False are predicted singlets
    doublet_scores, predicted_doublets = scrub.scrub_doublets()

    # convert True/False to string values
    predicted_doublets_str = [("doublet" if x else "singlet") for x in predicted_doublets]

    results = {"barcodes" : adata.obs_names, "scrublet_score" : doublet_scores, "scrublet_prediction" : predicted_doublets_str}
    results_df = pandas.DataFrame(results)

    return(results_df)

def main() -> None:

    parser = argparse.ArgumentParser(
        description="Detect doublets on a set of AnnData objects using scrublet.",
    )
    parser.add_argument(
        "--dataset_name",
        type=str,
        default="",
        help="Name of dataset to process, where the associated file is expected to be named `{name}_anndata.h5ad`."
    )
    parser.add_argument(
        "--data_dir",
        type=Path,
        help="The directory containing H5AD input file."
    )
    parser.add_argument(
        "--results_dir",
        type=Path,
        help="The directory to export TSV file with doublet inferences."
    )
    parser.add_argument(
        "--random_seed",
        type=int,
        default=2024,
        help="Random seed to ensure reproducibility."
    )
    args = parser.parse_args()

    # Prepare input arguments
    if not args.dataset_name:
        print(
            "Datasets must be provided with the `--dataset_name` flag.",
            file=sys.stderr
        )
        sys.exit(1)
    if not args.data_dir.exists():
        print(
            "A correct path to the input data must be provided with --data_dir.",
            file=sys.stderr
        )
        sys.exit(1)
    args.results_dir.mkdir(parents = True, exist_ok = True)

    # Run scrublet and export the results
    input_anndata = args.data_dir / args.dataset_name + "_anndata.h5ad"
    result_tsv = args.dataset_name + "_scrublet.tsv"

    adata = anndata.read_h5ad( args.data_dir / input_anndata )
    scrub_results = run_scrublet(adata, args.random_seed)
    scrub_results.to_csv( args.results_dir / result_tsv, sep="\t", index=False )

if __name__ == "__main__":
    main()
