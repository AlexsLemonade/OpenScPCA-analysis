#!/usr/bin/env python3


import argparse
import anndata
import scrublet
import pandas
from pathlib import Path
import pyprojroot
import sys

# Function to run scrublet
def run_scrublet(adata: anndata.AnnData) -> pandas.DataFrame:
    """
        Run scrublet on a counts matrix from an AnnData object, following https://github.com/swolock/scrublet?tab=readme-ov-file#quick-start
        Returns a dataframe with the following columns:
            - `barcodes` is the droplet's unique barcode
            - `scrublet_score` is the probability that each droplet is a doublet
            - `scrublet_prediction` is the predicted droplet status ("doublet" or "singlet"), based on an automatically-set threshold
    """

    scrub = scrublet.Scrublet(adata.X)
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
        "--datasets",
        type=str,
        default="",
        help="Names of datasets to process as a comma-separated list."
               " Datasets are expected to be named `{name}_anndata.h5ad`."
    )
    parser.add_argument(
        "--data_dir",
        type=Path,
        help="The directory containing input H5AD files."
    )
    parser.add_argument(
        "--results_dir",
        type=Path,
        help="The directory to export TSV files with doublet inferences."
    )

    args = parser.parse_args()

    # Prepare input arguments
    dataset_names = [p.strip() for p in args.datasets.split(",")] if args.datasets else []
    if len(dataset_names) == 0:
        print(
            "Datasets must be provided with the `--datasets` flag.",
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

    # Run scrublet on each dataset and export the results
    for dataname in dataset_names:
        input_anndata = dataname + "_anndata.h5ad"
        result_tsv = dataname + "_scrublet.tsv"

        adata = anndata.read_h5ad( args.data_dir / input_anndata )
        scrub_results = run_scrublet(adata)
        scrub_results.to_csv( args.results_dir / result_tsv, sep="\t", index=False )

if __name__ == "__main__":
    main()
