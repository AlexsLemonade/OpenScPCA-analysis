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

    # Predict doublets if possible. If it fails, return all NAs
    try:
        # predicted_doublets is boolean array, where True are predicted doublets and False are predicted singlets
        doublet_scores, predicted_doublets = scrub.scrub_doublets()
        # convert True/False to string values
        predicted_doublets_str = [("doublet" if x else "singlet") for x in predicted_doublets]
    except:
        print(
            "Could not run scrublet on the provided dataset."
            "The output TSV file will contain only NA values."
        )
        doublet_scores = "NA"
        predicted_doublets_str = "NA"

    results = {"barcodes" : adata.obs_names, "scrublet_score" : doublet_scores, "scrublet_prediction" : predicted_doublets_str}
    results_df = pandas.DataFrame(results)

    return(results_df)

def main() -> None:

    parser = argparse.ArgumentParser(
        description="Detect doublets on a set of AnnData objects using scrublet.",
    )
    parser.add_argument(
        "--input_anndata_file",
        type=Path,
        required=True,
        help="Path to the input AnnData file to process."
    )
    parser.add_argument(
        "--output_file",
        type=Path,
        required=True,
        help="Path to output TSV file with doublet inferences."
    )
    parser.add_argument(
        "--random_seed",
        type=int,
        default=2024,
        help="Random seed to ensure reproducibility."
    )
    args = parser.parse_args()

    # Define and check files, directories
    arg_error=False
    if not args.input_anndata_file.exists():
        print(
            f"The input AnnData file could not be found at: {args.input_anndata_file}.",
            file=sys.stderr
        )
        arg_error=True

    if not args.output_file.name.endswith(".tsv"):
        print(
            "The output TSV file must end in `.tsv`.",
            file=sys.stderr
        )
        arg_error=True

    if arg_error:
        sys.exit(1)

    # make output directory as needed
    args.output_file.parent.mkdir(parents = True, exist_ok = True)


    # Run scrublet and export the results
    adata = anndata.read_h5ad(args.input_anndata_file)
    scrub_results = run_scrublet(adata, args.random_seed)
    scrub_results.to_csv(args.output_file, sep="\t", index=False )

if __name__ == "__main__":
    main()
