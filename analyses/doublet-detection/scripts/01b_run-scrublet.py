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
        "--input_anndata_file",
        type=str,
        required=True,
        help="Path to the input AnnData file to process."
    )
    parser.add_argument(
        "--results_dir",
        type=Path,
        required=True,
        help="The directory to export TSV file with doublet inferences."
    )
    parser.add_argument(
        "--random_seed",
        type=int,
        default=2024,
        help="Random seed to ensure reproducibility."
    )
    args = parser.parse_args()

    # Define and check files, directories
    args.results_dir.mkdir(parents = True, exist_ok = True)

    input_anndata_file = Path(args.input_anndata_file)
    if not input_anndata_file.exists():
        print(
            "The input AnnData file could not be found at:",
            args.input_anndata_file,
            file=sys.stderr
        )
        sys.exit(1)
    result_tsv = args.results_dir / (input_anndata_file.name.replace(".h5ad", "_scrublet.tsv"))


    # Run scrublet and export the results
    adata = anndata.read_h5ad(input_anndata_file)
    scrub_results = run_scrublet(adata, args.random_seed)
    scrub_results.to_csv(result_tsv, sep="\t", index=False )

if __name__ == "__main__":
    main()
