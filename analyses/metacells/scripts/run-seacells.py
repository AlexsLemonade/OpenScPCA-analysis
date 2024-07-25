#!/usr/bin/env python3

"""
SEAcells Analysis Script
Joshua Shapiro
2024-07-18

This script runs the SEACells algorithm on a dataset.
"""


import argparse
import contextlib
import pathlib
import pickle
import sys

import anndata
import numpy as np
import scanpy as sc
import SEACells
import session_info


def convert_adata(adata: anndata.AnnData) -> anndata.AnnData:
    """
    Convert an ScPCA AnnData object to the formatting expected by scanpy and SEACells.
    """
    # put the highly variable genes in a column of adata.var
    adata.var["highly_variable"] = adata.var.gene_ids.isin(
        adata.uns["highly_variable_genes"]
    )

    # recompute principal components and umap (adds additional metadata)
    sc.tl.pca(adata, n_comps=50, mask_var="highly_variable")
    sc.pp.neighbors(adata)
    sc.tl.umap(adata)

    return adata


def run_seacells(
    adata: anndata.AnnData, cell_ratio: float = 75, verbose: bool = False
) -> tuple[anndata.AnnData, SEACells.core.SEACells]:
    """
    Run the SEACells algorithm on the given dataset.

    Parameters
    ----------
    adata : anndata.AnnData
        An AnnData object containing the data to run the SEACells algorithm on.
        Should contain an X_pca field with the PCA coordinates of the cells.
    cell_ratio : float
        The ratio of cells to metacells to use; i.e. number of cells per metacell
    verbose : bool
        Whether to print verbose output during the SEACells algorithm

    Returns
    -------
    anndata.AnnData
        The input AnnData object with the metacell assignments added to the obs table with the key "SEACell"
    SEACells.core.SEACells
        The SEACells model object
    """
    n_metacells = round(adata.n_obs / cell_ratio)
    n_eigs = 10  # number of eigenvalues for initialization
    # initialize the SEACells model
    model = SEACells.core.SEACells(
        adata,
        build_kernel_on="X_pca",
        n_SEACells=n_metacells,
        n_waypoint_eigs=n_eigs,
        convergence_epsilon=1e-5,
        verbose=verbose,
    )
    # initialize and fit model
    model.construct_kernel_matrix()
    model.initialize_archetypes()
    model.fit(min_iter=10, max_iter=50)

    return (adata, model)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Run the SEACell algorithm on the given dataset."
    )
    parser.add_argument(
        "adata_file", type=pathlib.Path, help="The input data in H5AD format."
    )
    parser.add_argument(
        "--adata_out",
        type=pathlib.Path,
        required=True,
        help="The output file path for the AnnData object (should end in .h5ad).",
    )
    parser.add_argument(
        "--model_out",
        type=pathlib.Path,
        required=False,
        help="The output file path for the SEACells model object.",
    )
    parser.add_argument(
        "--seed",
        type=int,
        help="The random seed to use for reproducibility.",
        default=2024,
    )
    parser.add_argument(
        "--logfile", type=pathlib.Path, help="File path for log outputs"
    )

    args = parser.parse_args()

    # check filenames
    if args.adata_out.suffix != ".h5ad":
        raise ValueError("Output file must end in .h5ad")
    if args.model_out and args.model_out.suffix != ".pkl":
        raise ValueError("Model output file must end in .pkl")

    if args.logfile:
        logs = open(args.logfile, "w")
        sys.stdout = logs

    # set seed for reproducibility
    np.random.seed(args.seed)

    adata = anndata.read_h5ad(args.adata_file)

    adata = convert_adata(adata)
    adata, seacell_model = run_seacells(adata, verbose=args.logfile is not None)

    # save the results
    print(f"Saving results to {args.adata_out}")
    adata.write_h5ad(args.adata_out, compression="gzip")

    if args.model_out:
        print(f"Saving results to {args.model_out}")
        with open(args.model_out, "wb") as f:
            pickle.dump(seacell_model, f)

    # only write session info if a logfile is provided
    if args.logfile:
        session_info.show(dependencies=True)
        logs.close()


if __name__ == "__main__":
    main()
