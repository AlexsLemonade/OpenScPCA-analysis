#!/usr/bin/env python3

"""
SEAcells Analysis Script
Joshua Shapiro

This script runs the SEACells algorithm on a dataset.
"""

import argparse
import datetime
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
    # check if PCA and metadata were already calculated with scanpy
    pca_present = (
        "X_pca" in adata.obsm and "pca" in adata.uns and "variance" in adata.uns["pca"]
    )

    # recalculate PCA if not present
    if not pca_present:
        if "X_pca" in adata.obsm:
            n_comps = adata.obsm["X_pca"].shape[1]
        elif "X_PCA" in adata.obsm:
            n_comps = adata.obsm["X_PCA"].shape[1]
        else:
            # set PCA to default 50, or min dimension if smaller
            n_comps = min(50, adata.n_vars, adata.n_obs)
        sc.tl.pca(adata, n_comps=n_comps, mask_var="highly_variable")

    # always recalculate UMAP (for now, at least)
    sc.pp.neighbors(adata)
    sc.tl.umap(adata)

    return adata


def run_seacells(
    adata: anndata.AnnData,
    cell_ratio: float = 75,
    min_cells: int = 20,
    verbose: bool = False,
) -> tuple[anndata.AnnData, SEACells.core.SEACells]:
    """
    Run the SEACells algorithm on the given dataset.

    Parameters
    ----------
    adata : anndata.AnnData
        An AnnData object containing the data to run the SEACells algorithm on.
    cell_ratio : float
        The ratio of cells to metacells to use; i.e. number of cells per metacell
        Default is 75, based on recommentations in https://github.com/dpeerlab/SEACells/blob/3462c624ffae0df6d3930490f345f00196c3503e/notebooks/SEACell_computation.ipynb
    min_cells : int
        The minimum number of cells for the SEACells algorithm to run,
        and the minimum number of metacells to create
    verbose : bool
        Whether to print verbose output during the SEACells algorithm

    Returns
    -------
    anndata.AnnData
        The input AnnData object with the metacell assignments added to the obs table with the key "SEACell"
    SEACells.core.SEACells
        The SEACells model object
    """

    if adata.n_obs < min_cells:
        raise ValueError(f"The dataset must have at least {min_cells} cells to run SEACells")
    # reformat the data for compatibility with SEACells and scanpy downstream
    adata = convert_adata(adata)

    n_metacells = max(round(adata.n_obs / cell_ratio), min_cells)
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
    model.fit(min_iter=10, max_iter=100)  # default iteration values

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
        start_time = datetime.datetime.now()
        print("Start time:", start_time)

    # set seed for reproducibility
    np.random.seed(args.seed)

    adata = anndata.read_h5ad(args.adata_file)
    try:
        adata, seacell_model = run_seacells(adata, verbose=args.logfile is not None)
    except ValueError as e:
        print(f"Error processing {args.adata_file}:", e, file=sys.stderr)
        seacell_model = None

    # save the results
    adata.write_h5ad(args.adata_out, compression="gzip")

    if args.model_out and seacell_model:
        with open(args.model_out, "wb") as f:
            pickle.dump(seacell_model, f)

    # only write session info if a logfile is provided
    if args.logfile:
        session_info.show(dependencies=True)
        end_time = datetime.datetime.now()
        print("End time:", end_time)
        print("Total time:", (end_time - start_time))
        logs.close()


if __name__ == "__main__":
    main()
