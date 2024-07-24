#!/usr/bin/env python3

"""
Seacells Analysis Script
Joshua Shapiro
2024-07-18

This script runs the SEACells algorithm on a dataset.
"""

# Load modules
#
# Load required Python modules at the top of your script
# We have included the standard `pathlib` module and the `session_info` module
# that we will be using at the bottom of this notebook to record the versions of
# the modules used in this analysis.
#
# Do not install modules here; only load them with `import` statements.
# Avoid renaming modules with `as` statements, unless you are performing a
# standard renaming (e.g., `import pandas as pd`).

import argparse
import contextlib
import pathlib
import pickle

import anndata
import numpy as np
import scanpy
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
    scanpy.tl.pca(adata, n_comps=50, mask_var="highly_variable")
    scanpy.tl.umap(adata)

    return adata


def run_seacells(
    adata: anndata.AnnData, cell_ratio: float = 75
) -> tuple[anndata.AnnData, SEACells.core.SEACells]:
    """
    Run the SEACells algorithm on the given dataset.

    Parameters
    ----------
    adata : anndata.AnnData
        An annData object containing the data to run the SEACells algorithm on.
        Should contain an X_pca field with the PCA coordinates of the cells.
    cell_ratio : float
        The ratio of cells to metacells to use; i.e. number of cells per metacell

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
        "seed",
        type=int,
        help="The random seed to use for reproducibility.",
        default=2024,
    )
    parser.add_argument("logfile", type=pathlib.Path, help="File path for log outputs")

    args = parser.parse_args()

    # check filenames
    if args.adata_out.suffix != ".h5ad":
        raise ValueError("Output file must end in .h5ad")
    if args.model_out and args.model_out.suffix != ".pkl":
        raise ValueError("Model output file must end in .pkl")

    # set seed for reproducibility
    np.random.seed(args.seed)

    adata = anndata.read_h5ad(args.datafile)

    adata = convert_adata(adata)
    adata, model = run_seacells(adata)

    # save the results
    adata.write_h5ad(args.output, compression="gzip")

    if args.model_out:
        with open("args.model_out", "wb") as f:
            pickle.dump(model, f)

    # As the last step, record the versions of the modules and dependencies
    # that were used in this analysis
    with open(args.logfile, "w") as f:
        with contextlib.redirect_stdout(f):  # direct the session_info output to a file
            session_info.show(dependencies=True)


if __name__ == "__main__":
    main()
