#!/usr/bin/env python3

"""
Seacells Analysis Script
Joshua Shapiro
2024-07-18

This script runs the SEACells algorithm on the given dataset.
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


def run_seacells(adata: anndata.AnnData, cell_ratio: float = 75) -> anndata.AnnData:
    """
    Run the SEACells algorithm on the given dataset.
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

    return adata


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Run the SEACell algorithm on the given dataset."
    )
    parser.add_argument(
        "datafile", type=pathlib.Path, help="The input data in H5AD format."
    )
    parser.add_argument(
        "seed",
        type=int,
        help="The random seed to use for reproducibility.",
        default=2024,
    )
    parser.add_argument("output", type=pathlib.Path, help="The output file path.")
    parser.add_argument("logfile", type=pathlib.Path, help="File path for log outputs")

    args = parser.parse_args()

    # set seed for reproducibility
    np.random.seed(args.seed)

    adata = scanpy.read(args.datafile)

    adata = convert_adata(adata)
    adata = run_seacells(adata)

    # As the last step, record the versions of the modules and dependencies
    # that were used in this analysis
    with open(args.logfile, "w") as f:
        with contextlib.redirect_stdout(f):  # direct the session_info output to a file
            session_info.show(dependencies=True)


if __name__ == "__main__":
    main()
