#!/usr/bin/env python3
"""
This script takes an AnnData object and checks for the `logcounts` in layers.
If present, `logcounts` is moved to `X` and `X` (which has the raw counts) is moved to `raw.X`

In addition, any DataFrames in `obsm` are converted to ndarrays, highly variable genes are converted to a `var` column.
If a pca metadata file is found, PCA variance statistics and standard creation values are in the format expected by scanpy.

Adapted from https://github.com/AlexsLemonade/scpca-nf/blob/v0.8.5/bin/reformat_anndata.py
Changed to work on a directory of files instead of a single file; only operates on the processed data files
"""

import argparse
import pathlib

import anndata
import pandas as pd


def reformat_anndata(anndata_file, pca_metafile):
    # read in anndata
    adata = anndata.read_h5ad(anndata_file)

    # if logcounts is present
    if "logcounts" in adata.layers:
        # move X to raw.X by creating the raw adata
        adata.raw = adata
        # move logcounts to X and rename
        adata.X = adata.layers["logcounts"]
        adata.uns["X_name"] = "logcounts"
        del adata.layers["logcounts"]

    # convert DataFrames in obsm to ndarrays
    for key, value in adata.obsm.items():
        if isinstance(value, pd.DataFrame):
            adata.obsm[key] = value.to_numpy()

    # add pca adata to uns if pca_meta_file is provided in the format created by scanpy
    if pathlib.Path(pca_metafile).exists():
        pca_meta = pd.read_csv(pca_metafile, sep="\t", index_col=0)
        if (
            "variance" not in pca_meta.columns
            or "variance_ratio" not in pca_meta.columns
        ):
            raise ValueError(
                "`pca_meta_file` must contain columns `variance` and `variance_ratio`"
            )
        pca_object = {
            "param": {
                "zero_center": True,
                "use_highly_variable": False,
                "mask_var": None,
            },
            "variance": pca_meta["variance"].to_numpy(),
            "variance_ratio": pca_meta["variance_ratio"].to_numpy(),
        }
        adata.uns["pca"] = pca_object

    # export adata
    adata.write_h5ad(anndata_file, compression="gzip")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-d",
        "--dir",
        help="directory containing H5AD files and PCA metadaa",
        required=True,
    )

    args = parser.parse_args()

    # find all processed rna h5ad files in the directory, recursively
    anndata_files = list(pathlib.Path(args.dir).rglob("*_processed_rna.h5ad"))
    pca_files = [str(f).replace("_rna.h5ad", "_pca.tsv") for f in anndata_files]

    for anndata_file, pca_file in zip(anndata_files, pca_files):
        reformat_anndata(anndata_file, pca_file)
