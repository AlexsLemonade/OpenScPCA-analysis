# infercnv-consensus-cell-type

## Description

The goal of this analysis module is to explore approaches for normal, in particular immune, consensus cell types as a normal reference for `inferCNV`.

## Usage


## Input files

* Results from the `cell-type-consensus` module for project `SCPCP000015`
* The merged SCE for project `SCPCP000015`

To obtain these files, run the following commands from the top-level of the repository after logging into your AWS account.

```sh
# Obtain the consensus cell type results
./download-results.py --module cell-type-consensus --project SCPCP000015

# Obtain the merged SCE
./download-results.py --module merge-sce --project SCPCP000015
```

## Output files


## Software requirements

This module uses `renv` and `docker` to manage software dependencies.

## Computational resources

This module can be run with the resources of a standard personal laptop.

