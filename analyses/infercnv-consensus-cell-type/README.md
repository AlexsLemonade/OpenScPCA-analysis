# infercnv-consensus-cell-type

## Description

The goal of this analysis module is to explore approaches for using normal, in particular immune, consensus cell types as a normal reference for `inferCNV`.

## Usage

To run the analysis, use the following command:

```sh
./run-analysis.sh
```

## Input files

* Results from the `cell-type-ewings` module
* Processed SCEs and the merged SCE for project `SCPCP000015`

To obtain these files, run the following commands from the top-level of the repository after logging into your AWS account.

```sh
# Download the cell-type-ewings results
./download-results.py --module cell-type-ewings

# Download the merged SCE object
./download-results.py --module merge-sce --project SCPCP000015

# Download the processed SCE objects
./download-data.py --project SCPCP000015
```

## Output files


## Software requirements

This module uses `renv` and `docker` to manage software dependencies.

## Computational resources

This module can be run with the resources of a standard personal laptop.
