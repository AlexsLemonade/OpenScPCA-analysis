# cell-type-neuroblastoma-04

The goal of this analysis module is to perform cell type annotation on samples from `SCPCP000004`.

## Description

TBD.

## Usage

TBD.

## Input files

This module requires the following input files:

* The processed `SCE` and `AnnData` objects for project `SCPCP000004`
* The merged `SCE` object for `SCPCP000004`
* Consensus cell types for `SCPCP000004`

These files can be obtained with the following code run from the root of the `OpenScPCA-analysis` repository.
You must be logged into your AWS account to download these files.

```sh
# Download processed SCE and AnnData objects
./download-data.py --projects SCPCP000004 --format SCE,AnnData

# Download merged SCE objects
./download-results.py --modules merge-sce --projects SCPCP000004

# Download consensus cell types
./download-results.py --modules cell-type-consensus --projects SCPCP000004
```


## Output files

TBD.

## Software requirements

This module uses both `renv` and `conda` to manage R and Python software environments, respectively.

## Computational resources

TBD.
