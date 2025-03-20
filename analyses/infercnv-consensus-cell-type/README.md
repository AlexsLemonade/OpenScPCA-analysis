# infercnv-consensus-cell-type

## Description

The goal of this analysis module is to explore approaches for using normal, in particular immune, consensus cell types as a normal reference for `inferCNV`.

## Usage


## Input files

* Results from the `cell-type-ewings` module
* Processed SCEs for project `SCPCP000015`

To obtain these files, run the following commands from the top-level of the repository after logging into your AWS account.

```sh
# Download the cell-type-ewings results
./download-results.py --module cell-type-ewings

# Download the processed SCE objects
./download-data.py --project SCPCP000015
```

## Output files


## Software requirements

This module uses `renv` and `docker` to manage software dependencies.

## Computational resources

This module can be run with the resources of a standard personal laptop.


## Analysis information

To validate `inferCNV` results, we use a set of expected CNVs across diagnoses:

### Ewing sarcoma

* Gain of chr1q [1,2]
* Gain of chr8 [1,2]
* Loss of chr16q [1]
* Gain of chr12 [1]

Sources:
1. https://doi.org/10.1158/2159-8290.CD-14-0622
2. https://doi.org/10.1158/2159-8290.CD-13-1037
