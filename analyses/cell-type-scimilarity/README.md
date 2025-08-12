# Annotating cell types with `SCimilarity`

## Description

This module will be used to annotate cell types in all ScPCA libraries using [`SCimilarity`](https://genentech.github.io/scimilarity/index.html). 

For cell type annotation, we use the model available on [Zenodo](https://zenodo.org/records/10685499). 

## Usage

You can run this module with the `run-analysis.sh` script present in this directory:

```sh
./run-analysis.sh
```

## Input files

This module requires the processed `SCE` objects for all projects as H5AD files. 

These files can be obtained with the following code run from the root of the `OpenScPCA-analysis` repository.
You must be logged into your AWS account to download these files.

```sh
# Download processed SCE objects as AnnData
./download-data.py --format AnnData
```

## Output files

TBD

## Software requirements

This module uses both `renv` and `conda` to manage R and Python software environments, respectively.

## Computational resources

Running `SCimilarity` requires at least 64 GB of RAM.  
