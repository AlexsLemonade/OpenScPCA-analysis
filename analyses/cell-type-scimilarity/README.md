# Annotating cell types with `SCimilarity`

## Description

This module will be used to annotate cell types in all ScPCA libraries using [`SCimilarity`](https://genentech.github.io/scimilarity/index.html). 

For cell type annotation, we use the model available on [Zenodo](https://zenodo.org/records/10685499). 

## Usage

Prior to running `SCimilarity` on all samples, the `setup-analysis.sh` script will need to be used to download the model and create a reference file with ontology identifiers for all cell type annotations. 

This can be done by running the following command: 
```sh
./setup-analysis.sh
```


## Input files

This module requires the processed objects for all projects as H5AD files. 

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
