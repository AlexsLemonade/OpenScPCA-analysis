# Metacell Analysis Module

This module explores the use of metacells in the ScPCA context.
This will include both processing of the data and some comparisons among different methods.

## Description

Initially, this module will contain some scripts and notebooks for applying analysis methods to the data.

## Usage

To run the module, first install and activate the `openscpca-metacells` conda environment.
From the module directory, run:

```sh
conda-lock install --name openscpca-metacells
conda activate openscpca-metacells
```

Then run the main bash script:

```sh
bash run-metacells.sh
```


## Input files

Initial analysis is expected use processed `SingleCellExperiment` objects from the ScPCA Portal, in AnnData format.

For initial testing, files were downloaded with:

```sh
./download-data.py --format "anndata" --projects SCPCP000001
```

## Output files

Currently, the script outputs three files per library: an updated AnnData object, a pickled SEACells model object, and a log file.

## Software requirements

Currently, software requirements are handled by conda/conda-lock, with the `environment.yml` and `conda-lock.yml` files in this directory and a conda environment named `openscpca-metacells`

## Computational resources

TBD
