# cell-type-neuroblastoma-04

The goal of this analysis module is to perform cell type annotation on samples from `SCPCP000004`.

## Description

TBD.

## Usage

You can run this module with the `run-analysis.sh` script present in this directory:

```sh
./run-analysis.sh
```

To run with test data and a subsetted version of NBAtlas, use this command:
```sh
testing=1 ./run-analysis.sh
```

By default, this script uses the following settings:
* Performs reference aggregation before training the SingleR model
* Runs on all samples in the project (vs a designated subset)
* Uses 4 threads

See documentation in the script for how to modify these settings.

## Input files

This module requires the following input files:

* The processed `SCE` objects for project `SCPCP000004`
* The merged `SCE` object for `SCPCP000004`

These files can be obtained with the following code run from the root of the `OpenScPCA-analysis` repository.
You must be logged into your AWS account to download these files.

```sh
# Download processed SCE objects
./download-data.py --projects SCPCP000004 --format SCE

# Download merged SCE objects
./download-results.py --modules merge-sce --projects SCPCP000004
```


## Output files

TBD.

## Software requirements

This module uses both `renv` and `conda` to manage R and Python software environments, respectively.

## Computational resources

TBD.
