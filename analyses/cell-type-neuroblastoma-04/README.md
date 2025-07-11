# cell-type-neuroblastoma-04

The goal of this analysis module is to perform cell type annotation on samples from `SCPCP000004`.

## Description

This module annotates cell types across samples in `SCPCP000004` using the [`NBAtlas` reference](https://doi.org/10.1016/j.celrep.2024.114804) (Bonine et al. 2024).

* First, the module performs cell type annotation using [`SingleR`](https://doi.org/10.1016/j.celrep.2024.114804) (Aran et al. 2019)

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

See documentation in the script for how to modify these or other settings.

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

The module outputs results to the `results` directory.
Please see `results/README.md` for additional details.

## Software requirements

This module uses both `renv` and `conda` to manage R and Python software environments, respectively.

## Computational resources

This module requires at least 16 GB available RAM for the `SingleR` portion of the analysis.
