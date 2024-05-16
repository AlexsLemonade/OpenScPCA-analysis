# `doublet-detection` analysis module

## Description

This module explores doublet detection across ScPCA datasets.

Methods used in this module include the following:

- [`scDBlFinder`](https://bioconductor.org/packages/release/bioc/html/scDblFinder.html)
- [`scds`](https://bioconductor.org/packages/release/bioc/html/scds.html)


## Usage

- `script/00_benchmark-methods.R` runs four methods for doublet detection across a selection of (non-multiplexed) ScPCA datasets of varying library sizes across all projects.
  - This script takes a single argument, `--cores`, specifying the number of cores to use when running `scDblFinder`.

## Input files

_Forthcoming._

Eventually, we'd like to run all ScPCA datasets through doublet detection, but this is still TBD for this specific module.

## Output files

- `results/benchmark-times_4-cores.tsv` is created by `script/00_benchmark-methods.R`, invoked as:
```bash
Rscript scripts/00_benchmark-methods.R --cores 4
```

## Software requirements

This module uses `renv` to manage software dependencies.
A Dockerfile created using [these guidelines](https://openscpca.readthedocs.io/en/latest/software-platforms/docker/docker-images/#r-based-images) is also provided.

## Computational resources

This analysis can be run on a laptop.
