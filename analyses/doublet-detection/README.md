# `doublet-detection` analysis module

## Description

This module explores doublet detection across ScPCA datasets.

Methods used in this module include the following:

- [`scDBlFinder`](https://bioconductor.org/packages/release/bioc/html/scDblFinder.html)
- [`scrublet`](https://github.com/swolock/scrublet)


## Usage

To run this module, first create the `openscpca-doublet-detection` conda environment, and then activate it:

```sh
# create the environment
conda-lock install --name openscpca-doublet-detection conda-lock.yml

# activate the environment
conda activate openscpca-doublet-detection
```

Then, run the following bash script:

```sh
./run_doublet-detection.sh
```

## Input files

This module currently uses input data from [a Zenodo repository](https://doi.org/10.5281/zenodo.4562782) to explore doublet detection methods.
Specifically, these datasets are used: `hm-6k`, `pbmc-1B-dm`, `pdx-MULTI`, and `HMEC-orig-MULTI`.

Eventually, we'd like to run all ScPCA datasets through doublet detection, but this is still TBD for this specific module.

## Output files

- `results/benchmark_results`
    - `{dataset_name}_sce.tsv`: TSV files with `scDblFinder` inferences
    - `{dataset_name}_scrublet.tsv`: TSV files with `scrublet` inferences

## Software requirements

This module uses both `renv` and `conda` to manage software dependencies.
TODO: NEEDS UPDATING! A Dockerfile created using [these guidelines](https://openscpca.readthedocs.io/en/latest/software-platforms/docker/docker-images/#r-based-images) is also provided.

## Computational resources

_Forthcoming._
