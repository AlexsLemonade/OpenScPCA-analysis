# `doublet-detection` analysis module

## Description

This module explores the use of doublet detection methods across ScPCA datasets.

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

Below is the results directory structure, annotated with file descriptions:
```
results
└── benchmark-results
    ├── {dataset}_scdblfinder.tsv # TSV files with `scDblFinder` inferences
    ├── {dataset}_scrublet.tsv    # TSV files with `scrublet` inferences
    └── rendered-notebooks
        └── {dataset}-results.nb.html # Exploration of doublet detection results for each individual dataset
```

## Software requirements

This module uses both `renv` and `conda` to manage software dependencies.
TODO: NEEDS UPDATING! A Dockerfile created using [these guidelines](https://openscpca.readthedocs.io/en/latest/software-platforms/docker/docker-images/#r-based-images) is also provided.

## Computational resources

_Forthcoming._
