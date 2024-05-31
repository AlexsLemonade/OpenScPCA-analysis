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

### Benchmarking input files

This module currently uses input data from [a Zenodo repository](https://doi.org/10.5281/zenodo.4562782) to explore doublet detection methods.
Specifically, these datasets are used (descriptions from [Xi and Li (2021)](https://doi.org/10.1016/j.cels.2020.11.008)):
- `hm-6k`
  - A mixture of human HEK293T and mouse NIH3T3 cells with 6806 droplets
  - Droplets were annotated as a doublets if the barcode was associated with both human and mouse
- `HMEC-orig-MULTI`
  - Human primary mammary epithelial cells (HMECs) with 26426 droplets
  - Doublets were annotated with the [`MULTI-seq` pipeline](https://github.com/chris-mcginnis-ucsf/MULTI-seq)
- `pdx-MULTI`
  - A mixture of human breast cancer cells and mouse immune cells from a PDX mouse model
  - Doublets were annotated with the [`MULTI-seq` pipeline](https://github.com/chris-mcginnis-ucsf/MULTI-seq)
- `pbmc-1B-dm`
  - PMBCs from a patient with systemic lupus erythematosus
  - Droplets were annotated with `demuxlet`


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
