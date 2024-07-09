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

There are two portions of this module, which can be run as follows:

- Benchmarking analysis: Compare several doublet detection methods on non-ScPCA datasets with annotated doublets
  ```sh
  ./run_doublet-detection-benchmarking.sh
  ```

- ScPCA analysis: Detect doublets in ScPCA data using `scDblFinder` for a given project id (`SCPCPXXXXXX`)
  - Note that the conda environment is not used for this analysis, since it is only R-based
  ```sh
  ./run_doublet-detection-scpca.sh {scpca project id}

  # for example, detect doublets on all libraries in ScPCA project SCPCP000001:
  ./run_doublet-detection-scpca.sh SCPCP000001
  ```



## Input files

### Benchmarking input files

The benchmarking portion of this module uses input data from [a Zenodo repository](https://doi.org/10.5281/zenodo.4562782) to explore doublet detection methods.
Specifically, these datasets are used (descriptions from [Xi and Li (2021)](https://doi.org/10.1016/j.cels.2020.11.008)):
- `hm-6k`
  - A mixture of human `HEK293T` and mouse `NIH3T3` cells with 6806 droplets
  - Droplets were annotated as a doublets if the barcode was associated with both human and mouse
- `HMEC-orig-MULTI`
  - Human primary mammary epithelial cells (`HMECs`) with 26426 droplets
  - Doublets were annotated with the [`MULTI-seq` pipeline](https://github.com/chris-mcginnis-ucsf/MULTI-seq)
- `pdx-MULTI`
  - A mixture of human breast cancer cells and mouse immune cells from a PDX mouse model
  - Doublets were annotated with the [`MULTI-seq` pipeline](https://github.com/chris-mcginnis-ucsf/MULTI-seq)
- `pbmc-1B-dm`
  - `PBMCs` from a patient with systemic lupus erythematosus
  - Droplets were annotated with `demuxlet`

### ScPCA input files

The ScPCA portion of this module uses `processed` SCE files from most recent OpenScPCA data release, which can be obtained with:

```sh
# to be run from the root of the repository
./download-data.py
```


## Output files

Below is the results directory structure, annotated with file descriptions, after both `run_doublet-detection-benchmarking.sh` and `run_doublet-detection-scpca.sh` have been run:

```
results
├── README.md
├── benchmark-results
│   ├── {dataset}_scdblfinder.tsv
│   ├── {dataset}_scrublet.tsv
│   └── rendered-notebooks # Exploration of doublet detection results for each individual dataset
│       └── {dataset}-doublet-results.nb.html
└── scpca-results # Results from running scDblFinder across ScPCA projects
    └── {project id}
        └── {sample id}
            └── {library id}_processed_scdblfinder.tsv
```

## Software requirements

This module uses both `renv` and `conda` to manage software dependencies.
A Dockerfile is also provided.

## Computational resources

This module does not require compute beyond what is generally available on a laptop.
