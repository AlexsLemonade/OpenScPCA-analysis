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

- Benchmarking analysis: Compare several doublet detection methods on publicly available (non-ScPCA) datasets with annotated doublets.
This analysis runs several doublet detection methods and assesses and compares their performance.
  ```sh
  ./run_doublet-detection-benchmarking.sh

  # to run in test mode using subsetted datasets:
  test=1 ./run_doublet-detection-benchmarking.sh

  # the script uses 4 cores by default. Use the cores argument to customize:
  cores=2 ./run_doublet-detection-benchmarking.sh
  ```

- ScPCA analysis: Detect doublets in ScPCA data using `scDblFinder` for a given project id (formatted as `SCPCPXXXXXX`).
  - The conda environment is not needed for this analysis, since it is only R-based
  ```sh
  ./run_doublet-detection-scpca.sh {scpca project id}

  # the script uses 4 cores by default, but you can use the cores argument to customize this.
  # for example, detect doublets on all libraries in ScPCA project SCPCP000001 using 2 cores:
  cores=2 ./run_doublet-detection-scpca.sh {scpca project id}
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

The ScPCA portion of this module uses `processed` SCE files from the most recent OpenScPCA data release, which can be obtained with:

```sh
# to be run from the root of the repository
./download-data.py
```


## Output files

Below is the results directory structure, annotated with file descriptions, after both `run_doublet-detection-benchmarking.sh` and `run_doublet-detection-scpca.sh` have been run.

```
results
├── README.md
├── benchmark-results
│   ├── {dataset}_scdblfinder.tsv # TSV files with `scDblFinder` inferences
│   ├── {dataset}_scrublet.tsv    # TSV files with `Scrublet` inferences
│   └── rendered-notebooks # Exploration of doublet detection results for each individual dataset
│       ├── compare-doublet-results.nb.html # Created from exploratory-notebooks/03_compare-benchmark-results.Rmd
│       └── {dataset}-doublet-results.nb.html # Created from template-notebooks/02_explore-benchmark-results.Rmd
└── scpca-results # Results from running `scDblFinder` across ScPCA projects
    └── {project id}
        └── {sample id}
            └── {library id}_processed_scdblfinder.tsv # TSV file with doublet results
```

### Result TSV files

The benchmarking TSV files have the following columns:

- `{dataset}_scdblfinder.tsv`: `barcodes`, `score`, `class`, `cxds_score`
    - `score` and `class` are the `scDblFinder` score and prediction columns, respectively, and `cxds_score` is a modified version of the [`scds::cxds` score](https://bioconductor.org/packages/devel/bioc/vignettes/scds/inst/doc/scds.html) as calculated by `scDblFinder`
- `{dataset}_scrublet.tsv`: `barcodes`, `scrublet_score`, `scrublet_prediction`

The ScPCA TSV files for `scDblFinder` results have only the columns `barcodes`, `score`, and `class`.
Note that SCEs with fewer than 10 droplets will not be run through `scDblFinder`, so their associated result TSVs will contain `NA` values in the `score` and `class` columns.

### Additional files

The directory `benchmark-test-data` contains the file `benchmark-test-data/data.zip`, a zipped directory of subsetted versions of the datasets from Xi and Li (2021) used for benchmarking.
Created with `benchmark-test-data/generate-benchmark-test-data.R`, these subsetted datasets are used when testing the benchmarking script `run_doublet-detection-benchmark.sh`.

## Software requirements

This module uses both `renv` and `conda` to manage software dependencies.
A Dockerfile is also provided.

## Computational resources

This module does not require compute beyond what is generally available on a laptop.
By default, `scDblFinder` is run with 4 cores.
To specify a different number of cores, use the `cores` argument when running either `run_doublet-detection-benchmark.sh` or `run_doublet-detection-scpca.sh`:

```sh
# For example, use 2 cores:

# benchmarking script
cores=2 ./run_doublet-detection-benchmarking.sh

# ScPCA script:
cores=2 ./run_doublet-detection-scpca.sh {scpca project id}
```

