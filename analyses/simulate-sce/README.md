# Simulate SCE objects

This module is for creating simulated SCE and AnnData files from existing data.

## Requirements

The R requirements for this module are described in the `renv.lock` file.
Before running the scripts here, open R in this directory and run the following command to install the required dependencies:

```r
renv::restore()
```

The AnnData conversion requires the `anndata` and `pandas` packages, which are included in the `openscpca-simulate-sce` conda environment.
To install and activate this environment, use the following commands:

```bash
conda-lock install -n openscpca-simulate-sce
conda activate openscpca-simulate-sce
```

### Data download

Data from projects to be processed should be downloaded from the OpenScPCA release bucket using the `download-data.py` script in the repository root directory.

You will want all types of data for this module, so you can run the following command to download all processing levels in both formats from the latest release:

```bash
./download-data.py --process-stage "filtered,unfiltered,processed" --format "sce,anndata"
```

## Running the scripts

The main entry point is the `simulate-project.sh` script.
This takes a single argument, the id of a project that is present in the `data/current` directory.

All files for the project will be simulated and saved to the `results/simulated` directory.

```bash
./simulate-project.sh SCPCP000015
```

This will first create a permuted metadata file for the project.
Then the `scripts/simulate-sce.R` script will be called for each sample to create simulated `.rds` files.

## Simulation description

Data is simulated using the [`splatter` package](https://bioconductor.org/packages/3.19/bioc/html/splatter.html)([Zappia _et. al._ 2017](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-017-1305-0)), using the simple simulation framework.

Most metadata is preserved from the original data to facilitate use in downstream analyses, though note that the metadata may not match the contents of the simulated data, leading to unusual results.
For details on the metadata that is recalculated, preserved, randomized or removed, see [Simulation and metadata in SCE objects](simulation-metadata.md).
