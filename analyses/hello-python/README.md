# Hello Python

This is an example analysis module using a Jupyter Notebook to perform a simple analysis of the OpenScPCA data.

## Description

This analysis counts the number of cells in each processed AnnData file in the current data directory, creating a table of cell counts and histograms of cell counts by project.

## Usage

The workflow is run using the following command from this analysis directory:

```bash
bash hello-python.sh
```

## Input files

- `data/current/{projects}/**/*_processed_rna.hdf5`: processed `SingleCellExperiment` objects from an ScPCA data release.

These input files can be obtained by running the `download-data.py` script in the `scripts` directory with the following arguments:

```bash
./download-data.py --format AnnData
```

Requires AWS credentials with access to the OpenScPCA data release bucket.

### Output files

- `hello.html`: a rendered version of the Jupyter Notebook
- `results/cell_counts.csv`: a table of counts for each library
- `plots/cell_counts.pdf`: histograms of cell counts by project

## Software requirements

General software requirements are listed in the conda `environment.yml` file for this module.

To run with locked dependencies, you will need to install and activate a conda environment using the `conda-lock` package.
From this analysis directory, run the following commands at the command line:

```bash
conda-lock install --name openscpca-hello-python conda-lock.yml
conda activate openscpca-hello-python
```

## Computational resources

Computational resources required are minimal, however the computer must have sufficient memory to open each AnnData file, and the analysis may take some time to run depending on the number of processed AnnData files in the data directory.
