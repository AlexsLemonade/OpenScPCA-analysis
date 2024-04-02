# Hello R

This is an example analysis module using an R notebook to perform a simple analysis of the OpenScPCA data.

## Description

This analysis counts the number of cells in each processed SCE file in the current data directory, creating a table of cell counts and histograms of cell counts by project.

## Usage

To run the analysis, run the following command from this analysis directory:

```bash
bash run_hello-R.sh
```

## Input files

- `data/current/{projects}/**/*_processed.rds`: processed `SingleCellExperiment` objects from an ScPCA data release

These input files can be obtained by running the `download-data.py` script in OpenScPCA repository root directory with default arguments.
Requires AWS credentials with access to the OpenScPCA data release bucket.

## Output files

- `hello.nb.html`: a rendered R Markdown notebook
- `results/cell_counts.csv`: a table of counts for each library
- `plots/cell_counts.pdf`: histograms of cell counts by project

## Software requirements

This module runs in R, with packages managed by `renv`.

Before running the analysis script, open an R session in the analysis module directory and run `renv::restore()` to install the required packages.
All packages required are listed in the `renv.lock` file.

## Computational resources

Computational resources required are minimal, however the computer must have sufficient memory to open each SCE file, and the analysis may take some time to run depending on the number of processed SCE files in the data directory.
