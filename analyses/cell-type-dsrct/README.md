# DSRCT Cell Type analysis module



## Description

This analysis aims to annotate the DSRCT (Desmoplastic small round cell tumor) samples in the SCPCP000013 (n=7) dataset.

To do so, we will:

- Assess the quality of the data set
- Curate a list of marker genes associated with DSRCT cells
- Detect the expression levels of DSRCT marker genes in the data set
- Identify clusters in the DSRCT samples
- Use the list of marker genes to identify tumor cells
- Perform copy number inference to identify tumor cells
- Annotate normal cells
- Perform clustering on tumor cells and identify tumor cell states

## Usage

The code for the module will be in the form of a notebook.

## Input files


The input is dependent on the output files run from `download-data.py`.

```
./download-data.py --projects SCPCP000013
```

This downloads the `SingleCellExperiment` files for DSRCT samples.


## Output files

Please include a description of the output from your analysis, including:

- Plots from each analysis
- RDS file containing the processed single cell data set

## Software requirements

The analysis will be done in R using the `Seurat`, `SingleCellExperiment`, and `scran` packages.

## Computational resources

This will be done on a local machine.
