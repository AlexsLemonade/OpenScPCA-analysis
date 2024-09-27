# Glioblastoma cell type identification

## Description

This module aims to identify the cell types present and provide labels for 16 pediatric Glioblastoma samples. This is split into the following objectives:
- [ ] Visualise the data quality of each sample using standard quality control metrics (gene count, mitochondrial reads etc.)
- [ ] Create a linear regression model using CellTypist with Glioblastoma atlas 'GBMap' (Core) as reference data
- [ ] Apply the linear regression model using CellTypist to each sample individually
- [ ] Use the pre-generated UMAP to cluster the cells and isolate clusters with low model confidence
- [ ] Manually label the low-confidence clusters using known marker genes and differential expression analysis
- [ ] Provide final cell labels for each sample 

## Usage

In progress

## Input files

The input is dependent on the output fles run from `download-data.py`. 

```
./download-data.py --projects SCPCP00001
```

This downloads the `SingleCellExperiment` files for Glioblastoma samples.


## Output files

In progress

## Software requirements

The analysis will be done in Python version 3.11 primarily using the `scanpy' and `CellTypist' packages. 

## Computational resources

TBC
