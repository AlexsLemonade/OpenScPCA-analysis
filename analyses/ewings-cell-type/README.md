# Cell typing of Ewing sarcoma samples

## Description

This module will include code to annotate cell types in the Ewing sarcoma samples from SCPCP000015 present on the ScPCA Portal.

Currently this module includes analysis notebooks used to explore cell typing methods.

1. `01-marker-gene-classification-cellassign.Rmd`: This notebook looks at marker gene expression in SCPCS000490.
Tumor cells are classified using manual annotation based on marker gene expression and clustering.
`CellAssign` is also used to classify cells with the references stored in `references/celassign_refs`.

## Usage

Please provide instructions on how to run the analysis module.
What commands are needed to execute all steps in the analysis?

## Input files

This module requires the processed `SingleCellExperiment` objects (`_processed.rds`) and processed `AnnData` objects (`.hdf5` files) from SCPCP0000015.
These files were obtained using the `download-data.py` script:

```sh
# download SCE objects
./download-data.py --projects SCPCP000015
# download AnnData objects
./download-data.py --projects SCPCP000015 --format "AnnData"
```
This module also requires `references/tumor_marker_genes.tsv` which contains a list of marker genes for identifying Ewing sarcoma tumor cells.

## Output files

Please include a description of the output from your analysis, including:

- What type of files are created?
- What are the contents of the output files?
- Where are the files stored?
- Are any intermediate files generated?
If so, where are they stored?

## Software requirements

Please describe the environment system used to manage software dependencies for the module (e.g., `conda`, `renv`, `docker`).
Include the locations of any files used to set up the environment.

## Computational resources

Please indicate the computational resources that are required to run the module.
If there are specific memory, CPU, and/or GPU requirements, please state that.
