# Seurat conversion

This analysis module is used to convert the SingleCellExperiment objects from the ScPCA Portal to Seurat objects.

## Description



## Usage

The `run_seurat-conversion.sh` script will convert all `_processed.rds` and `_filtered.rds` files in the repository's `data/current/` directory, placing output files in the module's `results/seurat/` directory.
To process specific sets of files, the `scripts/seurat-conversion.R` script can be run directly, specifying input and output directories as needed.


## Input files

Input files are SingleCellExperiment objects (as `.rds` files) from the ScPCA Portal.
By default, these are the files in the `data/current/` directory of the `OpenScPCA-analysis` repository, as downloaded by the `download-data.py` script.

## Output files

Seurat object files (as `.rds` files) will be generated as output, organized following the same directory structure as the input files.
For example, if an input file is `$REPO_ROOT/data/current/SCPCP000000/SCPCS000000/SCPCL000000_processed.rds`, the corresponding output file will be placed at `results/seurat/SCPCP000000/SCPCS000000/SCPCL000000_processed_seurat.rds`.

## Software requirements

The code for this module is all written in R, with all software dependencies managed by `renv`.

## Computational resources

The only computational requirements should be sufficient RAM to load the SingleCellExperiment objects and resulting Seurat objects into memory.
