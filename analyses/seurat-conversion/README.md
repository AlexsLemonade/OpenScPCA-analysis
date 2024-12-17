# Seurat conversion

This analysis module is used to convert the SingleCellExperiment objects from the ScPCA Portal to Seurat objects.

## Description

Please provide a description of your module, including:

- What type of analysis is included?
- What methods or tools are used?

If there are multiple steps in the module, please include an outline of the analysis steps and scripts used.

## Usage

TBD

## Input files

Input modules will be a SingleCellExperiment objects (as `.rds` files) from the ScPCA Portal.


## Output files

Seurat object files (as `.rds` files) will be generated as output, organized by project.

## Software requirements

The code for the project will be in R, so all software dependencies will be managed by `renv`.

## Computational resources

The only computational requirements should be sufficient RAM to load the SingleCellExperiment objects and resulting Seurat objects into memory.
