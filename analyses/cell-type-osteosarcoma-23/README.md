## Description

This analysis module will -

- Conduct cell typing for all 52 samples in Group ID SCPCAB0027, Project SCPCP000023
- Run InferCNV iteratively on each sample -
  - Using the Gene order file from Genome Reference Consortium Human Build 38 patch release 13 (GRCh38.p13)
  - Using the latest consensus cell type annotation, and the validation groups to construct normal references
  - Use KMeans clustering algorithm to computationally cluster cells based on CNV score into "tumor" and "normal"
  - Update the dataset with this annotation for each sample

## Usage

Please provide instructions on how to run the analysis module.
What commands are needed to execute all steps in the analysis?

## Input files

Please include a description of the required inputs to run the analysis module (e.g., processed `SingleCellExperiment` objects from the ScPCA Portal).
If the input to this module is dependent on running any other scripts (e.g., `download-data.py`) or on the output of another module, please state that.
If possible, include examples of specific commands used to obtain data.

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
