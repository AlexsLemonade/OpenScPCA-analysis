# `infercnv-consensus-cell-type`

## Description

The goal of this analysis module is to explore approaches for using normal consensus cell types as a normal reference for `inferCNV`.

## Usage

To run the full analysis, use the following command:

```sh
./run-analysis.sh
```

This script will run all [project-specific workflows](./project-workflows/README.md) in this module.
You can alternatively run a single project in the module using one of these project-specific scripts.

## Input files

* Results from the `cell-type-ewings` module
* Processed SCEs and the merged SCE for project `SCPCP000015`

To obtain these files, run the following commands from the top-level of the repository after logging into your AWS account.

```sh
# Download the cell-type-ewings results
./download-results.py --module cell-type-ewings

# Download the merged SCE object
./download-results.py --module merge-sce --project SCPCP000015

# Download the processed SCE objects
./download-data.py --project SCPCP000015
```

## Output files

The module outputs results to the `results` directory, as described in [`results/README.md`](./results/README.md).

## Software requirements

This module uses `renv` and `Docker` to manage software dependencies.

## Computational resources

This module can be run with the resources of a standard personal laptop.

## Analysis results

This section provides brief conclusions for different projects which this module processed.

### Ewing sarcoma (`SCPCP000015`)

* Among references `endo`, `immune`, `endo-immune`, the `endo`-only reference appears more able to distinguish between tumor and normal cells
  * Note that tumor cells were identified by the `cell-type-ewings` module
* Pooled references appear to do a similar or slightly better job than do internal references at distinguishing between tumor and normal cells
