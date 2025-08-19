# `infercnv-consensus-cell-type`

## Description

This analysis module holds various analyses that explore approaches for using consensus cell types as a normal reference for `inferCNV`.
Each analysis is described below.

### Pooled references

One analysis this module holds is an exploration of using normal references pooled across samples for a given project.
Project-specific scripts to run this analysis are available in [`pooled-workflows/`](./pooled-workflows) in this module.
You can run them all with:

```sh
./run-pooled-workflows.sh
```


#### Analysis results: Ewing sarcoma (`SCPCP000015`)

* Among references `endo`, `immune`, `endo-immune`, the `endo`-only reference appears more able to distinguish between tumor and normal cells, as observed in [this exploratory notebook](exploratory-notebooks/02_ewings-reference-cnv.Rmd)
  * Note that tumor cells were identified by the `cell-type-ewings` module
* Pooled references appear to do a similar or slightly better job compared to internal references at distinguishing between tumor and normal cells, as observed in [this exploratory notebook](exploratory-notebooks/03_ewings-pooled-internal.Rmd)


#### Input files

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
