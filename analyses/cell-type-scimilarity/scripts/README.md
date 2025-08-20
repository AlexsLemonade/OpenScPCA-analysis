This directory contains scripts used in the `cell-type-scimilarity` module.

* `01-assign-ontology-ids.R` is used to create a TSV file with all possible cell type annotations found in the `SCimilarity` model and an associated Cell Ontology identifier. 
The output from this script is saved in `references/scimilarity-mapped-ontologies.tsv`. 

* `02-run-scimilarity.py` is used to annotate all cells in a single library with `SCimilarity`. 
This script requires a `_processed_rna.h5ad` file as input and outputs a TSV with cell type annotations and the `min_dist` reported by `SCimilarity` for all cells. 
