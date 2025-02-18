# Identifying consensus cell types

This module explores creating rules that can be used to identify a consensus cell type label.
Specifically, the cell type annotations obtained from both `SingleR` and `CellAssign` will be used to create a single cell type label in an ontology aware manner.

## Creating a reference for consensus cell types 

The goal of this module is to create a reference that can be used to define an ontology aware consensus cell type label for all cells across all ScPCA samples. 
This module performs a series of steps to accomplish that goal: 

1. The cell type annotations present in the `PanglaoDB` reference file were assigned to an ontology term identifier, when possible.
See [`references/README.md`](./references/README.md) for a full description on how we completed assignments.  
2. We looked at all possible combinations of cell type labels between the `PanglaoDB` reference (used with `CellAssign`) and the `BlueprintEncodeData` reference (used with `SingleR`). 
We then explored using a set of rules used to define consensus cell types in [`exploratory-notebooks/01-reference-exploration.Rmd`](./exploratory-notebooks/01-reference-exploration.Rmd). 
3. We created a [reference table](./references/consensus-cell-type-reference.tsv) containing all combinations for which we were able to identify a consensus cell type label.  
The consensus cell type label corresponds to the [latest common ancestor (LCA)](https://rdrr.io/bioc/ontoProc/man/findCommonAncestors.html) between the `PanglaoDB` and `BlueprintEncodeData` terms. 

When creating the consensus cell type labels we implemented the following rules: 

- If the terms share more than 1 LCA, no consensus label is set. 
The only exception is if one of the LCA terms corresponds to `hematopoietic precursor cells`. 
If that is the case all other LCA terms are removed and `hematopoietic precursor cell` is used as the consensus label. 
- If the LCA has greater than 170 descendants, no consensus label is set, with some exceptions: 
  - When the LCA is `neuron`, `neuron` is used as the consensus label. 
  - When the LCA is `epithelial cell` and the annotation from `BlueprintEncodeData` is `Epithelial cells`, then `epithelial cell` is used as the consensus label. 
  - If the LCA is `bone cell`, `lining cell`, `blood cell`, `progenitor cell`, or `supporting cell`, no consensus label is defined. 

See the [`scripts/README.md`](./scripts/README.md) for instructions on running the individual scripts used to generate the reference. 

## Assigning consensus cell types for ScPCA samples

The `assign-consensus-celltypes.sh` script can be used to assign a consensus cell type for all samples in a single ScPCA project. 
This script outputs a single TSV file for each library in ScPCA with cell type annotations for all cells in that library. 
Cell type annotations assigned using `SingleR` with the `BlueprintEncodeData` reference and `CellAssign` using the `PanglaoDB` reference are included along side the assigned consensus cell type annotation and ontology identifier. 

To run this script for a given project use the following command: 

```sh
./assign-consensus-celltypes.sh "SCPCP000001"
```

### Input files


The `assign-consensus-celltypes.sh` script requires the processed `SingleCellExperiment` objects (`_processed.rds`) for all ScPCA samples.
These files were obtained using the `download-data.py` script:

```sh
# download SCE objects
./download-data.py
```

This script also requires four reference files, `blueprint-mapped-ontologies.tsv`, `panglao-cell-type-ontologies.tsv`, `consensus-cell-type-reference.tsv`, and `validation-markers.tsv`. 
See [Creating a reference for consensus cell types](#creating-a-reference-for-consensus-cell-types) and the [README.md in the references directory](./references/README.md) to learn more about the content of these files. 

### Output files

Running the `assign-consensus-celltypes.sh` script will generate two TSV files for each library. 
Output files will be in `results/cell-type-consensus` and organized as follows: 

```
results
└── cell-type-consensus
    └── <project id>
        └── <sample id>
            ├── <library id>_consensus-cell-type-assignments.tsv.gz
            └── <library id>_marker-gene-expression.tsv.gz

```

The `<library id>_consensus-cell-type-assignments.tsv.gz` file contains cell type annotations for all cells in a single library with the following columns: 

| | |
| --- | --- | 
| `project_id` | ScPCA project id |
| `sample_id` | ScPCA sample id | 
| `library_id` |  ScPCA library id |
| `barcodes` | cell barcode |
| `sample_type` | A string indicating the type of sample, with one of the following values: `"patient-derived xenograft"`, `"cell line"`, or `"patient tissue"`. If `cell line`, no cell type annotation is performed and all columns will have `NA` |
| `singler_celltype_ontology` | Cell type ontology term assigned by `SingleR` | 
| `singler_celltype_annotation` | Name associated with cell type ontology term assigned by `SingleR`; this term is equivalent to the `label.main` term in the `BlueprintEncodeData` reference | 
| `cellassign_celltype_annotation` | Cell type assigned by `CellAssign`; this term is the original term found in the `PanglaoDB` reference file | 
| `panglao_ontology` | Cell type ontology term associated with the term found in `cellassign_celltype_annotation` column | 
| `panglao_annotation` | Name associated with the cell type ontology term in `panglao_ontology` | 
| `blueprint_annotation_cl` | Name associated with the cell type ontology term in `singler_celltype_ontology` |
| `consensus_ontology` | Cell type ontology term assigned as the consensus cell type | 
| `consensus_annotation` | Name associated with the assigned consensus cell type in `consensus_ontology` | 


The `<library id>_marker-gene-expression.tsv.gz` file contains the `logcounts` for all marker genes in [`references/validation-markers.tsv`](./references/validation-markers.tsv). 
Only genes that are expressed in the library are included in the output. 

| | | 
| --- | --- | 
| `library_id` | ScPCA library id | 
| `barcodes` | cell barcode | 
| `ensembl_gene_id` | Ensembl gene identifier for marker gene | 
| `gene_symbol` | Gene symbol for marker gene | 
| `gene_expression` | Gene expression for marker gene obtained from the `logcounts` assay | 

## Software requirements

This module uses `renv` to manage software dependencies. 

## Computational resources

This module does not require compute beyond what is generally available on a laptop. 
