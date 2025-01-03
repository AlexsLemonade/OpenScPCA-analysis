# Identifying consensus cell types

This module explores creating rules that can be used to identify a consensus cell type label.
Specifically, the cell type annotations obtained from both `SingleR` and `CellAssign` will be used to create a single cell type label in an ontology aware manner.

## Description

The goal of this module is to create a reference that can be used to define an ontology aware consensus cell type label for all cells across all ScPCA samples. 
This module performs a series of steps to accomplish that goal: 

1. The cell type annotations present in the `PanglaoDB` reference file were assigned to an ontology term identifier, when possible. 
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


## Usage

See the [`scripts/README.md`](./scripts/README.md) for instructions on running the scripts in this module. 

## Input files

TBD

## Output files

TBD

## Software requirements

TBD

## Computational resources

TBD
