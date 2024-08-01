# Exploratory Analysis

This folder contains all notebooks used to conduct exploratory analysis while cell typing.
These notebooks all use data from the `2024-05-01` release.

1. `01-marker-gene-tumor-classification.Rmd`: This notebook looks at marker gene expression in `SCPCS000490`.
Tumor cells are classified using manual annotation based on marker gene expression and clustering.

2. `02-cellassign-tumor-cell-classification.Rmd`: This notebook looks at using `CellAssign` to identify tumor cells using a variety of marker gene references.
Annotations from `CellAssign` are also compared to manual annotations identified in `01-marker-gene-tumor-classification.Rmd`.

3. `03-copykat.Rmd`: This notebook looks at using `CopyKAT` to annotate tumor and normal cells in `SCPCS000490`.
`CopyKAT` is run using the `scripts/run-copykat.R` script and outputs are evaluated in this notebook.
This includes comparing annotations to those obtained from marker gene annotation in `01-marker-gene-tumor-classification.Rmd`.

4. `04-infercnv.Rmd`: This notebook looks at using `InferCNV` to identify CNVs and annotate tumor and normal cells in `SCPCS000490`.
`InferCNV` is run using `scripts/run-infercnv.R` and outputs are evaluated in this notebook.
This includes comparing annotations to those obtained from marker gene annotation in `01-marker-gene-tumor-classification.Rmd`.

5. `05-cluster-exploration.Rmd`: This notebook looks at clustering across different parameters to choose the most optimal clusters for `SCPCL000822` and `SCPCL000824`. 
Additionally, we look at expression of marker genes from `references/visser-all-marker-genes.tsv` across all clusters and use that to refine tumor cell annotations obtained from running `aucell-singler-annotation.sh`. 

## Annotation notebooks

The `annotation_notebooks` folder contains exploratory notebooks used to explore and validate annotations in individual samples. 
