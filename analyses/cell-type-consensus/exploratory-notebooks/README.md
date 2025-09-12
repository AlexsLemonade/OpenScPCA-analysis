# Exploratory notebooks

This folder contains exploratory notebooks for this module.

1. `01-reference-exploration.Rmd`: This notebook was used to explore possible consensus label assignments between cell types in the `PanglaoDB` and `BlueprintEncodeData` references.
Observations made in this notebook were used to define the set of possible consensus labels to be included in [`references/consensus-cell-type-reference.tsv`](../references/consensus-cell-type-reference.tsv).

2. `02-explore-consensus-results.Rmd`: This notebook summarizes the consensus labels assigned to all ScPCA samples.
Prior to rendering this notebook results from the `cell-type-consensus` module in `OpenScPCA-nf` using the `2024-11-25` were downloaded.

3. `03-osteosarcoma-consensus-celltypes.Rmd`: This notebook summarizes the consensus labels assigned to projects with Osteosarcoma samples.

4. `04-brain-and-cns-consensus-celltypes.Rmd`: This notebook summarizes the consensus labels assigned to all samples belonging to projects with Brain and CNS tumors, excluding any samples that are part of multiplexed libraries.

5. `05-marker-gene-validation.Rmd`: This notebook looks at marker gene expression of markers in [`references/validation-markers.tsv`](../references/validation-markers.tsv) across consensus cell types for `SCPCP000001`.

6. `06-update-consensus-scimilarity.Rmd`: This notebook looks at updating the consensus cell type reference to incorporate all possible cell type annotations from using [`SCimilarity`](https://genentech.github.io/scimilarity/index.html). 

7. `07-compare-consensus.Rmd`: This notebook compares the consensus cell type labels assigned without `SCimilarity` to the newer consensus cell types assigned using input from `SCimilarity`. 

## Cell type validation notebooks

The `cell-type-validation-notebooks` contains the rendered reports summarizing expression of cell type marker genes in the consensus cell types.
There is one file for each project and the reports are rendered from the template notebook in `template-notebooks/marker-gene-validation.Rmd`.

To render all reports you will need to first download the results from the `cell-type-consensus` module in `OpenScPCA-nf` using the download results script from the root of this repository:

```sh
./download-results.py --module cell-type-consensus
```

Then you can render the reports using the `render-validation-notebooks.sh` script.

```sh
./render-validation-groups.sh
```

## Utils

The `utils` folder contains scripts with any functions that are used by multiple notebooks.
