# `scripts` directory instructions

This directory contains all scripts used for cell typing Wilms tumor samples from SCPCP000006.


## `prepare-fetal-references.R`

This script is used to create the fetal kidney reference (Stewart et al) and the full fetal reference (Cao et al).

## `explore-cnv-methods.sh`

This script is used to run and explores 2 CNV inference methods `copykat` and `infercnv` for a selection of samples from SCPCP000006:

- `SCPCS000194`
- `SCPCS000179`
- `SCPCS000184`
- `SCPCS000205`
- `SCPCS000208`

It calls the following:

- `05_copyKAT.R`
- `06_inferCNV.R`
- `../cnv-exploratory-notebooks/05_copykat_exploration.Rmd`
- `../cnv-exploratory-notebooks/06_infercnv_exploration.Rmd`


### `05_copyKAT.R`

This script is used to run `copykat`.
Please note that `copykat` can take more than one hour per sample, even when using quite a lot of resources.
The script by default uses 16 cores, but the module script `00_run_workflow.sh` will run it with 32 cores unless otherwise specified.

For each sample and each condition (reference and distance), we saved in `results/{sample_id}/05_copykat/{reference}/{distance}/`:

- the final `copykat` object in `05_copykat_{sample_id}_{reference}_distance-{selection}.rds`
- the heatmap of CNV in `05_copykat_{sample_id}_{reference}_distance-{selection}_copykat_heatmap.jpeg`
- the prediction (aneuploid or diploid value per cell) in `05_copykat_{sample_id}_{reference}_distance-{selection}_copykat_prediction.txt`
- the CNA matrix `05_copykat_{sample_id}_{reference}_distance-{selection}_copykat_CNA_results.txt`


### `06a_build-geneposition.R`

This script builds the gene position file that will be used in `06_infercnv.R` for each of the samples.
We build the position file once via downloading the `gencode_v19_gen_pos.complete.txt` from the [Trinity/`CTAT`](https://data.broadinstitute.org/Trinity/CTAT/cnv/gencode_v19_gen_pos.complete.txt).

The `gencode_v19_gen_pos.complete.txt` is then saved in `results/references`.



### `06_infercnv.R`

This script is used to run `infercnv`.
We test for the sensitivity of `infercnv` in regards to the selection of normal cells used as reference.
Therefore, the parameter `reference` indicates if we want to select immune and/or endothelial cells as reference or no cell at all.

For each sample and each condition (reference), we saved in `results/{sample_id}/06_infercnv/reference-{selection}`:
- the final `infercnv` object in `06_infercnv_{sample_id}_reference-{selection}.rds`
- the heatmap of CNV in `06_infercnv_{sample_id}_reference-{selection}_heatmap.png`

The final `infercnv` object includes the following slots:

- `infercnv_obj@ expr.data` : contains the processed expression matrix as it exists at the end of that stage for which that `inferCNV` object represents.

- `infercnv_obj@reference_grouped_cell_indices` : list containing the expression matrix column indices that correspond to each of the normal (reference) cell types.

- `infercnv_obj@observation_grouped_cell_indices` : similar list as above, but corresponds to the tumor cell types.

Based on the above slots, it would be straightforward to extract info of interest and/or move data into other analysis frameworks.

## Deprecated scripts

- `06b_build-normal_reference.R`
  - This script builds a Seurat object with normal cells passing a fixed score threshold for use by inferCNV.
  - The output is saved to `results/references/06b_normal-cell-reference.rds`.
