# `scripts` directory instructions

This directory contains all scripts used for cell typing Wilms tumor samples from SCPCP000006.


## `download-reference-files.R`

This script is used to download the fetal kidney reference (Stewart et al) and the Azimuth homologs file for converting gene IDs.

## `prepare-fetal-references.R`

This script is used to create the fetal kidney reference (Stewart et al) and the full fetal reference (Cao et al).

## `explore-cnv-methods.R`

This script is used to run 2 CNV inference methods `copykat` and `infercnv` for a selection of samples from SCPCP000006.

We selected previously 5 samples to test for different parameters of `copykat` and `infercnv`:

- sample SCPCS000194 has > 85 % of cells predicted as kidney and 234 + 83 endothelium and immune cells.
- sample SCPCS000179 has > 94 % of cells predicted as kidney and 25 + 111 endothelium and immune cells.
- sample SCPCS000184 has > 96 % of cells predicted as kidney and 39 + 70 endothelium and immune cells.
- sample SCPCS000205 has > 89 % of cells predicted as kidney and 92 + 76 endothelium and immune cells.
- sample SCPCS000208 has > 95 % of cells predicted as kidney and 18 + 35 endothelium and immune cells.

### `06a_build-geneposition.R`

This script build the gene position file that will be used in `06_infercnv.R` for each of the samples.
We build the position file once via downloading the `gencode_v19_gen_pos.complete.txt` from the [Trinity/CTAT](https://data.broadinstitute.org/Trinity/CTAT/cnv/gencode_v19_gen_pos.complete.txt)
05_copyKAT.R and 06_infercnv.R for 5 samples.

The `gencode_v19_gen_pos.complete.txt` is then saved in `results/references`.

### `05_copyKAT.R`

This script is used to run `copykat`.
Please note that `copykat` can take more than one hour per sample, even when using quite a lot of resources.
By default, we run `copykat` with 16 cores but we are actually calling it with 32 cores (using 32 cpus).

For each sample and each condition (reference and distance), we saved in `results/{sample_id}/05_copykat/{reference}/{distance}/`:

- the final `copykat` rds object in `05_copykat_{sample_id}_{reference}_distance-{selection}.rds`
- the heatmap of CNV in `05_copykat_{sample_id}_{reference}_distance-{selection}_copykat_heatmap.jpeg`
- the prediction (aneuploid or diploid value per cell) in `05_copykat_{sample_id}_{reference}_distance-{selection}_copykat_prediction.txt`
- the CNA matrix `05_copykat_{sample_id}_{reference}_distance-{selection}_copykat_CNA_results.txt`

### `06_infercnv.R`

This script is used to run `infercnv`.
We test for the sensitivity of `infercnv` in regards to the selection of normal cells used as reference.
Therefore, the parameter `reference` indicates if we want to select immune and/or endothelial cells as reference or no cell at all.

For each sample and each condition (reference), we saved in `results/{sample_id}/06_infercnv/reference-{selection}`:
- the final `infercnv`rds object in `06_infercnv_{sample_id}_reference-{selection}.rds`
- the heatmap of CNV in `06_infercnv_{sample_id}_reference-{selection}_heatmap.png`

The final `infercnv` rds object includes the following slots:

- 'infercnv_obj@ expr.data' : contains the processed expression matrix as it exists at the end of that stage for which that inferCNV object represents.

- 'infercnv_obj@reference_grouped_cell_indices' : list containing the expression matrix column indices that correspond to each of the normal (reference) cell types.

- 'infercnv_obj@observation_grouped_cell_indices' : similar list as above, but corresponds to the tumor cell types.

Based on the above slots, it would be straightforward to extract info of interest and/or move data into other analysis frameworks.
