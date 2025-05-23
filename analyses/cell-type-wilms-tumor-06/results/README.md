# Azimuth compatible fetal references

To perform label transfer using code adapted from Azimuth, we prepare two references in [`scripts/prepare-fetal-references.R`](../scripts/prepare-fetal-references.R).
- First, we use the `fetal_full` object downloaded from CELLxGENE.
This is a fetal kidney reference from Stewart et al., and it is saved in `references/stewart_formatted_ref.rds`.
- Second, we format the Azimuth `fetusref` reference.
This is a fetal organ reference from Cao et al., and it is saved in `references/cao_formatted_ref.rds`.
- This script also downloads a file `references/homologs.rds` which the label transfer will use for converting Ensembl ids to gene names.

As part of the `00b_characterize_fetal_kidney_reference_Stewart.Rmd` notebook template, we characterized the fetal kidney reference and generated lists of marker genes for the compartment and cell types composing the reference.
The tables of marker genes can be found in the S3 bucket name `researcher-008971640512-us-east-2`:
- `references/00a_marker_cell-type_fetal_kidney_Stewart.csv`
- `references/00a_marker_compartment_fetal_kidney_Stewart.csv`


# Clustering and label transfer from fetal references

The S3 bucket name `researcher-008971640512-us-east-2` contains one folder for each of the samples of the Wilms tumor dataset SCPCP000006.

- `01-Seurat_{sample-id}.Rds` is the output of the [`01_seurat-processing.Rmd`](../notebook_template/01_seurat-processing.Rmd) notebook template.
In brief, the `_processed.rds` `sce object` is converted to `Seurat` and normalized using `SCTransform`.
Dimensionality reduction (`RunPCA` and `RunUMAP`) and clustering (`FindNeighbors` and `FindClusters`) are performed before saving the `Seurat` object.

- `02a-fetal_full_label-transfer_{sample-id}.Rds` is the output of the [`02a_label-transfer_fetal_full_reference_Cao.Rmd`](../notebook_template/02a_label-transfer_fetal_full_reference_Cao.Rmd) notebook template.
In brief, we used an `Azimuth`-adapted approach to transfer labels from the `Azimuth` fetal full reference (Cao et al.)

- `02b-fetal_kidney_label-transfer_{sample-id}.Rds` is the output of the [`02b_label-transfer_fetal_kidney_reference_Stewart.Rmd`](../notebook_template/02b_label-transfer_fetal_kidney_reference_Stewart.Rmd) notebook template.
In brief, we used an `Azimuth`-adapted approach to transfer labels from the fetal kidney reference (Stewart et al.)


Of note, as `01-Seurat_{sample-id}.Rds` is the input of `02a_label-transfer_fetal_full_reference_Cao.Rmd` and `02a-fetal_full_label-transfer_{sample-id}.Rds` the input of `02b-fetal_kidney_label-transfer_{sample-id}.Rds`, the latest `02b-fetal_kidney_label-transfer_{sample-id}.Rds` contains all metadata generated by the 3 notebook templates.

# Copy Number Alterations

## using `copykat`

We tried to infer aneuploidy in cancer cells using `copykat` and tested different parameters.
Mostly, we tested to run `copykat` with and without normal cells as reference and defined the clustering distance to be either `euclidean` or `spearman`.

We selected previously 5 samples to test for these parameters:
- sample `SCPCS000194`
- sample `SCPCS000179`
- sample `SCPCS000184`
- sample `SCPCS000205`
- sample `SCPCS000208`

For each of the samples, we ran `05_copyKAT.R` script, defining the `distance` parameters.

For each sample and each condition (reference and distance), we saved in `results/{sample_id}/05_copykat/{reference}/{distance}/`:

- the final `copykat` object in `05_copykat_{sample_id}_{reference}_distance-{selection}.rds`
- the heatmap of CNV in `05_copykat_{sample_id}_{reference}_distance-{selection}_copykat_heatmap.png`
- the prediction (aneuploid or diploid value per cell) in `05_copykat_{sample_id}_{reference}_distance-{selection}_copykat_prediction.txt`
- the CNA matrix `05_copykat_{sample_id}_{reference}_distance-{selection}_copykat_CNA_results.txt`

## using `infercnv`

We also tried to infer large CNV in cancer cells using `infercnv` and tested the sensibility of the output in regard to the definition of the normal cells.

We selected previously 5 samples to test for these parameters:
- sample `SCPCS000194`
- sample `SCPCS000179`
- sample `SCPCS000184`
- sample `SCPCS000205`
- sample `SCPCS000208`

`infercnv` requires a gene position file that we build in `06a_build-geneposition.R` and saved as `gencode_v19_gen_pos.complete.txt` in `results/references`.

For each sample and each condition (reference), we saved in `results/{sample_id}/06_infercnv/reference-{selection}`:
- the final `infercnv` object in `06_infercnv_{sample_id}_reference-{selection}.rds`
- the heatmap of CNV in `06_infercnv_{sample_id}_reference-{selection}_heatmap.png`

Of note, the final `infercnv` object includes the following slots:

- `infercnv_obj@expr.data` : contains the processed expression matrix as it exists at the end of that stage for which that `inferCNV` object represents.

- `infercnv_obj@reference_grouped_cell_indices` : list containing the expression matrix column indices that correspond to each of the normal (reference) cell types.

- `infercnv_obj@observation_grouped_cell_indices` : similar list as above, but corresponds to the tumor cell types.

Based on the above slots, it would be straightforward to extract info of interest and/or move data into other analysis frameworks.

In addition, for the condition `reference = "both"`, we ran `infercnv` with `HMM = TRUE`.
[HMM CNV prediction methods](https://github.com/broadinstitute/infercnv/wiki/inferCNV-HMM-based-CNV-Prediction-Methods) will allow us to explore the CNV results better, with an easy [merge](https://github.com/broadinstitute/infercnv/wiki/Extracting-features) of `infercnv` result with the `Seurat` object.
However, HMM CNV prediction methods uses a lot of resources, including time (~2h/sample/condition), and often causes the R session to crash.
This is why we only ran the HMM model for one `reference` condition. After selection of the best reference to use, we will run it for all samples.


## Annotations

Annotations in `SCPCP000006-annotations.tsv` were produced by `notebook/07_combined_annotation_across_samples_exploration.Rmd`.
