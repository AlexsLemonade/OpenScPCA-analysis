# Notebook directory instructions

The notebook directory holds a directory for each of the sample in the Wilms tumor dataset `SCPCP000006` and the fetal kidney reference that we used for label transfer.

## Fetal kidney reference

To perform label transfer using an Azimuth-adapted approach and the fetal kidney atlas, a reference is built via [`scripts/prepare-fetal-references.R`](../scripts/scripts/prepare-fetal-references.R).

As part of the `00b_characterize_fetal_kidney_reference_Stewart.Rmd` notebook, we characterized the fetal kidney reference and generated lists of marker genes for the compartment and cell types composing the reference.


## Analysis per sample

Each of the sample of the Wilms tumor dataset SCPCP000006 as been pre-processed and characterized as the following.
Reports for each of the steps are found in the notebook/{sample_id} directory:

- [x] `01_seurat_processing_{sample-id}.html` is the output of the [`01_seurat-processing.Rmd`](../notebook_template/01_seurat-processing.Rmd) notebook template.
In brief, the `_processed.rds` `sce object` is converted to `Seurat` and normalized using `SCTransform`.
Dimensionality reduction (`RunPCA` and `RunUMAP`) and clustering (`FindNeighbors` and `FindClusters`) are performed before saving the `Seurat` object.

- [x] `02a_fetal_full_label-transfer_{sample-id}.html` is the output of the [`02a_label-transfer_fetal_full_reference_Cao.Rmd`](../notebook_template/02a_label-transfer_fetal_full_reference_Cao.Rmd) notebook template.
In brief, we used an Azimuth-adapted approach to transfer labels from the Azimuth fetal full reference (Cao et al.)

- [x] `02b_fetal_kidney_label-transfer_{sample-id}.html` is the output of the [`02b_label-transfer_fetal_kidney_reference_Stewart.Rmd`](../notebook_template/02b_label-transfer_fetal_kidney_reference_Stewart.Rmd) notebook template.
In brief, we used an Azimuth-adapted approach to transfer labels from the fetal kidney reference (Stewart et al.)

- [x] `03_clustering_exploration_{sample-id}.html` is the output of the [`03_clustering_exploration.Rmd`](../notebook_template/03_clustering_exploration.Rmd) notebook template.
In brief, we explore the clustering results, we look into some marker genes, pathways enrichment and label transfer.


## Global analysis

The next step in analysis is to identify tumor vs. normal cells.

- [x] `04_annotation_Across_Samples_exploration.html` is the output of the [`04_annotation_Across_Samples_exploration.Rmd`](../notebook/04_annotation_Across_Samples_exploration.Rmd) notebook.
In brief, we explored the label transfer results across all samples in the Wilms tumor dataset `SCPCP000006` in order to identify a few samples that we can begin next analysis steps with.

One way to evaluate the label transfer is to look at the `predicted.score` for each label being transferred, which more or less correspond to the certainty for a label transfer to be `TRUE`. More information on the cell-level metric `predicted.score` can be found in the [mapping QC](https://azimuth.hubmapconsortium.org/#Mapping%20QC) section of `Azimuth` documentation.

We render the notebook with different thresholds for the `predicted.score` and evaluate the impact of filtering out cells with a `predicted.score` below 0.5, 0.75, 0.85 and 0.95.

Of important notes:

- The stroma compartment often has a poor `predicted.score`. This is for me an indication that these cells might be cancer cells and not normal stromal cells.

- We would rather use the `predicted.score` threshold to select normal cells for which we have a high confidence, i.e. immune and endothelial cells, but not to filter out all cells below the threshold.

- While a `predicted.score` of 0.5 is much too low (almost all cells having a higher `predicted.score`) and 0.95 is too high (so few cells pass the threshold), 0.75 and 0.85 looked both appropriate for our purpose.

--> We decided to go with the most stringent threshold of 0.85 as we want to be sure of our selection of normal cells (i.e. endothelial and immune cells) that we will use to run `inferCNV`.

- [x] `07_combined_annotation_across_samples_exploration.html` is the output of the [`07_combined_annotation_across_samples_exploration.Rmd`](../notebook/07_combined_annotation_across_samples_exploration.Rmd) notebook.
This notebook performs a draft annotation of samples using information from CNV inference and label transfer.
