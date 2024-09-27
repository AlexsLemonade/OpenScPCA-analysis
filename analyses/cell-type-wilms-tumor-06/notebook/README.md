# Notebook directory instructions

The notebook directory holds subdirectory for each of the sample in the Wilms tumor dataset SCPCP000006 and the fetal kidney reference that we used for label transfer. 

## Azimuth compatible fetal kidney reference

To perform label transfer using Azimuth and the fetal kidney atlas, a reference is built via [`scripts/download-and-create-fetal-kidney-ref.R`](../scripts/download-and-create-fetal-kidney-ref.R) using the fetal_full.Rds object download from:
"https://datasets.cellxgene.cziscience.com/40ebb8e4-1a25-4a33-b8ff-02d1156e4e9b.rds"

As part of the `00b_characterize_fetal_kidney_reference_Stewart.Rmd` notebook template, we characterized the fetal kidney reference and generated lists of marker genes for the compartment and cell types composing the reference. 


## Analysis per sample

Each of the sample of the Wilms tumor dataset SCPCP000006 as been pre-processed and characterized as the following. 
Reports for each of the steps are found in the notebook/{sample_id} directory:

- `01_seurat_processing_{sample-id}.html` is the output of the [`01_seurat-processing.Rmd`](../notebook_template/01_seurat-processing.Rmd) notebook template.
In brief, the `_processed.rds` `sce object` is converted to `Seurat` and normalized using `SCTransform`.
Dimensionality reduction (`RunPCA` and `RunUMAP`) and clustering (`FindNeighbors` and `FindClusters`) are performed before saving the `Seurat` object. 

- `02a_fetal_full_label-transfer_{sample-id}.html` is the output of the [`02a_label-transfer_fetal_full_reference_Cao.Rmd`](../notebook_template/02a_label-transfer_fetal_full_reference_Cao.Rmd) notebook template.
In brief, we used `Azimuth` to transfer labels from the `Azimuth` fetal full reference (Cao et al.) 

- `02b_fetal_kidney_label-transfer_{sample-id}.html` is the output of the [`02b_label-transfer_fetal_kidney_reference_Stewart.Rmd`](../notebook_template/02b_label-transfer_fetal_kidney_reference_Stewart.Rmd) notebook template.
In brief, we used `Azimuth` to transfer labels from the fetal kidney reference (Stewart et al.) 

- `03_clustering_exploration_{sample-id}.html` is the output of the [`03_clustering_exploration.Rmd`](../notebook_template/03_clustering_exploration.Rmd) notebook template.
In brief, we explore the clustering results, we look into some marker genes, pathways enrichment and label transfer.


## Global analysis

At some step during the analysis, a global picture was required to interprete correctly the results and take the best decision for the next steps. 

- `04_annotation_Across_Samples_exploration.html` is the output of the [`04_annotation_Across_Samples_exploration.Rmd`](../notebook/04_annotation_Across_Samples_exploration.Rmd) notebook. 
In brief, we explored the label transfer results across all samples in the Wilms tumor dataset SCPCP000006.