# Cell typing of Ewing sarcoma samples

This module will include code to annotate cell types in the Ewing sarcoma samples from SCPCP000015 present on the ScPCA Portal.

<!-- START doctoc generated TOC please keep comment here to allow auto update -->
<!-- DON'T EDIT THIS SECTION, INSTEAD RE-RUN doctoc TO UPDATE -->

- [Cell type annotation overview](#cell-type-annotation-overview)
  - [Annotation results](#annotation-results)
- [`AUCell` and `SingleR` annotation workflow](#aucell-and-singler-annotation-workflow)
  - [Usage](#usage)
  - [Input files](#input-files)
  - [Output files](#output-files)
  - [Computational resources](#computational-resources)
- [Clustering workflow](#clustering-workflow)
  - [Usage](#usage-1)
  - [Input files](#input-files-1)
  - [Output files](#output-files-1)
  - [Computational resources](#computational-resources-1)
- [Using `AUCell` to calculate gene set signatures](#using-aucell-to-calculate-gene-set-signatures)
  - [Usage](#usage-2)
  - [Input files](#input-files-2)
  - [Output files](#output-files-2)
  - [Computational resources](#computational-resources-2)
- [CNV annotation workflow](#cnv-annotation-workflow)
  - [Usage](#usage-3)
  - [Sample metadata](#sample-metadata)
  - [Input files](#input-files-3)
  - [Output files](#output-files-3)
  - [Annotation files](#annotation-files)
  - [Computational resources](#computational-resources-3)
- [Exploratory analyses](#exploratory-analyses)
- [Software requirements](#software-requirements)

<!-- END doctoc generated TOC please keep comment here to allow auto update -->

## Cell type annotation overview 

This module contains code and analysis used to annotate cell types in all Ewing sarcoma samples found on the ScPCA Portal in `SCPCP000015` that are part of the 2024-11-24 data release. 

The following steps were taken to assign cell types. 

Defining tumor cells: [`AUCell`](https://bioconductor.org/packages/release/bioc/html/AUCell.html) was run on the merged object containing all libraries with a set of MSigDB and custom gene sets specific to Ewing sarcoma using the `run-aucell-ews-signatures.sh` workflow.
See [below for detailed instructions on running this workflow](#using-aucell-to-calculate-gene-set-signatures). 

Tumor cells were classified as those that met the following criteria: 

- `tumor EWS-high`: 
  - AUC > 0.04 for [`aynaud-ews-targets`](./references/gene_signatures/aynaud-ews-targets.tsv)
  - AUC > 0.01 for [`STAEGE_EWING_FAMILY_TUMOR`](https://www.gsea-msigdb.org/gsea/msigdb/human/geneset/STAEGE_EWING_FAMILY_TUMOR.html)
- `tumor EWS-low`: 
  - AUC > 0.1 for [`wrenn-nt5e-genes`](./references/gene_signatures/wrenn-nt5e-genes.tsv)
  - AUC > 0.05 for [`HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION`](https://www.gsea-msigdb.org/gsea/msigdb/human/geneset/HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION.html)
- `tumor EWS-high proliferative`:
  - Cells that met the criteria for `tumor EWS-high` and had mean expression of marker genes for proliferative tumor cells (see [`references/tumor-cell-state-markers.tsv`](./references/tumor-cell-state-markers.tsv)) > 0  

Defining normal cell types: Consensus cell type labels were obtained by running the [`assign-consensus-celltypes.sh` workflow available in the `cell-type-consensus` module](../cell-type-consensus/assign-consensus-celltypes.sh). 
The consensus cell type label was used for all cells that did not meet the criteria for tumor cells mentioned above. 
If no consensus cell type is identified and the cell is not a tumor cell, it will have "Unknown" as the label. 

This script requires the processed `SingleCellExperiment` objects for all samples in `SCPCP000015` as input and was run using the following command: 

```sh
./assign-consensus-celltypes.sh "SCPCP000015"
```

Note that an alternative approach we used was to build a custom reference containing both tumor cells and normal cell types for use with `SingleR` to assign cell types. 
The results from this approachoch was concordant with the above outlined approach and is included in the TSV file with the final annotations (see the [results section below](#annotation-results)). 
See the [full instructions for the `SingleR` workflow below](#aucell-and-singler-annotation-workflow). 

### Annotation results 

See this [notebook](http://htmlpreview.github.io/?https://github.com/AlexsLemonade/OpenScPCA-analysis/blob/main/analyses/cell-type-ewings/exploratory_analysis/08-merged-celltypes.nb.html) for a summary of how the final cell type annotations were assigned and plots that validate these assignments. 

A TSV file containing the final annotations is stored in S3 at: 
`s3://researcher-211125375652-us-east-2/cell-type-ewings/results/final-annotations/SCPCP000015_celltype-annotations.tsv.gz`

The TSV file contains the following columns: 

| | | 
|--|--|
| `barcodes` | Unique cell barcode |
| `library_id` | ScPCA library id|
| `sample_id` | ScPCA sample id | 
| `sample_type` | Type of sample, either `patient tissue` or `patient-derived xenograft` |
| `singler_annotation` | The human readable term associated with the ontology term identifier for the associated cell type or `tumor-<library_id>`, output by `aucell-singer-annotation.sh` | 
| `singler_ontology` | The cell ontology identifier associated with the cell type or `tumor-<library_id>`, output by `aucell-singer-annotation.sh`  | 
| `consensus_annotation` | The human readable term associated with consensus cell type, output by `cell-type-consensus` | 
| `consensus_ontology` | The cell ontology identifier associated with the consensus cell type, output by `cell-type-consensus` | 
| `final_annotation` | The human readable term associated with the final annotation assigned using a combination of consensus cell types and `AUCell` results, output by `exploratory_analysis/08-merged-celltypes.Rmd` |
| `final_ontology` |The cell ontology identifier associated with the final annotation assigned using a combination of consensus cell types and `AUCell` results, output by `exploratory_analysis/08-merged-celltypes.Rmd`, if the `final_annotation` contains `tumor`, this column will match the `final_annotation` column instead of containing the ontology term | 

## `AUCell` and `SingleR` annotation workflow

The annotation workflow, `aucell-singler-annotation.sh`, can be used to identify potential tumor cells in a given sample using [`AUCell`](https://www.bioconductor.org/packages/release/bioc/html/AUCell.html) and [`SingleR`](https://bioconductor.org/packages/release/bioc/html/SingleR.html).
This workflow includes the following steps:

1. An AUC threshold is determined from running `AUCell` on the reference sample `SCPCS000490`.
   Tumor cells are classified as those that have an AUC score greater than or equal to the AUC threshold determined by `AUCell`.
   `SCPCS000490` was previously determined to have a bimodal distribution of AUC scores, where cells with AUC scores higher than the determined threshold corresponded to tumor cells classified using CNV methods in `cnv-annotation.sh`.
   This AUC threshold is then returned and used for classifying the remaining samples.
2. `AUCell` is then run on all samples and tumor cells are classified as those with an AUC score greater than or equal to the AUC threshold determined in the reference sample.
3. Gene set scores are calculated by taking the mean gene expression for all genes in the following gene sets:

- [`ZHANG_TARGETS_OF_EWSR1_FLI1_FUSION`](https://www.gsea-msigdb.org/gsea/msigdb/human/geneset/ZHANG_TARGETS_OF_EWSR1_FLI1_FUSION.html)
- [`RIGGI_EWING_SARCOMA_PROGENITOR_UP`](https://www.gsea-msigdb.org/gsea/msigdb/human/geneset/RIGGI_EWING_SARCOMA_PROGENITOR_UP.html?ex=1)
- [`SILIGAN_TARGETS_OF_EWS_FLI1_FUSION_DN`](https://www.gsea-msigdb.org/gsea/msigdb/cards/SILIGAN_TARGETS_OF_EWS_FLI1_FUSION_DN)

4. A summary report is generated that contains the results from `AUCell` and validation of tumor cells by looking at marker gene expression and gene set scores.
5. All cells classified as tumor by `AUCell` are combined to create a single tumor cell reference to use with `SingleR`.
6. `SingleR` is run on all samples using the combined tumor cell reference and two normal references, `BlueprintEncodeData` and `HumanPrimaryCellAtlasData` from [`celldex`](https://bioconductor.org/packages/devel/data/experiment/vignettes/celldex/inst/doc/userguide.html).
7. A summary report is generated that contains the results from `SingleR`, including a comparison to the results from running `AUCell`.

### Usage

The `aucell-singler-annotation.sh` workflow can be used to annotate tumor and normal cells in all Ewing sarcoma samples from SCPCP000015.
The workflow can be run using the following command:

```sh
./aucell-singler-annotation.sh
```

### Input files

The `aucell-singler-annotation.sh` workflow requires the processed `SingleCellExperiment` objects (`_processed.rds`) from SCPCP0000015.
These files were obtained using the `download-data.py` script:

```sh
# download SCE objects
./download-data.py --projects SCPCP000015
```

The workflow also requires the tumor marker gene reference file, [`references/tumor-marker-genes.tsv`](./references/tumor-marker-genes.tsv).
This file contains a list of marker genes for identifying Ewing sarcoma tumor cells.
The all marker gene reference file, [`references/all-visser-marker-genes.tsv`](./references/visser-all-marker-genes.tsv) is also required and contains a list of marker genes for both tumor and normal cell types obtained from [Visser, Beligs, _et al._ (2023)](https://doi.org/10.1158/2767-9764.CRC-23-0027).

### Output files

Running the `aucell-singler-annotation.sh` workflow will generate the following output files in `results/aucell_singler_annotation` for each sample/library combination.

```
aucell_annotation
├── <sample_id>
    ├── <library_id>_auc-classifications.tsv
    ├── <library_id>_aucell-report.html
    ├── <library_id>_gene-set-scores.tsv
    ├── <library_id>_marker-gene-classifications.tsv
    ├── <library_id>_singler-classifications.tsv
    ├── <library_id>_singler-report.html
```

The `aucell-report.html` file is a rendered report exploring the results from `AUCell`.
The `singler-report.html` file is a rendered report exploring the results from `SingleR`.
Both reports include validation of tumor cells by looking at marker gene expression and gene set scores.

The `auc-classifications.tsv` file contains the following columns:

|                      |                                                                                                                                                      |
| -------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------- |
| `barcodes`           | Unique cell barcode                                                                                                                                  |
| `auc`                | AUC score reported by running `AUCell`                                                                                                               |
| `auc_classification` | Either `Tumor` or `Normal`, where `Tumor` cells are those that have an AUC greater than or equal to the AUC threshold found in the reference sample. |

The `marker-gene-classifications.tsv` file contains the following columns:

|                              |                                                                                                  |
| ---------------------------- | ------------------------------------------------------------------------------------------------ |
| `barcodes`                   | Unique cell barcode                                                                              |
| `marker_gene_classification` | Either `Tumor` or `Normal`, where `Tumor` cells are those that express at least one marker gene. |

The `singler-classifications.tsv` file contains the following columns:

|                      |                                                                                                                           |
| -------------------- | ------------------------------------------------------------------------------------------------------------------------- |
| `barcodes`           | Unique cell barcode                                                                                                       |
| `singler_annotation` | The human readable term associated with the ontology term identifier for the associated cell type or `tumor-<library_id>` |
| `singler_ontology`   | The cell ontology identifier associated with the cell type or `tumor-<library_id>`                                        |
| `aucell_annotation`  | Either `tumor` or `normal` as determined by running `AUCell` in `aucell-annotation.sh`                                    |

The `gene-set-scores.tsv` file contains the scores (mean and sum) for all genes in three different `EWS-FLI1` target gene sets.

### Computational resources

The `aucell-singler-annotation.sh` uses 4 CPUs and can be run locally on a laptop or on a virtual computer on Lightsail for Research.
Note that `SingleR` does take up to 1 hour to run with 4 CPUs for some samples in this project.

## Clustering workflow 

The clustering workflow (`evaluate-clusters.sh`) can be used to perform and evaluate clusters across a range of parameters. 
For all processed `SingleCellExperiment` objects that are part of `SCPCP000015`, clustering is performed using the following options: 

- Algorithms: Louvain, Leiden with CPM, and Leiden with modularity, all using Jaccard for the weighting scheme
- Nearest neighbors: 5, 10, 15, 20, 25, 30, 35, 40
- Resolution (Louvain and Leiden with modularity): 0.5, 1, 1.5 
- Resolution (Leiden with CPM): 0.001, 0.005, 0.01

After clusters are calculated, a report with summary plots showing silhouette width, cluster purity, and cluster stability for each set of parameters used to calculate clusters is created. 

Clusters identified using the Leiden algorithm with the modularity objective function are then explored in a secondary report. 
A report is generated for resolution 0.5 and resolution 1.
This report includes: 

- A comparison of cluster assignments to cell type assignments (output from the `aucell-singler-annotation.sh` workflow) for all values of nearest neighbors used. 
- Marker gene set expression in each cluster for all marker gene sets in `references/visser-all-marker-genes.tsv`


### Usage

The `evaluate-clusters.sh` workflow can be used to perform clustering of all Ewing sarcoma samples from SCPCP000015.
The workflow can be run using the following command:

```sh
./evaluate-clusters.sh
```

Calculating the clusters and the metrics can be time consuming and takes ~ 24 hours on a laptop for all samples in the project. 
To skip this step and only render the second report to look at cell types and marker gene expression, use the `skip_metrics=1` option: 

```sh
skip_metrics=1 ./evaluate-clusters.sh
```

### Input files 

The `evaluate-clusters.sh` workflow requires the processed `SingleCellExperiment` objects (`_processed.rds`) from SCPCP000015.
These files were obtained using the `download-data.py` script:

```sh
# download SCE objects
./download-data.py --projects SCPCP000015
```

### Output files 

Running the `evaluate-clusters.sh` workflow will generate the following output files in `results/clustering` for each sample/library combination.

```
clustering
├── <sample_id>
    ├── <library_id>_cluster-results.tsv
    ├── <library_id>_cluster-summary-report.html
```

The `cluster-summary-report.html` file is a rendered report that compares clustering results from all clustering parameters tested. 
In particular, plots showing silhouette width, cluster purity, and cluster stability are shown.

The `cluster-results.tsv` file contains the following columns:

|                      |                                                                                                                                                      |
| -------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------- |
| `cell_id`           | Unique cell barcode                                                                                                                                  |
| `cluster_method`                | Algorithm/objective function combination used for clustering, one of `louvain`, `leiden_cpm` or `leiden_mod`                                                                                                            |
| `cluster` | Assigned cluster |
| `algorithm` | Algorithm used for clustering, one of `louvain` or `leiden` |
| `weighting` | Weighting scheme used for clustering |
| `nn` | Number of nearest neighbors used for clustering |
| `resolution` | Resolution used for clustering |
| `objective_function` | Objective function used for clustering, only applicable for leiden clustering |

### Computational resources 

The `evaluate-clusters.sh` uses 4 CPUs and can be run locally on a laptop or on a virtual computer on Lightsail for Research.


## Using `AUCell` to calculate gene set signatures 

The script, `run-aucell-ews-signatures.sh`, is used to run `AUCell` with a set of custom gene signatures on all samples in SCPCP000015. 
By default, AUC values are calculated for all gene signatures in [references/gene_signatures](references/gene_signatures/) and a set of `MSigDB` signatures associated with high and low EWS-FLI1 expression. 
The full list of gene signatures used can be found in [the references `README.md`](references/README.md#gene-signatures). 

`AUCell` is run for each gene signature, and AUC values along with the AUC threshold reported by `AUCell` are saved to a TSV file. 
By default, `AUCell` is run with an `aucMaxRank` value equal to 1% of the detected genes in the processed object. 
Results from `AUCell` will be generated for each library and for the merged object containing all samples in `SCPCP000015`. 

### Usage

The `run-aucell-ews-signatures.sh` script can be run using the following command:

```sh
./run-aucell-ews-gene-signatures.sh
```

To a different percentage of detected genes to determine `aucMaxRank` use the following command, specifying a max_rank_threshold between 0-1: 

```sh
max_rank_threshold=.05 ./run-aucell-ews-signatures.sh
```

### Input files

The `run-aucell-ews-signatures.sh` workflow requires the processed `SingleCellExperiment` objects (`_processed.rds`) from SCPCP0000015 and the merged object for SCPCP000015.
These files were obtained using the following commands:

```sh
# download SCE objects
./download-data.py --projects SCPCP000015

# download merged object
./download-results.py --projects SCPCP000015 --module merge-sce
```

The workflow also requires the custom marker gene sets present in [`references/gene_signatures`](./references/gene_signatures/). 

### Output files

Running the `run-aucell-ews-signatures.sh` script will generate the following output files in `results/aucell-ews-signatures` for each sample/library combination.

```
aucell-ews-signatures
├── <sample_id>
    ├── <library_id>_auc-ews-gene-signatures.tsv
```

The `auc-ews-gene-signatures.tsv` file contains the following columns:

|                      |                                                                                                                                                      |
| -------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------- |
| `barcodes`           | Unique cell barcode                                                                                                                                  |
| `gene_set`           | The name associated with the gene set used for calculating AUC values                                                                                |
| `auc`                | AUC value reported by running `AUCell`                                                                                                               |
| `auc_threshold` | The threshold AUC value used to classify cells as having high expression of genes in the gene signature as reported by `AUCell` |


### Computational resources 

`run-aucell-ews-signatures.sh` uses 4 CPUs and can be run locally on a laptop or on a virtual computer on Lightsail for Research.

## CNV annotation workflow

**NOTE:** This workflow is no longer used in the cell type annotation analysis, has been removed from CI, and is no longer maintained! 

The CNV annotation workflow (`cnv-annotation.sh`) can be used to identify potential tumor cells in a given sample.
Annotations are obtained by running the following methods within the workflow:

- Identify cells that express high levels of tumor marker genes.
  The full list of marker genes can be found in `references/tumor-marker-genes.tsv`.
- Use the list of tumor marker genes to classify tumor cells with [`CellAssign`](https://docs.scvi-tools.org/en/stable/user_guide/models/cellassign.html).
- Identify copy number variations and annotate tumor cells using [`CopyKAT`](https://github.com/navinlabcode/copykat).
- Identify copy number variations using [`InferCNV`](https://github.com/broadinstitute/inferCNV/wiki).
  This returns a proportion of each chromosome with a CNV detected.
  We then calculate the genomic CNV proportion for each cell across all chromosomes, weighted by the number of genes in a chromosome.
  Cells with a genomic CNV proportion greater than the mean CNV proportion across all cells are called as tumor cells.

### Usage

The `cnv-annotation.sh` workflow can be used to annotate tumor and normal cells in the Ewing's sarcoma samples from SCPCP000015.
**Note:** Before running this workflow be sure to run `renv::restore()` and activate the conda environment using the following commands:

```sh
Rscript -e "renv::restore()"
conda activate openscpca-cell-type-ewings
```

The following arguments are optional and can be used to run this workflow on additional samples (default sample is `SCPCS000490`):

- `sample_id`: Unique sample ID (name of folder containing library data)
- `normal_celltypes`: Comma separated list of cell types annotated with either `SingleR` or `CellAssign` to use as a reference list of normal cells.
  This should correspond to the value found in `sample-metadata.tsv` for this sample.
- `tumor_celltypes`: Comma separated list of cell types annotated with either `SingleR` or `CellAssign` that are expected to align with tumor cells.
  Any cell types used here will be used for comparing to tumor cells annotated in this workflow.
  This should correspond to the value found in `sample-metadata.tsv` for this sample.

**Note:** Before running the workflow, make sure that [`sample-metadata.tsv`](#sample-metadata) contains the information for any samples being processed.

Example of running the workflow with a different sample:

```sh
sample_id="SCPCS000491" ./cnv-annotation.sh
```

### Sample metadata

The `sample-metadata.tsv` file is a TSV file containing information about samples that have been run through the `cnv-annotation.sh` workflow.
This is manually updated as new samples are used in analysis.

Each row represents a unique sample and library combination.
The following columns are present in this file:

|                    |                                                                                                                                        |
| ------------------ | -------------------------------------------------------------------------------------------------------------------------------------- |
| `scpca_sample_id`  | Unique sample ID. The `sample_id` corresponds to the folder name containing data files for that sample after using `download-data.py`. |
| `scpca_library_id` | Unique library ID. The `library_id` will match the prefix of all data files (`.rds` and `.h5ad`) downloaded using `download-data.py`.  |
| `normal_celltypes` | A comma separated list of cell types annotated with either `SingleR` or `CellAssign` used to create a reference list of normal cells   |
| `tumor_celltypes`  | A comma separated list of cell types annotated with either `SingleR` or `CellAssign` that are expected to align with tumor cells.      |

**Note:** To identify the cell type annotations to use for `normal_celltypes` and `tumor_celltypes`, reference the plots found in `<library_id>_celltype-report.html`.
These can be downloaded using the `--include_reports` option in `download-data.py`.

### Input files

The `cnv-annotation.sh` workflow requires the processed `SingleCellExperiment` objects (`_processed.rds`) and processed `AnnData` objects (`.hdf5` files) from SCPCP0000015.
These files were obtained using the `download-data.py` script:

```sh
# download both SCE and AnnData objects
./download-data.py --projects SCPCP000015 --format "SCE,AnnData"
```

The workflow also requires the following reference files:

- `references/tumor-marker-genes.tsv`: This file contains a list of marker genes for identifying Ewing sarcoma tumor cells.
- `references/all-marker-genes.tsv`: This file contains a list of marker genes for identifying both normal and Ewing sarcoma tumor cells as described in [Visser, Beligs, _et al._ (2023)](https://doi.org/10.1158/2767-9764.CRC-23-0027).

### Output files

Running the `cnv-annotation.sh` workflow will generate the following output files in `results/cnv_annotation` for each sample/library combination.

```
cnv_annotation
├── <sample_id>
    ├── <library_id>_cellassign-classifications.tsv
    ├── <library_id>_cellassign-report.html
    ├── <library_id>_copykat-classifications.tsv
    ├── <library_id>_copykat-report.html
    ├── <library_id>_infercnv-classifications.tsv
    ├── <library_id>_infercnv-report.html
    ├── <library_id>_marker-gene-report.html
    ├── <library_id>_tumor-normal-classifications.tsv
    ├── annotations
    │   └── <library_id>_reference-cells.tsv
    ├── cellassign
    │   ├── <library_id>_panglao-endo-fibro_predictions.tsv
    │   ├── <library_id>_tumor-marker_predictions.tsv
    │   └── <library_id>_visser-all-marker_predictions.tsv
    ├── copykat
    │   ├── no_reference
    │   │   ├── <library_id>_copykat_heatmap.jpeg
    │   │   └── <library_id>_final-copykat.rds
    │   └── with_reference
    │       ├── <library_id>_copykat_heatmap.jpeg
    │       └── <library_id>_final-copykat.rds
    └── infercnv
        ├── <library_id>_cnv-metadata.tsv
        ├── <library_id>_cnv-obj.rds
        └── <library_id>_infercnv.png
```

All `.html` files are rendered reports summarizing use of each method (indicated in the filename) to classify tumor cells.
All `classifications.tsv` files contain the final annotation as reported by each method.
The `annotations` folder contains the reference table indicating which cells were used as "normal" or "tumor" cells in various analysis.
See [below](#annotation-files) for more information on this table.

### Annotation files

Additionally, a reference TSV will be generated for each library containing a table of cell types from `SingleR` and `CellAssign`.
This file contains any cell barcodes expected to line up with tumor or normal cell types specified with the `normal_celltypes` and `tumor_celltypes` arguments when running `cnv-annotation.sh`.
This table will be saved in `<sample_id>/annotations/<library_id>_reference-cells.tsv` and contains the following columns:

|                                  |                                                                                                |
| -------------------------------- | ---------------------------------------------------------------------------------------------- |
| `barcodes`                       | Cell barcode                                                                                   |
| `reference_cell_class`           | Indicates if the cell should be uses as a Normal or Tumor cell reference                       |
| `cellassign_celltype_annotation` | Original annotation as obtained by `CellAssign` in the processed `SingleCellExperiment` object |
| `singler_celltype_annotation`    | Original annotation as obtained by `SingleR` in the processed `SingleCellExperiment` object    |

### Computational resources

Currently, the `cnv-annotation.sh` can only be run using 1 CPU if running locally on a Mac.
Use `threads=1` when running the workflow to specify this.

To increase the speed of the workflow, we recommend running on a computer with at least 12 CPUs.
If this is the case, you may run the workflow using the default settings of 4 threads.

## Exploratory analyses

Any notebooks or scripts used specifically for exploratory analyses can be found in the `exploratory_analysis` folder.
The `exploratory_analysis/annotation_notebooks` folder contains notebooks used to generate and save cell type annotations for samples in `SCPCP000015`.
For more information on how these were generated, see the [README.md](./exploratory_notebooks/annotation_notebooks/README.md).

## Software requirements

All notebooks and R scripts in this module should be run with the `ewings-cell-type.Rproj`.
The `renv.lock` file contains all package and version information.

All Python scripts must be run using the conda environment `openscpca-cell-type-ewings`.
To create and activate this environment from the `.yml` file use:

```sh
conda env create -f environment.yml --name openscpca-cell-type-ewings
conda activate openscpca-cell-type-ewings
```
