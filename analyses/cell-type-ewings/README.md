# Cell typing of Ewing sarcoma samples

This module will include code to annotate cell types in the Ewing sarcoma samples from SCPCP000015 present on the ScPCA Portal.

<!-- START doctoc generated TOC please keep comment here to allow auto update -->
<!-- DON'T EDIT THIS SECTION, INSTEAD RE-RUN doctoc TO UPDATE -->
**Table of Contents**

- [`AUCell` annotation workflow](#aucell-annotation-workflow)
  - [Usage](#usage)
  - [Input files](#input-files)
  - [Output files](#output-files)
  - [Computational resources](#computational-resources)
- [CNV annotation workflow](#cnv-annotation-workflow)
  - [Usage](#usage-1)
  - [Sample metadata](#sample-metadata)
  - [Input files](#input-files-1)
  - [Output files](#output-files-1)
  - [Annotation files](#annotation-files)
  - [Computational resources](#computational-resources-1)
- [Exploratory analyses](#exploratory-analyses)
- [Software requirements](#software-requirements)

<!-- END doctoc generated TOC please keep comment here to allow auto update -->

## `AUCell` annotation workflow

The `AUCell` annotation workflow (`aucell-annotation.sh`) can be used to identify potential tumor cells in a given sample using [`AUCell`](https://www.bioconductor.org/packages/release/bioc/html/AUCell.html).
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


### Usage

The `aucell-annotation.sh` workflow can be used to annotate tumor and normal cells in all Ewing sarcoma samples from SCPCP000015.
The workflow can be run using the following command:

```sh
./aucell-annotation.sh
```

### Input files

The `aucell-annotation.sh` workflow requires the processed `SingleCellExperiment` objects (`_processed.rds`) from SCPCP0000015.
These files were obtained using the `download-data.py` script:

```sh
# download SCE objects
./download-data.py --projects SCPCP000015
```
The workflow also requires the tumor marker gene reference file, [`references/tumor-marker-genes.tsv`](./references/tumor-marker-genes.tsv).
This file contains a list of marker genes for identifying Ewing sarcoma tumor cells.

### Output files

Running the `aucell-annotation.sh` workflow will generate the following output files in `results/aucell_annotation` for each sample/library combination.

```
aucell_annotation
├── <sample_id>
    ├── <library_id>_auc-classifications.tsv
    ├── <library_id>_aucell-report.html
    ├── <library_id>_gene-set-scores.tsv
    ├── <library_id>_marker-gene-classification.tsv
```
The `.html` file is a rendered report exploring the results from `AUCell` and validating tumor cells by looking at marker gene expression and gene set scores.
The `auc-classifications.tsv` file contains the following columns:

| | |
|----|----|
| `barcodes` | Unique cell barcode |
| `auc` | AUC score reported by running `AUCell` |
| `auc_classification` | Either `Tumor` or `Normal`, where `Tumor` cells are those that have an AUC greater than or equal to the AUC threshold found in the reference sample.  |

The `marker-gene-classifications.tsv` file contains the following columns:

| | |
|----|----|
| `barcodes` | Unique cell barcode |
| `marker_gene_classification` | Either `Tumor` or `Normal`, where `Tumor` cells are those that express at least one marker gene.  |

The `gene-set-scores.tsv` file contains the scores (mean and sum) for all genes in three different EWS-FLI1 target gene sets.

### Computational resources

The `aucell-annotation.sh` uses 4 CPUs and can be run locally on a laptop.

## CNV annotation workflow

The CNV annotation workflow (`cnv-annotation.sh`) can be used to identify potential tumor cells in a given sample.
Annotations are obtained by running the following methods within the workflow:

- Identify cells that express high levels of tumor marker genes.
The full list of marker genes can be found in `references/tumor-marker-genes.tsv`.
- Use the list of tumor marker genes to classify tumor cells with [`CellAssign`](https://docs.scvi-tools.org/en/stable/user_guide/models/cellassign.html).
- Identify copy number variations and annotate tumor cells using [`CopyKAT`](https://github.com/navinlabcode/copykat).
- Identify copy number variations using [`InferCNV`](https://github.com/broadinstitute/inferCNV/wiki).
This returns a proportion of each chromosome with a CNV detected.
We then calculate the genomic CNV proportion for each cell across all chromosomes, weighted by the number of genes in a chromosome.
Cells with a genomic CNV proportion greater than the mean  CNV proportion across all cells are called as tumor cells.

### Usage

The `cnv-annotation.sh` workflow can be used to annotate tumor and normal cells in the Ewing's sarcoma samples from SCPCP000015.
**Note:** Before running this workflow be sure to run `renv::restore()`  and activate the conda environment using the following commands:

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

|  |   |
| --- | ---- |
|`scpca_sample_id`| Unique sample ID. The `sample_id` corresponds to the folder name containing data files for that sample after using `download-data.py`. |
|`scpca_library_id` | Unique library ID. The `library_id` will match the prefix of all data files (`.rds` and `.h5ad`) downloaded using `download-data.py`. |
|`normal_celltypes`| A comma separated list of cell types annotated with either `SingleR` or `CellAssign` used to create a reference list of normal cells |
|`tumor_celltypes`| A comma separated list of cell types annotated with either `SingleR` or `CellAssign` that are expected to align with tumor cells. |

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

|  |   |
| -------- | -------- |
| `barcodes` | Cell barcode |
| `reference_cell_class` | Indicates if the cell should be uses as a Normal or Tumor cell reference |
| `cellassign_celltype_annotation` | Original annotation as obtained by `CellAssign` in the processed `SingleCellExperiment` object |
| `singler_celltype_annotation` | Original annotation as obtained by `SingleR` in the processed `SingleCellExperiment` object |

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
