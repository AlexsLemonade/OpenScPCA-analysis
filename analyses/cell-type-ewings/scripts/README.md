# Scripts

This directory contains all scripts used for cell typing Ewing sarcoma samples from SCPCP000015.

## Scripts used in the CNV annotation workflow

The scripts in the `cnv-workflow` folder are implemented in `cnv-annotation.sh`:

1. `00-generate-cellassign-refs.R`: This script is used to generate the binary marker gene reference matrices required to run `CellAssign` in `cnv-annotation.sh`.
Running this script creates three different references found in `references/cellassign_refs`.
See [references/README.md](../references/README.md#cellassign-references) for a full description of all references that are created.

To run this script use the following command:

```sh
Rscript 00-generate-cellassign-refs.R
```

2. `00-make-gene-order-file.R`: This script is used to generate the [gene order file](https://github.com/broadinstitute/inferCNV/wiki/File-Definitions#gene-ordering-file) needed for running `InferCNV` with `04-run-infercnv.R`.
This script downloads the GTF file used to create the original index used in `scpca-nf` (Ensembl v104) from the public bucket `s3://scpca-references`.
This file is then converted to the required gene order formatted file for `InferCNV` and saved to `references/infercnv_refs`.

To run this script use the following command:

```sh
Rscript 00-make-gene-order-file.R
```

3. `01-select-cell-types.R`: This script is used to generate a table of normal and tumor cells to use as references for downstream analysis and copy number inference methods.
The output table contains the following columns:

|  |   |
| -------- | -------- |
| `barcodes` | Cell barcode |
| `reference_cell_class` | Indicates if the cell should be uses as a Normal or Tumor cell reference |
| `cellassign_celltype_annotation` | Original annotation as obtained by `CellAssign` in the processed `SingleCellExperiment` object |
| `singler_celltype_annotation` | Original annotation as obtained by `SingleR` in the processed `SingleCellExperiment` object

To create a table of cells use the `--normal_cells` and/or  `--tumor_cells` arguments.
Any cells that have the same annotation as the annotation provided with the `--normal_cells` argument will be labeled as "Normal".
Any cells that have the same annotation as the annotation provided with the `--tumor_cells` argument will be labeled as "Tumor".

In addition to providing a list of cell types, you must also provide an input `SingleCellExperiment` object with columns containing cell type annotations and a `output_filename` with the full path to save the reference file.

The following command specifies using cells annotated as endothelial cells as the normal reference and muscle cells as the tumor reference:

```sh
Rscript 01-select-normal-cells.R \
  --sce_file <path to processed sce file> \
  --normal_cells "endothelial cell" \
  --tumor_cells "smooth-muscle-cell" \
  --output_filename <full path to TSV file to save output table>
```

4. `02-run-cellassign.py`: This script runs [`CellAssign`](https://docs.scvi-tools.org/en/stable/user_guide/models/cellassign.html) on a processed `AnnData` object, requiring an `AnnData` object and a binary reference matrix with cell types as columns and genes as rows as input.
The output will be a TSV file containing a predictions matrix with all cells in the input `AnnData` object as rows and all possible cell types as columns.
The values in the predictions matrix represent the likelihood of each cell being assigned to the respective cell type.

Run this script with the following command:

```sh
python 02-run-cellassign.py \
  --anndata_file <path to anndata file> \
  --output_predictions <path to output tsv file> \
  --reference <path to marker gene reference>
```

5. `03-run-copykat.R`: This script is used to run [`CopyKAT`](https://github.com/navinlabcode/copykat) on a processed `SingleCellExperiment` object.
`CopyKAT` is run on the object using the default parameters specified in `copykat::copykat()` other than the `id.type`, which is set to `E` to account for Ensembl IDs used in the objects available on the Portal.

To run the script with the default parameters, use the following command:

```sh
Rscript run-copykat.R \
  --sce_file <path to processed sce file> \
  --results_dir <name of folder to save results> \
  --threads 4
```

`CopyKAT` can also accept a table output from `01-select-cell-types.R` that contains the `barcodes` and `reference_cell_class` columns using the `--reference_cell_file` argument.
Any cells where `reference_cell_class` is `Normal` will be used as the baseline reference for assigning cells as diploid or aneuploid.

The following command runs `CopyKAT` with a reference table indicating cells to use as the normal cell reference:

```sh
Rscript run-copykat.R \
  --sce_file <path to processed sce file> \
  --reference_cell_file <path to file with table of normal cell barcodes> \
  --results_dir <name of folder to save results> \
  --threads 4
```

6. `04-run-infercnv.R`: This script is used to run [`InferCNV`](https://github.com/broadinstitute/inferCNV/wiki) on a processed `SingleCellExperiment` object.
`InferCNV` is run with a gene cutoff of 0.1 and all other default settings.
The heatmap (saved as `.png`), full object from running `InferCNV` and a table with cell by CNV information are saved to a `output_dir` specified at the command line.

All other intermediate files are saved to a `scratch_dir`, which by default is `scratch/infercnv/{library_id}`.

To run `InferCNV`, an [annotations file](https://github.com/broadinstitute/inferCNV/wiki/File-Definitions#sample-annotation-file) containing all cell barcodes and associated annotations must be created.
This file is created as part of this script and saved to the specified path using `--annotations_file`.

A table containing normal cell barcodes can be provided using the `--reference_cell_file` argument (output from `01-select-normal-cells.R`).
If this is the case, all normal cells will be annotated as "reference" and all other cells will be denoted as "unknown".
If no normal cells are provided, then all cells will be labeled as "unknown" and `InferCNV` will be run with `ref_group_names = NULL`.

This script also requires a gene order file (created by `00-make-gene-order-file.R`).

To run this script use the following command:

```sh
Rscript 04-run-infercnv.Rmd \
  --annotations_file <path to save annotations file> \
  --reference_cell_file <path to file with table of normal cell barcodes> \
  --output_dir <full path to folder to save results> \
  --threads 4
```

## Scripts used to annotate tumor cells with `SingleR`

1. `run-singler.R`: This script is used to run `SingleR` on a processed `SingleCellExperiment` object using a previously annotated `SingleCellExperiment` object as the reference for identifying tumor cells.
This script requires both a reference `SingleCellExperiment` containing a column with previously annotated annotations from `SingleR`, `singler_celltype_annotation`, and a reference TSV file labeling cells as either `Tumor`, `Normal`, or `Ambiguous`.

`SingleR` is run with three different references:

1. The reference TSV file with cells labeled as `Tumor` or `Normal` and the `BlueprintEncodeData` from `celldex`.
Here any `Ambiguous` cells are removed from the reference.
2. The reference `SingleCellExperiment` object using the original annotations found in `singler_celltype_annotation` and replacing any cells labeled as `Tumor` in the reference TSV file with `Tumor`.
3. The same reference as in 2 with the addition of `BlueprintEncodeData` from `celldex`.

The full results from running `SingleR` with all three references will be saved as a single `.rds` file in the `--scratch_dir`.
A TSV file with annotations and cell ontology IDs output from `SingleR` will be saved to the specified file using the `--output_file` argument.

To run this script use the following command:

```sh
Rscript run-singler.R \
  --input_sce_file <path to processed sce file to be annotated> \
  --ref_sce_file <path to previously annotated sce file> \
  --ref_annotations_file <path to TSV file with tumor or normal cell type annotations for reference SCE> \
  --output_file <path to save TSV with SingleR annotations>\
  --threads 4
```

The `singler-workflow` folder contains scripts used in the workflow to run `SingleR` on all samples in SCPCP000015.

1. `00-create-annotation-table.R`: This script is used to create a TSV file with the tumor/normal cell type annotations obtained from three libraries that are used to build a reference for running `SingleR`.
The reference libraries are: `SCPCL000822`, `SCPCL000824`, and `SCPCL001114`.
The annotation table in `results/annotation_tables` is read in for each reference library and the desired column containing the annotations we would like to use are selected.

- `SCPCL000822`: Tumor cells were classified by taking the consensus between `CopyKAT` and `InferCNV`.
See `exploratory_analysis/annotation_notebooks/SCPCS000490/SCPCL000822_tumor-cell-validation.Rmd`.
- `SCPCL000824`: Tumor cells were classified by taking the consensus between two approaches - using a learned marker gene expression cutoff from `SCPCL000822` and `AUCell` results from running `aucell-annotation.sh`.
- `SCPCL001114`: Tumor cells were classified by taking the consensus between two approaches - using a learned marker gene expression cutoff from `SCPCL000822` and `AUCell` results from running `aucell-annotation.sh`.

The annotations across all samples are combined to create a single TSV file, `results/annotation_tables/reference-tumor-cell-annotations.tsv`.
This file contains the following columns:

|  |  |
|--|--|
| `library_id` | ScPCA library id |
| `cell_barcode` | Unique cell barcode |
| `tumor_cell_classification` | One of `Tumor`, `Normal`, or `Ambiguous` |
| `method` | Brief description of the method that was used to call tumor cells |

To run this script use the following command:

```sh
Rscript 00-create-annotation-table.R
```

## Scripts used to annotate tumor cells with `AUCell`

The scripts listed here are used to annotate tumor cells using [`AUCell`](https://www.bioconductor.org/packages/release/bioc/html/AUCell.html) and are implemented in the `auc-annotation.sh` workflow.

1. `01-run-aucell.R`: This script is used run `AUCell` with a list of tumor marker genes on a `SingleCellExperiment` object.
`AUCell` can be used to classify cells using a specified AUC value with `--auc_threshold`.
If no `--auc_threshold` is provided, the `AUCell::AUCell_exploreThresholds()` function will be used to determine the AUC value to use for defining tumor cells.

The output will be a TSV file containing the AUC value determined by `AUCell` and the classification (tumor or normal) for each cell barcode.
Optionally, the `--return_auc` flag can be used to print the AUC value used to classify tumor cells to `stdout`.

By default, tumor marker genes are identifed from the genes listed in `references/tumor-marker-genes.tsv`.
To use a different file, you can use the `--marker_genes_file` option.

To run this script using the AUC value determined by `AUCell` and print the determined AUC to `stdout`, use this command:

```sh
Rscript 01-run-aucell.R \
  --sce_file <path to reference sce file> \
  --output_file <path to TSV file to save results> \
  --return_auc
```

To run this script using a pre-defined AUC value, use this command:

```sh
Rscript 01-run-aucell.R \
  --sce_file <path to reference sce file> \
  --auc_threshold <AUC value required to be classified as a tumor cell> \
  --output_file <path to TSV file to save results>
```

2. `02-calculate-gene-set-scores.R`: This script is used to calculate gene set scores from three different EWS-FLI1 target gene lists for a given `SingleCellExperiment` object.
The sum, z-scaled sum, mean, and z-scaled mean is calculated for each gene set from the `logcounts` assay.
This script returns a TSV file with one row per cell and one column per gene set metric.

The following gene sets are used in this script:
  - [`ZHANG_TARGETS_OF_EWSR1_FLI1_FUSION`](https://www.gsea-msigdb.org/gsea/msigdb/human/geneset/ZHANG_TARGETS_OF_EWSR1_FLI1_FUSION.html)
  - [`RIGGI_EWING_SARCOMA_PROGENITOR_UP`](https://www.gsea-msigdb.org/gsea/msigdb/human/geneset/RIGGI_EWING_SARCOMA_PROGENITOR_UP.html?ex=1)
  - [`SILIGAN_TARGETS_OF_EWS_FLI1_FUSION_DN`](https://www.gsea-msigdb.org/gsea/msigdb/cards/SILIGAN_TARGETS_OF_EWS_FLI1_FUSION_DN)

To run this script use the following command:

```sh
Rscript 02-calculate-gene-set-scores.R \
  --sce_file <path to processed sce file> \
  --output_file <full path to TSV file to save results>
```

## Utils

The `utils` folder contains scripts with any helper functions that are used in `scripts` and notebooks for this module.

1. `jaccard-functions.R`: These functions are used to calculate the Jaccard similarity index between groups of cells.
These functions are taken directly from https://github.com/AlexsLemonade/scpca-nf/blob/main/templates/qc_report/celltypes_supplemental_report.rmd.

2. `tumor-validation-helpers.R`: These functions are used in the template notebooks in `template_notebooks/auc-workflow` and the notebooks found in the `exploratory_analysis/annotation_notebooks` folder.
The functions are helpful for exploring and validate tumor cell annotations in individual samples.
