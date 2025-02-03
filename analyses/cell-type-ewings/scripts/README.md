# Scripts

This directory contains all scripts used for cell typing Ewing sarcoma samples from SCPCP000015.

<!-- START doctoc generated TOC please keep comment here to allow auto update -->
<!-- DON'T EDIT THIS SECTION, INSTEAD RE-RUN doctoc TO UPDATE -->
**Table of Contents** 

- [Scripts used to annotate tumor cells with `AUCell` and `SingleR`](#scripts-used-to-annotate-tumor-cells-with-aucell-and-singler)
- [Scripts used to annotate tumor cells with `SingleR`](#scripts-used-to-annotate-tumor-cells-with-singler)
- [Scripts used in the clustering workflow](#scripts-used-in-the-clustering-workflow)
- [Scripts used for running `AUCell` with `EWS-FLI` gene signatures](#scripts-used-for-running-aucell-with-ews-fli-gene-signatures)
- [Scripts used in the CNV annotation workflow](#scripts-used-in-the-cnv-annotation-workflow)
- [Utils](#utils)

<!-- END doctoc generated TOC please keep comment here to allow auto update -->


## Scripts used to annotate tumor cells with `AUCell` and `SingleR`

The scripts listed here are used to annotate tumor cells using [`AUCell`](https://www.bioconductor.org/packages/release/bioc/html/AUCell.html) and `SingleR` and are implemented in the `auc-singler-annotation.sh` workflow.

1. `01-run-aucell.R`: This script is used run `AUCell` with a list of tumor marker genes on a `SingleCellExperiment` object.
   `AUCell` can be used to classify cells using a specified AUC value with `--auc_threshold`.
   If no `--auc_threshold` is provided, the `AUCell::AUCell_exploreThresholds()` function will be used to determine the AUC value to use for defining tumor cells.

The output will be a TSV file containing the AUC value determined by `AUCell` and the classification (tumor or normal) for each cell barcode.
Optionally, the `--return_auc` flag can be used to print the AUC value used to classify tumor cells to `stdout`.

By default, tumor marker genes are identified from the genes listed in `references/tumor-marker-genes.tsv`.
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

2. `02-calculate-gene-set-scores.R`: This script is used to calculate gene set scores from three different `EWS-FLI1` target gene lists for a given `SingleCellExperiment` object.
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

3. `03-generate-tumor-ref.R`: This script is used to create a merged `SingleCellExperiment` object containing all high-confidence tumor cells from all samples in `SCPCP000015`.
   This script reads in all processed `SingleCellExperiment` objects and all results from running `01-run-aucell.R` for `SCPCP000015`.
   Any cells that are labeled as "Tumor" in the `auc_classification` column of the results are kept and all other cells are discarded.
   The tumor cells from all objects are then merged into a single `SingleCellExperiment` object that will be used as a reference for running `SingleR`.

Running the script using the default options will save the merged reference object to `scratch/tumor-ref-singler.rds`.
A different path for the output file can be specified using the `--output_reference_file` argument.

To run this script use the following command:

```sh
Rscript 02-generate-tumor-ref.R
```

4. `04-run-singler.R`: This script is used to run `SingleR` on any non-confident tumor cells from a processed `SingleCellExperiment` object.
   The input `SingleCellExperiment` object is run through `SingleR` using the reference from `03-generate-tumor-ref.R`, `BlueprintEncodeData` from `celldex`, and `HumanPrimaryCellAtlasData` from `celldex` as the references.
   Note that any tumor cells from the same participant as the library being annotated are removed from the tumor reference object prior to running `SingleR`.

The output is a TSV file with one row for each of the cells in the original `SingleCellExperiment` object being annotated and the following columns:

|                      |                                                                                                                           |
| -------------------- | ------------------------------------------------------------------------------------------------------------------------- |
| `barcodes`           | Unique cell barcode                                                                                                       |
| `singler_annotation` | The human readable term associated with the ontology term identifier for the associated cell type or `tumor-<library_id>` |
| `singler_ontology`   | The cell ontology identifier associated with the cell type or `tumor-<library_id>`                                        |
| `aucell_annotation`  | Either `tumor` or `normal` as determined by running `AUCell` in `aucell-annotation.sh`                                    |

If the cell is a tumor cell it will be labeled with `tumor-library_id` in both the ontology and annotation columns, where `library_id` represents the library that the tumor cell resembled the most.

Running the script using the default options will use the merged object output from `03-generate-tumor-ref.R` saved in `scratch/tumor-ref-singler.rds`.
A different path to this file can be specified using the `--tumor_reference_file` argument.

To run this script use the following command:

```sh
Rscript 04-run-singler.R \
  --sce_file <path to processed sce file to be annotated> \
  --output_annotations_file <path to TSV file to save annotations> \
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

## Scripts used in the clustering workflow 

1. `01-clustering.R`: This script is used to calculate clusters for a given library across a set of clustering parameters. 
The script supports both Louvain and Leiden clustering (either CPM or modularity objective function). 
The nearest neighbors and resolution parameters can be tested over a range of values, specified at the command line.

When running the script, you must specify which clustering algorithm(s) to consider by using the `--louvain`, `--leiden_mod`, and/or `--leiden_cpm` flags at the command line. 

By default the script is run with the following default parameters: 
- Nearest neighbors: 5, 10, 15, 20, 25, 30, 35, 40
- Resolution (Louvain): 0.5, 1, 1.5
- Resolution (Leiden modularity): 0.5, 1, 1.5
- Resolution (Leiden CPM): 0.001, 0.005, 0.01

To run this script using the default parameter ranges use the following command: 

```sh
Rscript 01-clustering.R \
  --sce_file <path to processed SCE file> \
  --output_file <path to TSV file to save clustering results> \
  --louvain --leiden_mod --leiden_cpm \ #flags to indicate using all three algorithms
  --threads 4 \
  --seed 2024
```

To modify the range of values for any parameters use the following command (note that resolution ranges are specified separately for louvain, leiden modularity, and leiden CPM): 

```sh
Rscript 01-clustering.R \
  --sce_file <path to processed SCE file> \
  --output_file <path to TSV file to save clustering results> \
  --louvain --leiden_mod --leiden_cpm \ #flags to indicate using all three algorithms
  --nn_range "5,10,15" \
  --louvain_res_range "1,1.5" \ # resolution to test for louvain
  --mod_res_range "1,1.5" \ # resolution to test for leiden modularity
  --cpm_res_range ".01,.001" \ # resolution to test for leiden cpm 
  --threads 4 \
  --seed 2024
```

## Scripts used for running `AUCell` with `EWS-FLI` gene signatures 

1. `01-aucell.R`: This script is used to run `AUCell` with a set of custom gene signatures on a single processed object. 
By default, all gene signatures in [references/gene_signatures](../references/gene_signatures/) are used alongside a set of gene signatures from `MSigDB` associated with high and low EWS-FLI1 expression. 
The full list of gene signatures used can be found in [the references `README.md`](../references/README.md#gene-signatures). 

`AUCell` is run for each gene signature and AUC values along with the AUC threshold reported by `AUCell` are saved to a TSV file. 
`AUCell` is run with an `aucMaxRank` value equal to 1% of the detected genes in the processed object. 
This can be changed using the `--max_rank_threshold` parameter. 

By default, this script uses 4 CPUs. 

To run this script on a single library use the default parameters use the following command: 

```sh
Rscript 01-aucell.R \
  --sce_file <path to processed SCE file> \
  --output_file <path to TSV file to save AUC results>
```

To run this script with a merged object use the following command: 

```sh
Rscript 01-aucell.R \
  --sce_file <path to processed SCE file> \
  --output_file <path to TSV file to save AUC results> \
  --is_merged
```

## Scripts used in the CNV annotation workflow

**NOTE:** This workflow is no longer used in the cell type annotation analysis, has been removed from CI, and is no longer maintained! 

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

|                                  |                                                                                                |
| -------------------------------- | ---------------------------------------------------------------------------------------------- |
| `barcodes`                       | Cell barcode                                                                                   |
| `reference_cell_class`           | Indicates if the cell should be uses as a Normal or Tumor cell reference                       |
| `cellassign_celltype_annotation` | Original annotation as obtained by `CellAssign` in the processed `SingleCellExperiment` object |
| `singler_celltype_annotation`    | Original annotation as obtained by `SingleR` in the processed `SingleCellExperiment` object    |

To create a table of cells use the `--normal_cells` and/or `--tumor_cells` arguments.
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
Rscript 04-run-infercnv.R \
  --annotations_file <path to save annotations file> \
  --reference_cell_file <path to file with table of normal cell barcodes> \
  --output_dir <full path to folder to save results> \
  --threads 4
```


## Utils

The `utils` folder contains scripts with any helper functions that are used in `scripts` and notebooks for this module.

1. `jaccard-functions.R`: These functions are used to calculate the Jaccard similarity index between groups of cells.
   These functions are taken directly from https://github.com/AlexsLemonade/scpca-nf/blob/main/templates/qc_report/celltypes_supplemental_report.rmd.

2. `tumor-validation-helpers.R`: These functions are used in the template notebooks in `template_notebooks/auc-singler-workflow` and the notebooks found in the `exploratory_analysis/annotation_notebooks` folder.
   The functions are helpful for exploring and validate tumor cell annotations in individual samples.
   
3. `clustering-functions.R`: These functions are used in the template notebook for exploring clustering, `template_notebooks/clustering-workflow/03-clustering.Rmd`. 
The functions included are used to calculate and plot clusters and cluster statistics. 

