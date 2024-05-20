# Scripts

This directory contains all scripts used for cell typing Ewing sarcoma samples from SCPCP000015.

1. `run-cellassign.py`: This script runs [`CellAssign`](https://docs.scvi-tools.org/en/stable/user_guide/models/cellassign.html) on a processed `AnnData` object, requiring an `AnnData` object and a binary reference matrix with cell types as columns and genes as rows as input.
The output will be a TSV file containing a predictions matrix with all cells in the input `AnnData` object as rows and all possible cell types as columns.
The values in the predictions matrix represent the likelihood of each cell being assigned to the respective cell type.

Run this script with the following command:

```sh
python run-cellassign.py \
  --anndata_file <path to anndata file> \
  --output_predictions <path to output tsv file> \
  --reference <path to marker gene reference>
```

2. `select-cell-types.R`: This script is used to generate a table of normal and tumor cells to use as references for downstream analysis and copy number inference methods.
The files produced from this script are saved in `references/cell_lists/<sample_id>/<library_id>_reference-cells.tsv`.
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
Rscript select-normal-cells.R \
  --sce_file <path to processed sce file> \
  --normal_cells "endothelial cell" \
  --tumor_cells "smooth-muscle-cell" \
  --output_filename "references/cell_lists/<sample_id>/<library_id>_reference-cells.tsv"
```

3. `run-copykat.R`: This script is used to run [`CopyKAT`](https://github.com/navinlabcode/copykat) on a processed `SingleCellExperiment` object.
`CopyKAT` is run on the object using the default parameters specified in `copykat::copykat` other than the `id.type`, which is set to `E` to account for Ensembl IDs used in the objects available on the Portal.
The files output from running `CopyKAT` are saved in `results/copykat/{library_id}/{copykat_results_prefix}`, where `copykat_results_prefix` is specified at the command line and `library_id` is taken from the `library_id` stored in the processed `SingleCellExperiment` object.

To run the script with the default parameters, use the following command:

```sh
Rscript run-copykat.R \
  --sce_file <path to processed sce file> \
  --copykat_results_prefix <name of folder to save results> \
  --threads 4
```

`CopyKAT` can also accept a vector of cells to use as the baseline reference for assigning cells as diploid or aneuploid.
To provide a list of cell types to use as a baseline, provide a `.txt` file with the list of barcodes corresponding to normal cells using the `--normal_cells` argument.
This file can be created by running `select-normal-cells.R`.

The following command runs `CopyKAT` with an input list of cells to use as the normal cell reference:

```sh
Rscript run-copykat.R \
  --sce_file <path to processed sce file> \
  --normal_cells <path to file with list of normal cell barcodes> \
  --copykat_results_prefix <name of folder to save results>
  --threads 4
```

4. `make-gene-order-file.R`: This script is used to generate the [gene order file](https://github.com/broadinstitute/inferCNV/wiki/File-Definitions#gene-ordering-file) needed for running `InferCNV` with `run-infercnv.R`.
This script downloads the GTF file used to create the original index used in `scpca-nf` (Ensembl v104) from the public bucket `s3://scpca-references`.
This file is then converted to the required gene order formatted file for `InferCNV` and saved to `references/infercnv_refs`.

To run this script use the following command:

```sh
Rscript make-gene-order-file.R
```

5. `run-infercnv.R`: This script is used to run [`InferCNV`](https://github.com/broadinstitute/inferCNV/wiki) on a processed `SingleCellExperiment` object.
`InferCNV` is run with a gene cutoff of 0.1 and all other default settings.
The files output from running `InferCNV` are saved in `results/infercnv/{library_id}/{infercnv_results_prefix}`, where `infercnv_results_prefix` is specified at the command line and `library_id` is taken from the `library_id` stored in the processed `SingleCellExperiment` object.

To run `InferCNV`, an [annotations file](https://github.com/broadinstitute/inferCNV/wiki/File-Definitions#sample-annotation-file) containing all cell barcodes and associated annotations must be created.
This file is created as part of this script and saved to the specified path using `--annotations_file`.
If a list of normal cell barcodes are provided (output from `select-normal-cells.R`), then all normal cells will be annotated as "reference" and all other cells will be denoted as "unknown".
If no normal cells are provided, then all cells will be labeled as "unknown" and `InferCNV` will be run with `ref_group_names = NULL`.

This script also requires a gene order file (created by `make-gene-order-file.R`).

To run this script use the following command:

```sh
Rscript run-infercnv.Rmd \
  --annotations_file <path to save annotations file> \
  --normal_cells <path to file with list of normal cell barcodes> \
  --results_dir <name of folder to save results> \
  --threads 4
```

## Utils

The `utils` folder contains Rscripts with any helper functions that are used in `scripts` and notebooks for this module.

1. `jaccard-functions.R`: These functions are used to calculate the Jaccard similarity index between groups of cells.
These functions are taken directly from https://github.com/AlexsLemonade/scpca-nf/blob/main/templates/qc_report/celltypes_supplemental_report.rmd.
