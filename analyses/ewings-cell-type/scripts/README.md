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

2. `run-copykat.R`: This script is used to run [`CopyKAT`](https://github.com/navinlabcode/copykat) on a processed `SingleCellExperiment` object.
`CopyKAT` is run on the object using the default parameters specified in `copykat::copykat` other than the `id.type`, which is set to `E` to account for Ensembl IDs used in the objects available on the Portal.
The files output from running `CopyKAT` are saved in `output_dir/{library_id}`, where `output_dir` is specified at the command line and `library_id` is taken from the `library_id` stored in the processed `SingleCellExperiment` object.
The default `output_dir` will save results in `results/copykat`.

To run the script with the default parameters, use the following command:

```sh
Rscript run-copykat.R \
  --sce_file <path to processed sce file> \
  --threads 4
```

`CopyKAT` can also accept a vector of cells to use as the baseline reference for assigning cells as diploid or aneuploid.
To provide a list of cell types to use as a baseline, use the `--singler_normal_cells` and/or  `--cellassign_normal_cells` arguments.
If either of these arguments are provided, then `CopyKAT` will be run twice - once without any normal cells and once with normal cells as the reference.
If using both `--singler_normal_cells` and `--cellassign_normal_cells`, the cells present in both groups will be used as the normal reference.

The following command specifies using cells annotated as endothelial cells in both `SingleR` and `CellAssign` as the normal reference:

```sh
Rscript run-copykat.R \
  --sce_file <path to processed sce file> \
  --singler_normal_cells "endothelial cell" \
  --cellassign_normal_cells "Endothelial cells" \
  --threads 4
```

## Utils

The `utils` folder contains Rscripts with any helper functions that are used in `scripts` and notebooks for this module.

1. `jaccard-functions.R`: These functions are used to calculate the Jaccard similarity index between groups of cells.
These functions are taken directly from https://github.com/AlexsLemonade/scpca-nf/blob/main/templates/qc_report/celltypes_supplemental_report.rmd.
