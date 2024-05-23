# Scripts

This directory contains all scripts used in the module.

1. `00_format-benchmark-data.R` formats ground-truth datasets used to explore doublet detection methods into SCE and AnnData versions.

Use this command to run the bash script:

```sh
00_format-benchmark-data.R \
    --dataset_name <name of dataset to format> \
    --input_dir <directory with raw data file from Zenodo repository> \
    --output_dir <directory to export formatted files>
```

2. `01a_detect-doublets.R` and `01b_detect-doublets.py` detect doublets using R-based (`scDblFinder`) and Python-based (`scrublet`) methods, respectively.
- `01a_detect-doublets.R` takes an SCE file and exports a TSV file with `scDblFinder` results
- `01b_detect-doublets.py` takes an AnnDAta file and exports a TSV fileswith `scrublet` results

Use these commands to run each script:

```sh
# Run 01a_detect-doublets.R
Rscript 01a_detect-doublets.R \
    --dataset_name <name of dataset to format> \
    --data_dir <directory to SCE file> \
    --results_dir <directory to export TSV result file> \
    --cores <optionally, the number of cores to use. Default is 4> \
    --random_seed <optionally, a random seed. Default is 2024>

# Run 01b_detect-doublets.py
python3 01b_detect-doublets.py \
    --dataset_name <name of dataset to format> \
    --data_dir <directory to AnnData file> \
    --results_dir <directory to export TSV result file> \
    --random_seed <optionally, a random seed. Default is 2024>
```
