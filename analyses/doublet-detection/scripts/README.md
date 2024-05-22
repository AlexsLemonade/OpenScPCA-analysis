# Scripts

This directory contains all scripts used in the module.

1. `00a_download-benchmark-data.sh` and `00b_format-benchmark-data.R` are used to download and format, respectively, ground-truth datasets used to explore doublet detection methods.
Both the original datasets and SCE and AnnData versions are saved in `../scratch/benchmark_datasets`.
Note that `00a_download-benchmark-data.sh` calls `00b_format-benchmark-data.R`; `00b_format-benchmark-data.R` is not expected to be invoked on its own.

Use this command to run the bash script:
```sh
bash 00a_download-benchmark-data.sh
```

2. `01a_detect-doublets.R` and `01b_detect-doublets.py` detect doublets using R-based (`scDblFinder`) and Python-based (`scrublet`) methods, respectively.
Each script exports results to `results/benchmark_results`:
- `01a_detect-doublets.R` exports SCE files with `scDblFinder` results
- `01b_detect-doublets.py` exports TSV files with `scrublet` results

Use these command to run each script:
```sh
Rscript 01a_detect-doublets.R
python3 01a_detect-doublets.py
```