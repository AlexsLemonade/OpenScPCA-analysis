# Scripts

This directory contains all scripts used in the module.

1. `00a_download-benchmark-data.sh` and `00b_format-benchmark-data.R` are used to download and format, respectively, ground-truth datasets used to explore doublet detection methods.
Both the original datasets and SCE and AnnData versions are saved in `../scratch/benchmark_datasets`.
Note that `00a_download-benchmark-data.sh` calls `00b_format-benchmark-data.R`; `00b_format-benchmark-data.R` is not expected to be invoked on its own.

Use this command to run the bash script:
```sh
bash 00a_download-benchmark-data.sh
```
