# Results directory

Results files are available in `s3://researcher-654654257431-us-east-2/doublet-detection/results`

- `benchmark-results` contains TSV files with `scDblFinder` and `scrublet` inferences on ground-truth benchmarking datasets
- `benchmark-results/exploratory-notebooks` contains knitted notebooks exploring the benchmarking results
- `scpca-results` contains TSV files with `scDblFinder` results for ScPCA data
    - This directory is organized by project and sample IDs as:`scpca-results/{project id}/{sample id}/{library id}_processed_scdblfinder.tsv`
    - Currently, the directory contains only results for `SCPCP000001`
