# Results directory

- `benchmark-times_{x}-cores.tsv` is created by `scripts/00_benchmark-methods.R`
  - This file contains time (in seconds) benchmarks for each of four doublet detect methods as run on `{x}` cores across 4 processed samples per ScPCA project of varying library sizes
  - Sample diagnoses and processed cell counts are also included
    - Note that only `scDblFinder` uses parallel processing

