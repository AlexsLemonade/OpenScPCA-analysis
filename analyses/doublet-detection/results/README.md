# Results directory

- Files in `benchmark-results` are created by `scripts/00_benchmark-methods.R`, using 4 cores for `scDblFinder`
  - `benchmark_runtimes.tsv` contains time (in seconds) benchmarks for each of four doublet detect across libraries, as well as some library metadata
  - `benchmark_results.tsv` contains output from those four doublet detection methods where each line is a cell from a given library
