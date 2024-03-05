# test-module

This analysis module is intended for testing during `OpenScPCA-analysis` setup.
It contains two notebooks:

- `scripts/01_plot-iris.Rmd`
  - Exports a density plot of `iris` sepal lengths to `plots`
- `scripts/02_save-coffee.ipynb`
  - Reads from URL and exports a TSV to `results` (not under VC)

The `renv.lock` file has been updated with dependencies used by `scripts/01_plot-iris.Rmd`.

The `environment.yaml` file has currently _not_ been updated with dependencies, but note that these dependencies were installed:
```
conda install nb_conda # provides notebooks with ability to choose conda environment
conda install pandas
```
