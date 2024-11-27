This directory contains exploratory notebooks that are only run on a subset of samples, generated when the module workflow script is run as:

```sh
cd /path/to/module
RUN_EXPLORATORY=1 CNV_RUN_EXPLORATORY=1 ./00_run_workflow.sh
```

- `05_copykat_exploration.Rmd`
  - Explores and compares `CopyKAT` results run on a subset of samples
- `06_infercnv_exploration.Rmd`
  - Explores and compares `inferCNV` results run on a subset of samples


Rendered notebooks are stored in sample-specific subfolders.
Note that the rendered exploratory HTML notebooks may not be fully up-to-date with their R Markdown counterparts.
