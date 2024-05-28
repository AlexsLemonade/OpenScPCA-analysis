This directory contains template files and directories used by OpenScPCA maintainers.

- `analysis-module`: The template directory for a new analysis module
- `jupyter`: Holds a Jupyter notebook template and an `environment.yml` file to use when creating a new analysis module with the flag `--use_jupyter`
- `python`: Holds a Python script template and an `environment.yml` file to use when creating a new analysis module with the flag `--use_python`
- `rmarkdown`: Holds an R Markdown notebook template and an `environment.yml` file to use when creating a new analysis module with the flag `--use_r` or `--use_renv`
- `workflows`
  - `run_renv-module.yml`: The template GitHub Action workflow for running an R-based analysis module with an `renv` environment
  - `run_conda-module.yml`: The template GitHub Action workflow for running an analysis module with a conda environment
  - `docker_module.yml`: The template GitHub Action workflow for test building a module's Docker image
