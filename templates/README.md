This directory contains template files and directories used by OpenScPCA maintainers.

- `analysis-module`: The template directory for a new analysis module
- `jupyter`: Holds a template Jupyter notebook template and an `environment.yml` file to use when creating a new analysis module with the flag `--use_jupyter`
- `python`: Holds a python script template and an `environment.yml` file to use when creating a new analysis module with the flag `--use_python`
- `rmarkdown`: Holds an R Markdown notebooks template and an `environment.yml` file to use when creating a new analysis module with the flag `--use_r` or `--use_renv`
- `run-analysis-module-R.yml`: The template GitHub Action workflow for running an R-based analysis module with an `renv` environment
- `run-analysis-module-python.yml`: The template GitHub Action workflow for running a python-based analysis module with a conda environment
- [Forthcoming] `run-analysis-module-docker.yml`: The template GitHub Action workflow for running an analysis module using the module's Dockerfile
