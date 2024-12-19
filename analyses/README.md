This directory contains the [OpenScPCA analysis modules](https://openscpca.readthedocs.io/en/latest/contributing-to-analyses/analysis-modules/).

Each directory represents a single analysis module.
Each module contains a `README.md` file with information on its scientific goals and how to run it.

## Example modules

Several modules are provided as examples for how to use certain frameworks in OpenScPCA:

| Module name | Purpose
|-------------|---------
| `hello-R` | Demonstrates an example structure of an R-based module
| `hello-python` | Demonstrates an example structure of a Python-based module
| `hello-cluster` | Provides notebooks demonstrating the use of clustering and cluster evaluations functions in the [`rOpenScPCA` package](https://github.com/AlexsLemonade/rOpenScPCA/)
| `seurat-conversion`| Provides notebooks demonstrating the use of functions related to converting SCE to Seurat objects in the [`rOpenScPCA` package](https://github.com/AlexsLemonade/rOpenScPCA/), which also includes functions to convert Ensembl ids to gene symbols