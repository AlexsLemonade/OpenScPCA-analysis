# Managing software environments

A key way the OpenScPCA project ensures reproducibility is through the use of module-specific software environments.
Using software environments has several benefits:

- Software dependencies, including both packages and stand-alone command line tools, can be pinned to specific versions.
As different versions of the same software often have different behaviors/usage and can produce different results, pinning versions ensures that results are much more reproducible.
- Environments make it easier for other researchers to run a module.
Analysis modules will have files specifying the specific software dependencies and versions which researchers can use to re-create the environment to reproducibility run the module's code.

We specifically recommend these two package managers to specify module-specific software environments:

- The R package manager [`renv`](https://rstudio.github.io/renv/) for modules containing R code.
  - Before working with `renv`, we encourage you to read [the excellent introduction](https://rstudio.github.io/renv/articles/renv.html) for more information.
- The package manager [conda](https://docs.conda.io/en/latest/), along with [`conda-lock`](https://conda.github.io/conda-lock/), for modules containing Python code.
Conda can also be used to manage standalone software packages that do not depend on a specific language.

When you [create a module with `create-analysis-module.py`](../../contributing-to-analyses/analysis-modules/creating-a-module.md), you can use one (or more!) of several [flags](../../contributing-to-analyses/analysis-modules/creating-a-module.md#module-creation-script-flags) to establish your module with an initialized software environment.

## Updating software environments

As analysis modules mature, the Data Lab will activate several [GitHub Action workflows](../workflows/index.md) that ensure module reproducibility.
Once these workflows are activated, you will need to file _two PRs_ every time you wish to update your software environment, including adding or updating a package used in R or Python: 

1. First, file a PR that updates your software environment, e.g. `renv.lock`, `conda.lock`, or `Dockerfile` files
1. After that first PR is been merged, you can then file a second PR that contains your code changes that use the updated environment

For more details on this process, please refer to our [documentation on updating Docker images](../docker/updating-images.md).