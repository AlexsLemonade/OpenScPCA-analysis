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

When you [create a module with `create-analysis-module.py](../../contributing-to-analyses/analysis-modules/creating-a-module.md), you can use one (or more!) of several [flags](../../contributing-to-analyses/analysis-modules/creating-a-module.md#module-creation-script-flags) to establish your module with an initialized software environment.
