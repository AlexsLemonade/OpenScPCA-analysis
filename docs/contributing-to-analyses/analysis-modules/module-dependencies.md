# Managing packages and dependencies

When developing an analysis module, you will generally need to install packages and/or software tools.
These are known as your module-specific dependencies, i.e., the collection of software languages, packages, or libraries required for your analysis to run.

The OpenScPCA project uses package management tools to track and manage module-specific dependencies, as well as pin those dependency versions to ensure long-term reproducibility:

- Use the R package manager [**`renv`**](https://rstudio.github.io/renv/) to manage and add new R package dependencies.
- Use [**conda**](https://docs.conda.io/en/latest/) to manage and add new Python package dependencies, as well as standalone software tools.

<div class="grid" markdown>

[Read more about using `renv`](../../ensuring-repro/managing-software/using-renv.md){ .md-button .md-button-- }

[Read more about using conda](../../ensuring-repro/managing-software/using-conda.md){ .md-button }

</div>

When you [create a module with `create-analysis-module.py`](../../contributing-to-analyses/analysis-modules/creating-a-module.md), you can use one (or more!) of several [flags](../../contributing-to-analyses/analysis-modules/creating-a-module.md#module-creation-script-flags) to establish your module with an initialized `renv` and/or conda software environment, depending on your needs.

As you develop your module, you'll also need to [add documentation](./documenting-analysis.md) about how to activate your module's software environment.

