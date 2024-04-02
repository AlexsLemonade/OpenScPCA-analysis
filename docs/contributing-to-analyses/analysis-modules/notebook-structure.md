# Structuring your notebooks

Notebook formats such as [R Markdown](https://rmarkdown.rstudio.com/index.html) or [Jupyter](https://jupyter-notebook.readthedocs.io/en/latest/) are a great way to document an analysis and share it with others.
We have found that following some common patterns can enhance sharing and reproducibility, so we suggest the following structure when you create a notebook.


## R Markdown notebooks

We will illustrate the structure of an R Markdown notebook using the example notebook found in the [`hello-R`](https://github.com/AlexsLemonade/OpenScPCA-analysis/tree/main/analyses/hello-R) analysis module: [`hello.Rmd`](https://github.com/AlexsLemonade/OpenScPCA-analysis/tree/main/analyses/hello-R/hello.Rmd).


### Headers

R Markdown notebooks start with a YAML header that should contain a `title`, `author`, and `date`, as shown in this example:

```yaml
---
title: Hello OpenScPCA
author: Childhood Cancer Data Lab, ALSF
date: "`r Sys.Date()`"
output: html_notebook
---
```

The date is here set dynamically in the output to the date that the analysis was run, using the `Sys.Date()` function, but you can set it manually if preferred.

The `output` field specifies the output format of the notebook, which will usually be `html_notebook` for the analyses in the OpenScPCA project.


### Introduction

The remainder of the notebook will be a mix of Markdown and R code chunks.
We suggest starting with an **Introduction** section that briefly describes the analysis and its purpose, to help orient readers.

You should include a bit of information about the input data, the analysis steps, and the expected output.


### Setup

Following the Introduction, a **Setup** section should handle loading R packages and defining paths for input and output files.


#### Loading packages

The notebook should not download any R packages.
All packages should already be installed on the system running the notebook, and can be separately tracked using [`renv`](../determining-requirements/determining-software-requirements.md#using-renv).

We generally recommend keeping the number of packages loaded with `library()` to a minimum and using the `package::function()` syntax to call functions, to make it clear which package each function comes from.
For some common packages, however, this can become burdensome.

It is often also convenient to suppress the sometimes verbose messages that packages print on loading using `suppressPackageStartupMessages()`, as shown below:

```r
# load required libraries
suppressPackageStartupMessages({
  library("SingleCellExperiment")
  library("ggplot2")
})
```


#### Setting paths

Defining paths to all input and output files at the start of the notebook makes it much easier to users to understand the structure of the analysis and to modify the paths if needed.

R notebooks by default use paths relative to the location of the notebook file, but this can lead to errors if the notebook file itself is moved, so we recommend defining paths relative to either the OpenScPCA project root or the analysis module root.

In the `hello.Rmd` notebook, the paths are defined using the [`rprojroot` package](https://rprojroot.r-lib.org), which has a number of convenient functions for finding the root of a project based on its contents.
In the example below, we find the OpenScPCA project root by looking for the `.git` directory, and the analysis module root by looking for the `renv` files that were set up in that project.

```r
# Find the repository and module root directories
repo_root <- rprojroot::find_root(rprojroot::is_git_root)
module_root <- rprojroot::find_root(rprojroot::is_renv_project)
```

Alternatively, the module root could be found relative to the repository root, as shown below:

```r
repo_root <- rprojroot::find_root(rprojroot::is_git_root)
module_root <- file.path(repo_root, "analyses", "hello-R")
```

Note that we use `file.path()` to construct paths, which is a platform-independent way to concatenate directories and filenames for safe navigation.

All file paths should be defined as relative paths, with the exception of any input data that might be downloaded as part of the analysis, which may be referenced by a URL.

!!! note
    When working on an analysis, it is quite common to find that you need a file that you had not anticipated at the start.
    No problem!
    Just add that path to the Setup section with the other paths, rather than defining it later in the notebook.


### Defining custom functions

Followith the Setup section, it is often useful to have a **Functions** section where all custom functions that will be used later in the notebook are defined.
Keeping the functions in a central place in the notebook makes it easier to find and modify them later, and also makes it easier to reuse them in other analyses.
Each function should be defined in a separate code chunk, and should be documented with comments to explain what the function does does and how to use it, including inputs and outputs.

In `hello.Rmd`, we define a simple function count cells in a `SingleCellExperiment` object:

```r
count_sce <- function(sce_file) {
  # Count cells in an sce file
  # Args:
  # - sce_file: path to the SingleCellExperiment file
  # Returns a data frame with the following columns:
  # - project_id
  # - sample_id
  # - library_id
  # - n_cells: number of cells in the library

  ...
}
```

We generally do _not_ recommend using `source()` to load functions from external files, as this can make it harder to keep track of where functions are defined and can lead to errors if the file is moved or renamed.
It also means that the notebook is not self-contained, and functions will not be present in the output html files, which can make it harder to share and reproduce the analysis.


### Analysis steps

Once all of the inctions are defined, the remainder of the notebook should be a series of code chunks that perform the analysis steps, with Markdown text to explain what each step is doing and why.
Use headings and subheadings as appropriate to break up the analysis into logical sections, and include plots in the notebook to help illustrate the results.

Code chunks should still contain comments as needed to explain logic and implementation, with the Markdown text providing a higher-level overview of the analysis steps, including any interpretation of results.


### Session info

The final section of every notebook should be a **Session info** section that uses the [`sessionInfo()` function](../determining-requirements/determining-software-requirements.md#using-sessioninfo) to print out the versions of all packages used in the analysis.

!!! tip
    In `hello.Rmd`, we use a package with a slightly nicer output format: `sessioninfo::session_info()`; either is fine!


## Jupyter notebooks
