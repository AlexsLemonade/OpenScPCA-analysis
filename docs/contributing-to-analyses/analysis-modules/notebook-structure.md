# Structuring analysis notebooks

Notebook formats such as [R Markdown](https://rmarkdown.rstudio.com/index.html) or [Jupyter](https://jupyter-notebook.readthedocs.io/en/latest/) are a great way to document an analysis and share it with others.
We have found that following some common patterns can enhance sharing and reproducibility, so we suggest the following overall structure when you create a notebook:

- **Headers**: Start with a header that includes a title, author, and date.
- **Introduction**: Briefly describe the analysis and its purpose.
- **Setup**: Load required packages and define paths for input and output files.
- **Functions**: Define any custom functions used in the analysis.
- **Analysis**: This is where the bulk of the analysis should appear.
    - Intersperse code, to do the work, and text, to explain the steps and interpret results.
    - Use headings and subheadings to break up the analysis into logical sections.
- **Session info**: Print out the versions of all packages used in the analysis.

Below we provide more detail about each of these sections for [R Markdown](#r-markdown-notebooks) and [Jupyter](#jupyter-notebooks) notebooks.

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

The `date` field is set here to show the date the notebook was run using the `Sys.Date()` function, but you can set the value manually if preferred.
If you do set the date manually, be sure to keep it current as the notebook is modified.

The `output` field specifies the output format of the notebook, which we recommend setting to `html_notebook` for the analyses in the OpenScPCA project.


### Introduction

The remainder of the notebook should be a mix of Markdown and R code chunks.
We suggest starting with an **Introduction** section that briefly describes the analysis and its purpose, to help orient readers.

You should include a bit of information about the input data, the analysis steps, and the expected output.


### Setup

Following the Introduction, a **Setup** section should load R packages and define paths for input and output files.


#### Loading packages

The notebook should not download any R packages; in other words, you should never include the code `install.packages()` in your notebook.
All packages should already be installed on the system running the notebook, and can be separately tracked using [`renv`](../determining-requirements/determining-software-requirements.md#using-renv).

We generally recommend keeping the number of packages loaded with `library()` to a minimum and using the `package::function()` syntax to call functions, to make it clear which package each function comes from.
However, for some commonly invoked packages such as `SingleCellExperiment` and `ggplot2` this can become burdensome, making loading the package with `library()` more convenient.

It is often also convenient to suppress the sometimes verbose messages that packages print on loading using `suppressPackageStartupMessages()`, as shown below:

```r
# load required libraries
suppressPackageStartupMessages({
  library("SingleCellExperiment")
  library("ggplot2")
})
```


#### Setting paths

Defining paths to all input and output files at the start of the notebook makes it much easier for users to understand the analysis structure and to modify the paths if needed.

Do not use absolute file paths, as this can make it harder to share and reproduce the analysis.
Instead, define paths relative to the root of the OpenScPCA project or the root of the analysis module.
While R notebooks allow you to set paths relative to location of the notebook file, this can be brittle if the notebook is ever moved to a new location.


In the `hello.Rmd` notebook, we define root paths using the [`rprojroot` package](https://rprojroot.r-lib.org), which has a number of convenient functions for finding the root of a project based on its contents.
In the example below, we find the OpenScPCA project root by looking for the `.git` directory, and the analysis module root by looking for the `renv` files that were set up in that module.

```r
# Find the repository and module root directories
repo_root <- rprojroot::find_root(rprojroot::is_git_root)
module_root <- rprojroot::find_root(rprojroot::is_renv_project)
```

Alternatively, if your module does not use `renv`, the module root could be found relative to the OpenScPCA repository root, as shown below:

```r
repo_root <- rprojroot::find_root(rprojroot::is_git_root)
module_root <- file.path(repo_root, "analyses", "hello-R")
```

Note that we use `file.path()` to construct paths, which is a platform-independent way to define directories and filenames for safe navigation.


!!! note
    When working on an analysis, it is quite common to find that you need a file that you had not anticipated when you started.
    No problem!
    Just add that path to the Setup section with the other paths, rather than defining it later in the notebook.


### Defining custom functions

Following the Setup section, it is often useful to have a **Functions** section where you define any custom functions you write for use later in the notebook.
Keeping the functions in a central place in the notebook makes it easier to find and modify them later.
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

!!! note
    We generally do _not_ recommend using `source()` to load functions from external files, as this can make it harder to keep track of where functions are defined and can lead to errors if the file is moved or renamed.

    It also means that the notebook is not self-contained, and functions will not be present in the output html files, which can make it harder to share and reproduce the analysis.


### Analysis steps

Once all of the functions are defined, the remainder of the notebook should be a series of code chunks that perform the analysis steps, with [Markdown text](../../software-platforms/general-tools/writing-in-markdown.md) to explain what each step is doing and why.
Use headings and subheadings as appropriate to break up the analysis into logical sections, and include plots in the notebook to help illustrate the results.

Code chunks should still contain comments to explain logic and implementation.
The Markdown text should focus on providing a higher-level overview of the analysis steps, including any interpretation of results.


### Session info

The final section of every R Markdown notebook should be a **Session info** section that uses the [`sessionInfo()` function](../determining-requirements/determining-software-requirements.md#using-sessioninfo) to print out the versions of all packages used in the analysis.

!!! tip
    In `hello.Rmd`, we use [the `sessioninfo` package](https://sessioninfo.r-lib.org) to generate a slightly nicer output format using the `sessioninfo::session_info()` function; either approach is fine!


## Jupyter notebooks

We will illustrate the structure of a Jupyter notebook using the example notebook found in the [`hello-python`](https://github.com/AlexsLemonade/OpenScPCA-analysis/tree/main/analyses/hello-python) analysis module: [`hello.ipynb`](https://github.com/AlexsLemonade/OpenScPCA-analysis/tree/main/analyses/hello-python/hello.ipynb).


### Headers

Jupyter notebooks do not have a YAML header like R Markdown notebooks, but you should still begin with a Markdown cell at the top of the notebook that includes a title, author, and date, as shown in this example:

```markdown
# Hello OpenScPCA

Childhood Cancer Data Lab, ALSF

2024-03-19
```

The date should be kept current as the notebook is modified.


### Introduction

Following the header cell, an **Introduction** section should briefly describe the analysis and its purpose, to help orient readers.

You should include a bit of information about the input data, the analysis steps, and the expected output.


### Setup

Following the Introduction, a **Setup** section should include cells to load Python packages and define paths for input and output files.


#### Loading Python packages

Load all Python packages in a cell or cells at the start of the notebook using `import` statements.
For example, you will always want to import the [`session_info` package](#session-info-1) to document the versions of all packages used in the analysis:

```python
import session_info
```

Avoid renaming packages at import with `as` statements, unless you are performing a standard renaming (e.g., `import pandas as pd`).

Do not install new Python packages within a notebook.
All packages should instead be [installed and tracked using `conda`](../determining-requirements/determining-software-requirements.md#module-specific-conda-environments).


#### Setting paths

Defining paths to all input and output files at the start of the notebook makes it much easier for users to understand the analysis structure and to modify the paths if needed.

Do not use absolute file paths, as this can make it harder to share and reproduce the analysis.
Instead, define paths relative to the root of the OpenScPCA project or the root of the analysis module.

In the `hello.ipynb` notebook, we find the repository root path using the [`GitPython` package](https://gitpython.readthedocs.io/en/stable/intro.html) (imported with `import git`), as shown below:

```python
# find the repository root directory by looking for the .git directory
repo_root = git.Repo(".", search_parent_directories=True).working_dir
```

You can then use the [`os` package](https://docs.python.org/3/library/os.html) or the [`pathlib` package](https://docs.python.org/3/library/pathlib.html) to construct paths relative to the repository root, as shown below to find the module root:

```python
# use os.path.join to construct paths relative to the repository root
module_root = os.path.join(repo_root, "analyses", "hello-python")

# alternatively, use pathlib to construct paths
module_root_path = pathlib.Path(repo_root) / "analyses" / "hello-python"
```

Either of these packages can be used to construct paths in a platform-independent way, and both are part of the Python standard library, so you should use whichever you are more comfortable with.

### Defining custom functions

Following the Setup section, it is often useful to have a **Functions** section where you define any custom functions you write for use later in the notebook.
Keeping the functions in a central place in the notebook makes it easier to find and modify them later.

Each function should be defined in a separate code cell, and should be documented with comments to explain what the function does does and how to use it, including inputs and outputs.

For example, in `hello.ipynb`, we define a simple function to count cells in a `SingleCellExperiment` object:

```python
def count_anndata(anndata_file):
    """
    Count the number of cells in an anndata file.
    Returns a tuple with the Project ID, Sample ID, Library ID, and the number of cells.
    """
    ...
```

!!! note
    If you look at the `hello.ipynb` notebook itself, you will see that the function there has additional notation to take advantage of [Python's type hints](https://mypy.readthedocs.io/en/stable/cheat_sheet_py3.html).
    This can be a useful way to additionally document the expected inputs and outputs of a function, but it is not required.

If you have a large number of functions, you may consider moving them to a separate Python file and importing them into the notebook, but note that this makes the notebook less self-contained.
Imported function definitions will no longer appear in the output files, which can make it harder to share and reproduce the analysis.


### Analysis steps

Once all of the functions are defined, the remainder of the notebook should be a series of cells with code to perform the analysis steps, and [Markdown text](../../software-platforms/general-tools/writing-in-markdown.md) to explain what each step is doing and why.
Use headings and subheadings as appropriate to break up the analysis into logical sections, and include plots in the notebook to help illustrate the results.

Code cells should still contain comments to explain logic and implementation.
Markdown cells should focus on providing a higher-level overview of the analysis steps, including any interpretation of results.


### Session info

The final section of every Jupyter notebook should be a **Session info** section that uses the [`session_info.show()` function](../determining-requirements/determining-software-requirements.md#using-session_infoshow) to print out the versions of all packages used in the analysis.
