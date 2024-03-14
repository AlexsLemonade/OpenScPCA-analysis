# Software Requirements

## Determining and managing software dependencies in R

### Using `sessionInfo()`

The R programming language includes a helpful function called [`sessionInfo()`](https://stat.ethz.ch/R-manual/R-patched/library/utils/html/sessionInfo.html) that prints information about the R version and loaded packages.

#### In a notebook

Include `sessionInfo()` in the final chunk or cell of any computational notebook using R to ensure that the report output will include information about R and package versions.

#### In a script

Not using a notebook?
No problem.
You can capture the output of `sessionInfo()` in a text file by including the following in your script:

```r
sink("sessionInfo.txt")
sessionInfo()
sink()
```

This will result in the contents being written to a file called `sessionInfo.txt` in your working directory, so be sure to use unique filenames if multiple scripts in your analysis are writing out the contents of `sessionInfo()`.

You can post the contents of these files on [your pull requests](STUB_LINK index page for PRs) or commit them to the repository.

### Using `renv`

[`renv`](https://rstudio.github.io/renv/articles/renv.html) is a package manager for R.
We encourage you to read [the excellent introduction](https://rstudio.github.io/renv/articles/renv.html) for more information.
Here, we will limit our discussion to the most common commands we expect you to use when writing an analysis and assume that you have already [installed `renv`](STUB_LINK renv installation instructions).

#### Initializing `renv` in a module

!!! note
    If you [used `--use-renv` when creating your analysis module](STUB LINK), this step has been taken care of.

To start using `renv` in your analysis module, you can run the following R command from the root directory of your analysis module:

```r
renv::init()
```

[This command](https://rstudio.github.io/renv/reference/init.html) will create and selectively ignore a series of files, including the lockfile `renv.lock` that keeps track of the packages and versions being used.

#### Taking snapshots

As you develop your analysis, you may install packages using `install.packages()`, `renv::install()`, or `BiocManager::install()`.
You should periodically update the lockfile to make sure all dependencies are captured by using the following R command in the root directory of your module ([reference](https://rstudio.github.io/renv/reference/restore.html)):

```r
renv::snapshot()
```

When prompted, respond `y` to save the new packages in your lockfile.
Commit the changes to your lockfile to the repository.

You can use [`renv::status()`](https://rstudio.github.io/renv/reference/status.html) at any time to check if there are inconsistencies between the module dependencies and the lockfile.

#### Pinning dependencies that are not captured automatically

Taking [a snapshot using the default arguments will only capture packages that are used in your module and their required dependencies](https://rstudio.github.io/renv/reference/snapshot.html), but there may be some other _recommended_ package that you want to include and pin to a specific version.

For example, `ggplot2` needs the `svglite` package to save `.svg` files, but that package is not listed as a _requirement_, so `renv` may not know to track it, even if you have the package installed. 

You can make `renv` include a package by loading it in a file called `dependencies.R` in a directory called `components` within your analysis.

!!! note
    If you [used `--use-renv` when creating your analysis module](STUB LINK), `components/dependencies.R` was already created.
    
For instance, if you wanted to make sure `renv` was keeping track of the `scuttle` and `scRNAseq` packages, your module's `components/dependencies.R` would include the following:

```r
library(scuttle)
library(scRNAseq)
```

#### Restoring dependencies

If you are switching between your computer and Amazon Web Services or running someone else's module, it can be helpful to restore dependencies using a module's lockfile.
To do so, use the following command ([reference](https://rstudio.github.io/renv/reference/restore.html)):

```r
renv::restore()
```
