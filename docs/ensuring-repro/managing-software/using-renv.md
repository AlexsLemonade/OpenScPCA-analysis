# Managing dependencies with renv

[`renv`](https://rstudio.github.io/renv/) is a package manager for R.
We encourage you to read [the excellent introduction](https://rstudio.github.io/renv/articles/renv.html) for more information.
Here, we will limit our discussion of `renv` to the most common commands we expect you to use when writing an analysis and assume that you have already [installed `renv`](../../technical-setup/environment-setup/install-r-rstudio.md#install-r-package-dependencies).

## Initializing a module-specific `renv` environment

When creating an R module using [`create-analysis-module.py` with the `--use-renv` flag](../../contributing-to-analyses/analysis-modules/creating-a-module.md#flags-to-create-an-r-module), a `renv` environment will be initialized in your module directory.

If you did not use the `--use-renv` flag but need to add an `renv` environment, you can run the following R command from the root directory of your analysis module:

```r
renv::init()
```

[The `renv::init()` command](https://rstudio.github.io/renv/reference/init.html) will create and selectively ignore a series of files, including the lockfile `renv.lock` that keeps track of the packages and versions being used.

You can commit any changes introduced by adding and committing `renv::init()` to the repository.

The `.Rprofile` file created by `renv::init()` will cause `renv` to automatically activate its library if you launch `R` or `Rscript` from your module's root directory.
As long as you run a module's code from its root directory, the `renv` environment will be active.

## Using `renv` with RStudio Projects

If you are using [the RStudio IDE](../../technical-setup/environment-setup/install-r-rstudio.md#install-the-rstudio-ide) for developing your module, you might find it helpful to use [an RStudio Project](https://docs.posit.co/ide/user/ide/guide/code/projects.html) for your module in conjunction with `renv`.

Following the linked directions, we recommend creating an RStudio Project from an existing directory (your module).
You can commit the `.Rproj` file that is created to the repository.

When you open the Project in RStudio, it will automatically source the `.Rprofile` file created by `renv::init()` and activate `renv` in your R session.

## Adding packages to the environment

### Taking snapshots

As you develop your analysis, you may install packages using `install.packages()`, `renv::install()`, or `BiocManager::install()`.
You should periodically update the lockfile to make sure all dependencies are captured by using the following R command in the root directory of your module ([reference](https://rstudio.github.io/renv/reference/snapshot.html)):

```r
renv::snapshot()
```

When prompted, respond `y` to save the new packages in your `renv.lock` file.
Commit the changes to your lockfile to the repository.

You can use [`renv::status()`](https://rstudio.github.io/renv/reference/status.html) at any time to check if there are inconsistencies between the module dependencies and the lockfile.

### Pinning dependencies that are not captured automatically

Taking [a snapshot using the default arguments will only capture packages that are used in your module and their required dependencies](https://rstudio.github.io/renv/reference/snapshot.html), but there may be some other _recommended_ package that you want to include and pin to a specific version.

For example, `ggplot2` needs the `svglite` package to save `.svg` files, but that package is not listed as a _requirement_, so `renv` may not know to track it, even if you have the package installed.

You can make `renv` include a package by loading it in a file called `dependencies.R` in a directory called `components` within your analysis.

!!! note
    If you [used `--use-renv` when creating your analysis module](../../contributing-to-analyses/analysis-modules/creating-a-module.md#use-renv), `components/dependencies.R` was already created.

For instance, if you wanted to make sure `renv` was keeping track of the `scuttle` and `svglite` packages, your module's `components/dependencies.R` would include the following:

```r
library(scuttle)
library(svglite)
```

## Restoring dependencies

If you are switching between your computer and Amazon Web Services or running another contributor's module, it can be helpful to restore dependencies using a module's lockfile.
To do so, use [`renv::restore()` function](https://rstudio.github.io/renv/reference/restore.html):

```r
renv::restore()
```
