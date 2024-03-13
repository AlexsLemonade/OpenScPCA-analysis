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

To start using `renv` in your analysis module, you can run the following from the root directory of your analysis module:

```r
renv::init()
```

[Which will create and selectively ignore a series of files, including the lockfile `renv.lock` that keeps track of the packages and versions being used.](https://rstudio.github.io/renv/reference/init.html)

#### Taking snapshots

As you develop your analysis, you may install packages using `install.packages()` or `renv::install()`.
You should periodically update the lockfile to make sure all dependencies are captured by using the following command in the root directory of your module ([reference](https://rstudio.github.io/renv/reference/restore.html)):

```r
renv::snapshot()
```

Commit the changes to your lockfile.

You can use [`renv::status()`](https://rstudio.github.io/renv/reference/status.html) to check if there are inconsistencies between the module dependencies and the lockfile.

#### Pinning dependencies that are not captured automatically

Taking [a snapshot using the default arguments will only capture packages that are used in your module](https://rstudio.github.io/renv/reference/snapshot.html), but there may be some other recommended package that you want to pin to a specific version.

You can make `renv` include a packages by loading it in a file called `dependencies.R` in a directory called `components` within your analysis.

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

## Managing software dependencies in Python with conda

We strongly recommend using conda to manage dependencies for any module written primarily in Python.
These instructions assume you have already [installed conda and set up your base environment](STUB_LINK conda technical setup).

### Module-specific conda environments

#### Creating environments

When [creating a Python module using `create-analysis-module.py`](STUB_LINK module creation), a conda environment will be created for that module by default with two effects:

* A basic `environment.yaml` file will already be present in the module's directory
* A conda environment named `openscpca-{module_name}` will already exist

You can activate the environment by running the following command:

```bash
conda activate openscpca-{module_name}
```

If no `environment.yaml` file is present, you can create and activate a new environment for the module by running the following command from the module's root directory, replacing `{module_name}` with the name of the module you are working on:

```bash
conda create --file environment.yml --name openscpca-{module_name}
conda activate openscpca-{module_name}
```

#### Activating existing environments

If you are working with an existing module that already has an  `environment.yaml` file or have switched computers, you can create and activate the environment on the computer you are working on by running the following commands:

```bash
# Navigate to the module's root directory
cd analyses/{module_name}/
conda env create --file environment.yml --name openscpca-{module_name}
conda activate openscpca-{module_name}
```

#### Note for ARM (Apple Silicon) computers

While most conda packages are available for ARM-based computers (such as macOS computers with M-series processors), some software is only available for Intel architectures.
However, it is still usually possible to run an Intel-based package on macOS, but you will need to create an environment that uses the `osx-64` architecture for all of its software.
To do this, modify the creation command for the environment and add a setting to the environment, as shown below:

```bash
cd analyses/{module_name}/
CONDA_SUBDIR=osx-64 conda env create --file environment.yaml --name openscpca-{module_name}
conda activate openscpca-{module_name}
conda config --env --set subdir osx-64
```

Moving forward, you should be able to install any Intel-based package into the environment as usual.

### Adding software to the environment and tracking installed software

If there is additional software you intend to use for a module or discover is required, you can add it with the `conda install` command when the environment is activated.
Follow installing with an update to the `environment.yml` file with that package and the dependencies that were installed.

For example, to install the `pandas` package and record it in the `environment.yml` file, you would run the following commands:

```bash
# Navigate to the module's root directory
cd analyses/{module_name}
# Activate that module's envirionment
conda activate openscpca-{module_name}
# Install the package
conda install pandas
# Update environment.yml
conda env export --no-builds | grep -v "^prefix:" > environment.yml
```

(The `grep` command in the final line is there to remove user-specific paths that `conda` includes in its export.)

### Finding available software

To find software available through conda, you can search the conda repositories with the following command:

```bash
conda search {package_name}
```

Alternatively, you can search [anaconda.org](https://anaconda.org) for packages and channels.
