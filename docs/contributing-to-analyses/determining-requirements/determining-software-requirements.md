# Software requirements

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

You can post the contents of these files on [your pull requests](../creating-pull-requests/index.md) or commit them to the repository.

### Using `renv`

[`renv`](https://rstudio.github.io/renv/) is a package manager for R.
We encourage you to read [the excellent introduction](https://rstudio.github.io/renv/articles/renv.html) for more information.
Here, we will limit our discussion to the most common commands we expect you to use when writing an analysis and assume that you have already [installed `renv`](../../technical-setup/environment-setup/install-r-rstudio.md#install-r-package-dependencies).

#### Initializing `renv` in a module

!!! note
    If you [used `--use-renv` when creating your analysis module](../analysis-modules/creating-a-module.md#use-renv), this step has been taken care of.

To start using `renv` in your analysis module, you can run the following R command from the root directory of your analysis module:

```r
renv::init()
```

[This command](https://rstudio.github.io/renv/reference/init.html) will create and selectively ignore a series of files, including the lockfile `renv.lock` that keeps track of the packages and versions being used.

You can commit any changes introduced by running `renv::init()` to the repository.

The `.Rprofile` file created by `renv::init()` will cause `renv` to automatically activate its library if you launch `R` or `Rscript` from your module's root directory.

#### Using `renv` with RStudio Projects

If you are using [the RStudio IDE](../../technical-setup/environment-setup/install-r-rstudio.md#install-the-rstudio-ide) for developing your module, you might find it helpful to use [an RStudio Project](https://docs.posit.co/ide/user/ide/guide/code/projects.html) for your module in conjunction with `renv`.

Following the linked directions, we recommend creating an RStudio Project from an existing directory (your module).
You can commit the `.Rproj` file that is created to the repository.

When you open the Project in RStudio, it will automatically source the `.Rprofile` file created by `renv::init()` and activate `renv` in your R session.

#### Taking snapshots

As you develop your analysis, you may install packages using `install.packages()`, `renv::install()`, or `BiocManager::install()`.
You should periodically update the lockfile to make sure all dependencies are captured by using the following R command in the root directory of your module ([reference](https://rstudio.github.io/renv/reference/snapshot.html)):

```r
renv::snapshot()
```

When prompted, respond `y` to save the new packages in your `renv.lock` file.
Commit the changes to your lockfile to the repository.

You can use [`renv::status()`](https://rstudio.github.io/renv/reference/status.html) at any time to check if there are inconsistencies between the module dependencies and the lockfile.

#### Pinning dependencies that are not captured automatically

Taking [a snapshot using the default arguments will only capture packages that are used in your module and their required dependencies](https://rstudio.github.io/renv/reference/snapshot.html), but there may be some other _recommended_ package that you want to include and pin to a specific version.

For example, `ggplot2` needs the `svglite` package to save `.svg` files, but that package is not listed as a _requirement_, so `renv` may not know to track it, even if you have the package installed.

You can make `renv` include a package by loading it in a file called `dependencies.R` in a directory called `components` within your analysis.

!!! note
    If you [used `--use-renv` when creating your analysis module](../analysis-modules/creating-a-module.md#use-renv), `components/dependencies.R` was already created.

For instance, if you wanted to make sure `renv` was keeping track of the `scuttle` and `svglite` packages, your module's `components/dependencies.R` would include the following:

```r
library(scuttle)
library(svglite)
```

#### Restoring dependencies

If you are switching between your computer and Amazon Web Services or running another contributor's module, it can be helpful to restore dependencies using a module's lockfile.
To do so, use [`renv::restore()` function](https://rstudio.github.io/renv/reference/restore.html):

```r
renv::restore()
```

## Managing software dependencies with conda

We strongly recommend using [conda](https://docs.conda.io/en/latest/) and [`conda-lock`](https://conda.github.io/conda-lock/) to manage dependencies for any module written primarily in Python, as well as for standalone software packages that do not depend on a specific language.
Many useful bioinformatics tools are available through [Bioconda](https://bioconda.github.io/), which should already be configured as a conda channel.
These instructions assume you have already [installed conda and set up your base environment](../../technical-setup/environment-setup/setup-conda.md).

### Module-specific conda environments

#### Creating environments

When [creating a Python module using `create-analysis-module.py`](../analysis-modules/creating-a-module.md#flags-to-create-a-python-module), a conda environment will be created for that module by default with two effects:

* A basic `environment.yml` file will already be present in the module's directory
* A conda environment named `openscpca-{module_name}` will already exist

You can activate the environment by running the following command:

```bash
conda activate openscpca-{module_name}
```

If no `environment.yml` file is present, you can create and activate a new environment for the module by running the following command from the module's root directory, replacing `{module_name}` with the name of the module you are working on:

```bash
conda create --file environment.yml --name openscpca-{module_name}
conda activate openscpca-{module_name}
```

#### Freezing dependencies with `conda-lock`

Conda is a powerful package manager, but it has some limitations when trying to perfectly reproduce an environment across different computer platforms (e.g., between Linux and macOS).
In particular, some packages have dependencies that vary depending on the platform that is being used, which makes tracking dependencies across platforms difficult with the built-in conda tools.
To address this, we use [`conda-lock`](https://conda.github.io/conda-lock/), which provides the following advantages:

- The `conda-lock.yml` files created by `conda-lock` capture exact versions of all software packages in the environment, including both packages listed in the `environment.yml` file and all of their dependencies.
- Dependencies required for each platform are discovered and recorded separately by `conda-lock`, ensuring that the environment can be recreated on every platform where the software is available.
- This means that we can leave the `environment.yml` file as a "high-level" list of required software which is easy to read and update, while also having a "low-level" lockfile that contributors can use to recreate the environment exactly.


#### Creating and updating `conda-lock.yml` files

To create a `conda-lock.yml` file from an `environment.yml` file, run the following command from the module directory:

```bash
conda-lock --file environment.yml
```

This will create a `conda-lock.yml` file in the module directory that contains the exact versions of all software packages in the environment, including any dependencies that may be specific to a given platform.

You should perform this step before [filing a pull request](../creating-pull-requests/index.md) to ensure that the `conda-lock.yml` file is up-to-date with the current state of your environment.

!!! note
    If the `conda-lock` command fails, it may be because a package is not available for one of the platforms listed in the `environment.yml` file.
    Usually this will be a package that is not available for the `osx-arm64` (Apple Silicon) platform.

    If this happens, see the [Software not available on a specific platform](#software-not-available-on-a-specific-platform) section below for instructions on how to handle this situation.

#### Activating existing environments in a module

If you are working with an existing analysis module, it should have a `conda-lock.yml` file present in the module directory.
You can create and activate the environment on the computer you are using by running the following commands:

```bash
# Navigate to the module's root directory
cd analyses/{module_name}/
conda-lock install --name openscpca-{module_name} conda-lock.yml
conda activate openscpca-{module_name}
```

If the `conda-lock.yml` file has been updated since the last time you worked with the repository, you should rerun the `conda-lock install` command to update your environment to the latest version.

Occasionally you may encounter a Python-based module that does not yet have a `conda-lock.yml` file.
This is most likely to occur early in development of a module when the requirements may still be in flux.
(If you do find a module that is missing a `conda-lock.yml` file, you may want to [file an issue](../../communications-tools/github-issues/index.md) to let us know.)
As long as there is an `environment.yml` file, you can still create an environment to work with the module.
To do so, run the following commands from the module's root directory, replacing `{module_name}` with the name of the module you are working on:

```bash
# Navigate to the module's root directory
cd analyses/{module_name}/
conda env create --name openscpca-{module_name} --file environment.yml
conda activate openscpca-{module_name}
```

!!! note
    During the `conda env create` step, you may see the following warning:

    ```{.console .no-copy}
    EnvironmentSectionNotValid: The following section on 'environment.yml' is invalid and will be ignored:
     - platforms
    ```

    This warning can be safely ignored; proceed with the installation and activate the environment.


### Adding software to the environment and tracking installed software

If there is additional software you intend to use or discover is required for a module, you can add it with the `conda install` command when the environment is activated.
Follow installing with an update to the `environment.yml` file with that package and the dependencies that were installed.

For example, to install the `pandas` package and record it in the `environment.yml` file, you would run the following commands:

```bash
# Navigate to the module's root directory
cd analyses/{module_name}
# Activate that module's environment
conda activate openscpca-{module_name}
# Install the package
conda install pandas
```

You should then immediately add the newly installed package to the module's `environment.yml` file.

To do this, first check which version was installed, which you can find by either:

- Looking at the output from the previous `conda install` command, or
- Running a command like the following, which will print out the version number of the package in your current environment:

    ```bash
    conda list pandas
    ```

Then add the package to the `dependencies:` section of the `environment.yml` file with a version number, as shown below:

<div class="grid" markdown>

``` { .yaml .no-copy title="Before"}
dependencies:
  - python=3.11
  - awscli=2.15
  - conda-lock=2.5
  - jq=1.7
  - pre-commit=3.6
  - session-info=1.0
```

``` { .yaml .no-copy title="After"}
dependencies:
  - python=3.11
  - awscli=2.15
  - conda-lock=2.5
  - jq=1.7
  - pre-commit=3.6
  - session-info=1.0
  - pandas=2.2.1
```

</div>

#### Finding available software

To find software available through conda, you can search the conda repositories with the following command:

```bash
conda search {package_name}
```

Alternatively, you can search [anaconda.org](https://anaconda.org) for packages and channels.


### Software not available on a specific platform

While most conda packages are available for all platforms, there may be some cases where a particular platform does not have a version of a package.

Most often this will occur for ARM-based computers, such as macOS computers with Apple Silicon M-series processors.
If you encounter an error with `conda-lock --file environment.yml`, it may be because a package is not available for the `osx-arm64` platform.

In this case, you should edit the `environment.yml` file to *remove* the `- osx-arm64` line from the `platforms:` section.
Then you will want to rerun `conda-lock`, but this time creating platform-specific lockfiles:

```bash
conda-lock --file environment.yml --kind explicit
```

This will result in separate files for each supported platform with names like `conda-linux-64.lock` or `conda-osx-64.lock`.

You can then use `conda-lock` to install or update your environment from the lockfile specific to the platform that you are using.
For example, to install and activate and environment from the `conda-linux-64.lock` file, you would run the following commands:

```bash
conda-lock install --name openscpca-{module_name} conda-linux-64.lock
conda activate openscpca-{module_name}
```

If you are using an Apple Silicon computer, you will need to use the `conda-osx-64.lock` file, with one extra option to allow installation of the non-native platform lockfile, and one extra command to ensure that any future software you install in the environment will use the same architecture:

```bash
conda-lock install --no-validate-platform --name openscpca-{module_name} conda-osx-64.lock
conda activate openscpca-{module_name}
conda config --env --set subdir osx-64
```

!!! tip "Adding packages that are not available for your current platform"
    If you find you need to add a package to an existing environment that is not available for your _current_ platform (i.e., a package that is only available for Intel on an Apple Silicon computer), the easiest method is to follow a slightly different order of operations:

    - First, add the package name and version to the `environment.yml` file.
    - Next, remove the `osx-arm64` line from the `platforms:` section of the `environment.yml` file.
    - Then, run `conda-lock --file environment.yml --kind explicit` to create the platform specific lockfiles.
    - Finally, use the `conda-lock install --no-validate-platform` command above to update your environment using the newly-created lockfile.


### Using `session_info.show()` in Python

The [`session-info`](https://pypi.org/project/session-info/) Python package can be used to report version information about Python and loaded modules.
If you [created a module using `--use-jupyter` or `--use-python`](../analysis-modules/creating-a-module.md#use-jupyter), `session-info` was automatically included in the module's conda environment.

Import `session-info` by placing the following in [the Setup section of your Jupyter notebook](../analysis-modules/notebook-structure.md#jupyter-notebooks) or in the [load packages section of your script](../analysis-modules/script-structure.md):

```python
import session_info
```

#### In Jupyter notebooks

Include the following in the final cell of any Jupyter notebook using Python to ensure that the notebook reports version information:

```python
session_info.show()
```

#### In Python scripts

Using a script instead?

Import `contextlib` by placing the following in the [load packages section of your script](../analysis-modules/script-structure.md):

```python
import contextlib
```

Then, you can use the following to write the output of `session_info.show()` to a text file called `session_info.txt` in your current directory:

```python
session_info_path = "session_info.txt"

with open(session_into_path, "w") as f:
    with contextlib.redirect_stdout(f):  # direct the session_info output to a file
        session_info.show(dependencies=True)
```

If you're using this approach with multiple scripts in your module, be sure to use unique file names when naming your session info text files.

You can post the contents of these files on [your pull requests](../creating-pull-requests/index.md) or commit them to the repository.
