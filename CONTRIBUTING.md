# Contributing to OpenScPCA

**Table of Contents**
- [Conda environment setup](#conda-environment-setup)
  - [Installing Miniconda](#installing-miniconda)
  - [Creating and activating the base OpenScPCA environment](#creating-and-activating-the-base-openscpca-environment)
  - [Module-specific environments](#module-specific-environments)
  - [Adding software to the environment and tracking installed software](#adding-software-to-the-environment-and-tracking-installed-software)
  - [Updating `conda-lock.yml`](#updating-conda-lockyml)
  - [Finding available software](#finding-available-software)
    - [Note for ARM (Apple Silicon) computers](#note-for-arm-apple-silicon-computers)
- [Setting up pre-commit](#setting-up-pre-commit)
  - [Adding additional local hooks](#adding-additional-local-hooks)
    - [Code formatting and linting](#code-formatting-and-linting)
    - [Spell checking](#spell-checking)
    - [Other pre-commit hooks](#other-pre-commit-hooks)
- [Creating a new analysis module](#creating-a-new-analysis-module)

## Conda environment setup

To facilitate software setup and reproducibility, we have provided a basic conda environment file, `environment.yml`, that you can use to create a virtual environment with all the necessary dependencies.
The software included in this file includes:

- Python 3.10
- `pre-commit` for managing code quality checks
- `aws` command line tool for interacting with AWS
- `jq` for parsing JSON files (useful with `aws`)

### Installing Miniconda

If you do not already have a working conda installation, we encourage you to install Miniconda, which includes only the conda package manager and its dependencies, rather than the full Anaconda distribution.
This will save disk space and make it easier to manage your conda installation.

Install Miniconda by following the instructions in the [Miniconda documentation](https://docs.anaconda.com/free/miniconda/#quick-command-line-install).
Note that the installation instructions differ by operating system and architecture, so be sure to select the correct installation instructions for your system.

Once you have completed the basic installation, define the default conda channels and priorities by running the following commands in your terminal:

```bash
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority strict
```


### Creating and activating the base OpenScPCA environment

Once Miniconda is installed on your system, you can update your base environment to include the packages that you will need for OpenScPCA development.
These packages are listed in the repository `environment.yml` file, and you can install them into your current environment by running the following command in the root directory of the repository:

```bash
conda env update --file environment.yml
```

### Module-specific environments

When you are working on an analysis module, you can create and/or use a conda environment specific to that module.
If you used the `create-analysis-module.py` script to create the module with the `--use-conda` option, a basic `environment.yaml` file should already be present in the module's directory and a conda environment named `openscpca-{module_name}` should already exist.
See the [Creating a new analysis module](#creating-a-new-analysis-module) section for more information on creating a new module.
You can activate the environment by running the following command:

```bash
conda activate openscpca-{module_name}
```

If no `environment.yaml` file is present, you can create and activate a new environment for the module by running the following command from the repository root directory, replacing `{module_name}` with the name of the module you are working on:

```bash
conda create --file environment.yml --name openscpca-{module_name}
conda activate openscpca-{module_name}
```

If you are working with a pre-existing module that already has an  `conda-lock.yml` file, you can create and activate the environment from that file by running the following commands:

```bash
conda-lock install --name openscpca-{module_name} analyses/{module_name}/conda-lock.yml
conda activate openscpca-{module_name}
```

### Adding software to the environment and tracking installed software

If there is additional software you intend to use for a module, you can add it with the `conda install` command when the environment is activated.
When you install new software into your environment, you will also want to update the `environment.yml` file with that package and the dependencies that were installed.
For example, to install the `pandas` package and record it in the `environment.yml` file, you would run the following commands:

```bash
cd analyses/{module_name}
conda activate openscpca-{module_name}
conda install pandas
```

You should then add the newly installed module to the module's `environment.yml` file.
First check which version was installed, either by looking at the output from the install command, or by running a command like the following:

```bash
conda list pandas
```

Then add the module to the `dependencies` section of the `environment.yml` file with a version number, as shown below:

```yaml
  - pandas=2.2.1
```


### Updating `conda-lock.yml`

Before submitting a pull request, you should also update the `conda-lock.yml` file in the module directory.

To do this, run the following command from the module directory:

```bash
conda-lock --file environment.yml
```

To be sure that this exactly matches your environment, it is a good idea to follow locking the environment with:

```bash
conda-lock install --name openscpca-{module_name} conda-lock.yml
conda activate openscpca-{module_name}
```



### Finding available software

To find software available through conda, you can search the conda repositories with the following command:

```bash
conda search {package_name}
```

Alternatively, you can search [anaconda.org](https://anaconda.org) for packages and channels.

#### Note for ARM (Apple Silicon) computers

While most conda packages are available for ARM-based computers (such as macOS computers with M-series processors), some software is only available for Intel architectures.
However, it is still usually possible to run an Intel-based package on macOS, but you will need to create an environment that uses the `osx-64` architecture for all of its software.
To do this, you can slightly modify the creation command for the environment and add a setting to the environment, as shown below:

```bash
CONDA_SUBDIR=osx-64 conda env create --file environment.yaml --name {module_name}
conda activate {module_name}
conda config --env --set subdir osx-64
```

After that point, you should be able to install any Intel-based package into the environment as usual.

## Setting up pre-commit

[`pre-commit`](https://pre-commit.com) is a small software package that makes it easy to manage and run code quality checks.
All contributors should use pre-commit as part of their workflow, installing the package as described below.
`pre-commit` checks code quality by defining a set of "hooks" that will run every time you commit changes to a repository.
We have used it in this project to set up some pre-commit hooks to manage basic code security and other common errors, such as the following:

- Large data files that should not be committed to the repository
- Credential files and other sensitive information
- Merge conflicts that have not yet been resolved

Pre-commit is installed via the conda environment file, so you should already have it installed if you have followed the instructions above.
With pre-commit installed, you can activate the hooks by running the following command in the root directory of the repository:

```bash
# run in the root directory of the repository
pre-commit install
```

After that point, any time you commit a change to the repository, the hooks will run and check for errors.
If any errors are found, the commit will be aborted and you will be prompted to fix the errors, after which you can retry the commit.
For some hooks, the errors will be automatically fixed, and you will only need to stage the updated files and retry the commit.

Note that the first time you commit after installing `pre-commit` or updating the hooks, the hooks may take a while to complete, as software may need to be downloaded and installed, but the hooks will run much faster on subsequent commits.

### Adding additional local hooks

While we have taken a limited approach to the required pre-commit hooks in this project, there are a number of other pre-commit hooks that you might find useful for your own development.

Because the `.pre-commit-config.yaml` file is shared by all repository users, it should not be modified.
Instead, you should make a copy of this file in the root directory of the repository named `.pre-commit-local.yaml`.
We have configured Git to ignore this file in the repository `.gitignore` file, so you are free to modify it as you wish without affecting other users.
To switch your local activation of `pre-commit` to use this file, you will need to re-activate `pre-commit` using the the following code:

```bash
# make and activate a local pre-commit configuration
cp .pre-commit-config.yaml .pre-commit-local.yaml
pre-commit install --config .pre-commit-local.yaml
```

You can then add your own pre-commit hooks by editing the `.pre-commit-local.yaml` file.

#### Code formatting and linting

One example is code formatting and linting tools, which can help you write error-free, consistent, and readable code.
For more on the value of these tools, see [this article about linters and formatters](https://www.freecodecamp.org/news/using-prettier-and-jslint/).
While the article focuses on JavaScript, the same principles apply to other languages.

We have not included any code formatting or linting checks in the pre-commit hooks we require, but you may find it helpful to use these tools in your own copy of the repository.

Note that these tools will often directly modify your files when run.
If they are run as a pre-commit hook the initial commit will fail, and you will then need to check and stage the changes that were made by the tool before re-trying the commit.


Some formatters that we recommend are [`ruff-format`](https://docs.astral.sh/ruff/formatter/) for Python and the [`style-files` hook from the `precommit` package](https://lorenzwalthert.github.io/precommit/articles/available-hooks.html#style-files) for R.
You can add those with the following code added to the `.pre-commit-config.yaml` file in the `repos:` section:

```yaml
  # ruff formatter for Python
  - repo: https://github.com/astral-sh/ruff-pre-commit
    rev: v0.2.1
    hooks:
      - id: ruff-format
  # code styling with the {styler} package for R
  - repo: https://github.com/lorenzwalthert/precommit
    rev: v0.4.0
    hooks:
      - id: style-files
  # prettier formatter for many other languages
  - repo: https://github.com/pre-commit/mirrors-prettier
    rev: v3.1.0
    hooks:
      - id: prettier
```

Code linting tools are often more intrusive, enforcing not only general formatting but also particular style standards and "best" practices.
This can make them more likely to find errors and inconsistencies, but also more likely to require manual intervention to fix those errors.
They also might complain about things that are not actually errors, but are simply not to the linter's taste.
For python, we recommend [`ruff`](https://docs.astral.sh/ruff/) (which goes along with `ruff-format`, above), and for R we recommend the [`lintr`](https://lintr.r-lib.org) package.


#### Spell checking

Spell checking is another useful class of tools that can be added as a pre-commit hook.
Because biological and computational words are often not in default dictionaries, it may be more helpful in a pre-commit context to use a spell checker that looks for common errors rather than one that checks every word against a dictionary.
One such tool is [`typos`](https://github.com/crate-ci/typos), which runs quickly, checking both text and code for common mistakes.
 `typos` can be installed as a pre-commit hook with the following code added to the `.pre-commit-config.yaml` file:

```yaml
  - repo: https://github.com/crate-ci/typos
    rev: v1.18.2
    hooks:
      - id: typos
```

#### Other pre-commit hooks

There are many other pre-commit hooks available, and you can find a fairly extensive list of them at the pre-commit ["Supported hooks" page](https://pre-commit.com/hooks.html).
Just note that the more hooks you add, the longer each commit to your repository will take!


## Creating a new analysis module

To facilitate creating new analysis modules with recommended file structure and documentation, we have provided a script, `create-analysis-module.py`, that you can use to create a new module.
This script will create a new directory within the `analyses` directory for the new module and populate it with a recommended file structure and some template files.

If you plan to use `conda` or `renv` for package management, you can also use the script to create a conda environment file or initialize an `renv` project.
To use the script, run one of the following commands in the root directory of the repository:

```bash
# a basic module with no environment setup
./create-analysis-module.py {module_name}

# a module with a conda environment
./create-analysis-module.py {module_name} --use-conda

# a module with an renv environment
./create-analysis-module.py {module_name} --use-renv
```

If you have created a module with a conda environment, the conda environment created will be named `openscpca-{module_name} and you can then activate it with the following command:

```bash
conda activate openscpca-{module_name}
```

See the [Conda environment setup](#conda-environment-setup) section for more information on managing conda environments.

If you have created a module with an `renv` environment, that environment will automatically be activated when you open `R` in the module directory.
