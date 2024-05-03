# Install and set up conda

## What is conda?

The OpenScPCA project uses [conda](https://docs.anaconda.com/free/miniconda/index.html) to manage your software environment.
Conda is a command-line software management tool which helps you install and track specific versions of software.
It also allows you to have multiple software environments with different sets of packages on the same computer.

This page provides instructions on how to install conda and use it to install certain software you will need to contribute to OpenScPCA.


### Why use conda?

There are two main reasons we use conda for OpenScPCA:

- conda provides a "one stop shop" for installing lots of different software
    - Rather than having to figure out how to install every new program and its dependencies, conda can handle it all for you.
    - You'll use conda to install the software dependencies you'll need to contribute to OpenScPCA.
- conda allows you set up different software environments for different projects
    - For example, you may have two projects that require different versions of the same package.
    With conda, you can create separate, fully isolated software environments for each project with different package versions.
    - Python-based OpenScPCA analysis modules as well as those that require other external software will use different conda environments to prevent conflicts and improve reproducibility.
    Therefore, installing conda is also part of setting up your computer to be able to contribute to Python-based modules.


## Install conda

We recommend installing [Miniconda](https://docs.anaconda.com/free/miniconda/index.html) to obtain conda.
Miniconda is lightweight version of the full conda platform and includes the conda tool itself, Python, and a few other commonly-used packages.

If you already have conda on your system, you do not need to re-install it.


### Install Miniconda on macOS

To install Miniconda, [download the graphical installer for macOS](https://docs.anaconda.com/free/miniconda/miniconda-install/), and follow all instructions.

  - If you are on a macOS computer, be sure to download one of the links ending in `pkg`, _not `bash`_:
    - Apple Silicon (M1-3) Mac users should download `Miniconda3 macOS Apple M1 64-bit pkg`
    - Intel Mac users should download `Miniconda3 macOS Intel x86 64-bit pkg`

### Install Miniconda on Windows with WSL 2

To install Miniconda onto the WSL 2 side of your Windows machine, you will need to follow the [Quick command line install instructions for Linux](https://docs.anaconda.com/free/miniconda/#quick-command-line-install).

- Copy and paste all commands, in order, into the [Ubuntu terminal](../../software-platforms/general-tools/using-the-terminal.md).
  - Be sure to also take the step to initialize conda for the Bash shell.
- Then, as prompted, close your terminal and open a new terminal window to complete installation.

## Set up conda

Next, you will need to set certain conda settings and install a few packages that will allow you to contribute to OpenScPCA in general.

1. [Open a terminal (command line prompt)](../../software-platforms/general-tools/using-the-terminal.md) to interact with conda.
    - We recommend opening a [terminal from GitKraken](../../software-platforms/general-tools/using-the-terminal.md#gitkraken) since you'll need to run some of these steps from the `OpenScPCA-analysis` folder.

1. Copy and paste the following code into the terminal, and hit enter.
These commands will set the [recommended channels](https://docs.conda.io/projects/conda/en/latest/user-guide/concepts/channels.html) conda should use to search for software.
    - If you get the error `conda: command not found`, you may need to try again in a new terminal window.
    If this still doesn't work, [you can always ask us for help](../../troubleshooting-faq/index.md).

    ```sh
    conda config --add channels defaults
    conda config --add channels bioconda
    conda config --add channels conda-forge
    conda config --set channel_priority strict
    ```


### Create an `openscpca` conda environment

The last step is to create an `openscpca` conda environment and install the packages needed for OpenScPCA development.
These are specified in the `environment.yml` in the root of the repository, and include:

- The [`awscli` package](https://docs.aws.amazon.com/cli/latest/userguide/cli-chap-welcome.html) will allow you to interact with [data stored in the Amazon Web Services (AWS) S3 bucket](../../software-platforms/aws/index.md)
- The [`conda-lock` package](https://conda.github.io/conda-lock/) is used to create fully reproducible, cross-platform [conda environments for analysis modules](../../contributing-to-analyses/determining-requirements/determining-software-requirements.md#managing-software-dependencies-with-conda)
- The [`jq` package](https://jqlang.github.io/jq/) provides JSON parsing capabilities
- The [`pre-commit` package](https://pre-commit.com) will allow you to use [pre-commit hooks when contributing to analysis modules](../../contributing-to-analyses/working-with-git/making-commits.md#pre-commit-checks)

<!-- Comment to force above to be bullets, next to be numbers -->


1. To create this environment, navigate in the terminal with `cd` to the `OpenScPCA-analysis` repository.
    - Again, if you open a [terminal in GitKraken](../../software-platforms/general-tools/using-the-terminal.md#gitkraken), you will automatically be in the repository folder.

1. Enter the following command in the terminal to create the `openscpca` environment and install the packages into it:

    ```sh
    conda env create -f environment.yml
    ```

    You will see a number of status messages flow by as packages are downloaded and installed.
    When the process is complete, you'll see the following message in the terminal, which includes the commands to activate and deactivate the environment:

    ```{ .console .no-copy title="Output message after conda environment install"}
    Preparing transaction: done
    Verifying transaction: done
    Executing transaction: done
    # To activate this environment, use
    #
    #     $ conda activate openscpca
    #
    # To deactivate an active environment, use
    #
    #     $ conda deactivate
    ```

    ??? info "Prefer to modify your base environment instead?"
        While we recommended creating an `openscpca` environment, this is not strictly necessary.
        If you prefer to add these packages to your base environment, use this command instead:
        ```sh
        conda env update -n base -f environment.yml
        ```
        If you experience package conflicts when running this command, we recommend you go back to setting up a dedicated `openscpca` environment.



1. As indicated in the message, activate the environment with the following command:

    ```sh
    conda activate openscpca
    ```

    At this point, your terminal prompt will show the prefix `(openscpca)`, which lets you know that you are in that conda environment.
    In general, this is how conda indicates which environment, if any, the terminal is working in.

    !!! tip "Activating the `openscpca` conda environment"
        Whenever you are developing for OpenScPCA, you should work from either the `openscpa` conda environment or [a module-specific conda environment](../../contributing-to-analyses/determining-requirements/determining-software-requirements.md#managing-software-dependencies-with-conda).
        At the start of any coding session, run `conda activate openscpca` to activate this environment.
        If you wish to deactivate the environment, you can run `conda deactivate`.
