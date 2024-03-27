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
    - Python-based OpenScPCA analysis modules will use different conda environments to prevent conflicts and improve reproducibility.
    Therefore, installing conda is also part of setting up your computer to be able to contribute to Python-based modules.


## Install conda

We recommend installing [Miniconda](https://docs.anaconda.com/free/miniconda/index.html) to obtain conda.
Miniconda is lightweight version of the full conda platform and includes the conda tool itself, Python, and a few other commonly-used packages.

To install Miniconda, [download the installer for your operating system](https://docs.anaconda.com/free/miniconda/), and follow all instructions.

  - If you are on a macOS computer, be sure to download one of the links ending in `pkg`, _not `bash`_:
    - If you are on an Apple Silicon (M1-3) mac, download `Miniconda3 macOS Apple M1 64-bit pkg`
    - If you are on an Intel mac, download `Miniconda3 macOS Intel x86 64-bit pkg`

If you already have conda on your system, you do not need to re-install it.

## Set up conda

Next, you will need to set certain conda settings and install a few packages that will allow you to contribute to OpenScPCA in general.

1. [Open a terminal (command line prompt)](../../software-platforms/general-tools/using-the-terminal.md) to interact with conda.

1. Copy and paste the following code into the terminal, and hit enter.
These commands will set the [recommended channels](https://docs.conda.io/projects/conda/en/latest/user-guide/concepts/channels.html) conda should use to search for software.
    - If you get the error `conda: command not found`, you may need to try again in a new terminal window.
    If this doesn't help, you can [get help here](STUB_LINK somewhere to get help).

    ```sh
    conda config --add channels defaults
    conda config --add channels bioconda
    conda config --add channels conda-forge
    conda config --set channel_priority strict
    ```


1. The last step is to install the packages that you will need to contribute to OpenScPCA into your `base` conda environment.
`base` is the default conda environment, the one that will be active when you first open your terminal.
Copy and paste the following command into the terminal, and hit enter.

    ```sh
    conda install awscli conda-lock jq pre-commit
    ```
    <!-- Do we want to suggest this instead? `conda env update --name base --file environment.yml`? -->

    - The [`awscli` package](https://docs.aws.amazon.com/cli/latest/userguide/cli-chap-welcome.html/) will allow you to interact with [data stored in the Amazon Web Services (AWS) S3 bucket](STUB_LINK)
    - The [`conda-lock` package](https://conda.github.io/conda-lock/) is used to create fully reproducible, cross-platform [conda environments for analysis modules](../../contributing-to-analyses/determining-requirements/determining-software-requirements.md#managing-software-dependencies-in-Python-with-conda)
    - The [`jq` package](https://jqlang.github.io/jq/) provides JSON parsing capabilities
    - The [`pre-commit` package](https://pre-commit.com) will allow you to use [pre-commit hooks when contributing to analysis modules](STUB_LINK)

    !!! note
        You may be prompted to enter **`y`** or **`n`** (yes or no) during this setup.
        If/when this prompt appears, you should hit **`y`** to give conda permissions to proceed.


1. At the end of the installation, you should see these messages in the terminal which indicate successful installation:

    <figure markdown="span">
        ![Conda completed install output.](../../img/conda-success.png){width="275"}
    </figure>


All set! 🎉
