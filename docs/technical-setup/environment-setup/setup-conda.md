# Install and set up conda

## What is conda?

The OpenScPCA project uses [conda](https://docs.anaconda.com/free/miniconda/index.html) to set up your software environment.
conda is a command-line software management tool which helps you install and track specific versions of software.
It also allows you to have multiple software environments with different sets of packages on the same computer.

This page provides instructions on how to install conda and use it to install certain software you will need to contribute to OpenScPCA.


### Why use conda?

There are two main reasons we use conda for OpenScPCA:

- conda provides a "one stop shop" for installing lots of different software
    - Rather than having to figure out how to install every new program and its dependencies on its own, conda can handle it all for you.
    - You'll use conda to install the software dependencies you'll need to contribute to OpenScPCA.
- conda allows you set up different software environments for different projects
    - For example, you may have two projects that require different versions of the same package.
    With conda, you can create separate, fully isolated software environments for each project with different package versions.
    - Python-based OpenScPCA analysis modules will use different conda environments to prevent conflicts and improve reproducibility.


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

1. Open a terminal (command line prompt) to interact with conda.
The application to open depends on your operating system:
    - _If you are on a macOS machine_, open the Terminal application.
    Please refer to [this documentation](STUB_LINK) about how to open and use the Terminal.

    - _If you are on a Windows machine_, open the conda prompt.
    To launch the prompt, SOMETHING SOMETHING START MENU.
    TODO: Add screenshot.

1. Copy and paste the following code into the terminal, and hit enter.
These commands will set the [recommended channels](https://docs.conda.io/projects/conda/en/latest/user-guide/concepts/channels.html) conda should use to search for software.
    - If you are on a macOS machine and get the error `conda: command not found`, you may need to try again in a new terminal window.
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
    conda install awscli jq pre-commit
    ```

    - The [`awscli` package](https://pypi.org/project/awscli/) will allow you to interact with [data stored in the Amazon Web Services (AWS) S3 bucket](STUB_LINK)
    - The [`jq` package](https://pypi.org/project/jq/) provides JSON parsing capabilities
    - The [`pre-commit`](https://pypi.org/project/pre-commit/) package will allow you to use [pre-commit hooks when contributing to analysis modules](STUB_LINK)

    !!! note
        You may be prompted to enter **`y`** or **`n`** (yes or no) during this setup.
        If/when this prompt appears, you should hit **`y`** to give conda permissions to proceed.


1. At the end of the installation, you should see these messages in the terminal which indicate successful installation:

    <figure markdown="span">
        ![Conda completed install output.](../../img/conda-success.png){width="275"}
    </figure>


All set! 🎉
