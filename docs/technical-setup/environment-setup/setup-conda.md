# Install and setup `conda`

## What is `conda`?

The `OpenScPCA` project uses [`conda`](https://docs.anaconda.com/free/miniconda/index.html) to setup your software environment.
`conda` is a command-line software management tool which helps you install and track specific versions of software in many programming languages.
We also use `conda` to handle [software dependencies for Python-based analysis modules](STUB_LINK).
This page provides instructions on how to install `conda` and use it to install certain software you will need to contribute to `OpenScPCA`.

With `conda`, you can also set up different software environments on your computer.
For example, you may have one project that uses Python version 3.6, and another that uses Python version 3.12.
You can use `conda` to set upip differe

## Install `conda`

We recommend installing [Miniconda](https://docs.anaconda.com/free/miniconda/index.html) to obtain `conda`.
Miniconda is lightweight version of the full `conda` platform and includes the `conda` tool itself, Python, and a few other commonly-used packages.

To install Miniconda, [download the installer for your operating system](https://docs.anaconda.com/free/miniconda/miniconda-install/), and follow all instructions.

## Setup `conda`

To setup `conda` after installation, follow the instructions based on your operating system:

### Setup `conda` on MacOS

1. You can use `conda` via the `Terminal` application.
To launch `Terminal`, search for "terminal" in spotlight and open the application.



1. Copy and paste the following code into `Terminal`, and hit enter.
These commands will set the [recommended channels](https://docs.conda.io/projects/conda/en/latest/user-guide/concepts/channels.html) `conda` should use to search for software.


    !!! info
        If you get an error `conda: command not found`, this means `conda` was not properly installed.
        Please DO SOMETHING TO GET HELP? IS IT OPEN A DISCUSSION? IS IT DM US? IS IT SEE OTHER DOCS FOR WAYS TO GET HELP?

1. In `Terminal`, use `cd` to navigate to the `OpenScPCA-analysis` repository you cloned. To make this easier, you can take the following steps:
    - Type into the `Terminal`
    - Open a `Finder` window
    - Navigate to the repository folder
    - In `Terminal`, type `cd ` (don't forget the space!), but do not press enter yet
    - Click and drag the folder from `Finder` into the `Terminal`
    - Finally, hit enter
    <!-- TODO: Insert movie/gif here!! -->


1. To confirm you are in the correct folder in `Terminal`, type the `ls` command.
You should see this output:
<!-- img of ls output -->

1. install environment packages to create `openscpca` environment

1.



### Setup `conda` on Windows

!!! note
    This option is best for contributors who have previous experience with Git and the command line.
    The OpenScPCA documentation will generally not present how to use Git via the command line.<br><br>
    If you are new to Git, we do not recommend this option.

Follow these steps to install and setup `conda`:

1.



