# Technical setup

This section provides instructions for technical aspects of setup.

All of these steps are required setup when working locally on your computer (or server).

If you are planning to work on a Lightsail for Research instance, please instead refer to the [Lightsail for Research documentation](../software-platforms/aws/index.md#lightsail-for-research-virtual-computing-with-aws) for setup instructions.


## Local setup steps

To prepare to contribute to OpenScPCA, please do the following in order:

1. If you are on a Windows machine, you will need to [install WSL 2](./install-wsl-2.md), as described below.
1. Download and set up a [Git client](./install-a-git-client.md)
1. [Fork the `AlexsLemonade/OpenScPCA-analysis` repository](./fork-the-repo.md)
1. [Clone your fork](./clone-the-repo.md) to your computer
1. Set up [additional dependencies](./environment-setup/index.md) that you'll need to contribute to OpenScPCA:
      - [Download, install, and set up conda](./environment-setup/setup-conda.md)
        - OpenScPCA uses conda to manage your software environment and dependencies.
        - Installing conda also provides you with Python and other tools you need to write or contribute to Python-based modules
      - [Set up `pre-commit`](./environment-setup/setup-precommit.md)
        - This will allow you to make contributions in Git
      - [Configure the AWS command line interface (CLI)](./environment-setup/configure-aws-cli.md)
        - You can only complete this step if the Data Lab has [created an Amazon Web Services account for you](../getting-started/accessing-resources/index.md)
      - [Optional] If you plan on writing or contributing to R-based modules, you should also [download and install R and RStudio, along with a few handy R packages](./environment-setup/install-r-rstudio.md)
      - [Optional] [Install Docker](../software-platforms/docker/index.md#how-to-install-docker)
        - This allows you to run your analysis modules inside a Docker container, ensuring reproducibility across machines

## Special considerations for Windows users

Due to the OpenScPCA project's technical environment needs, Windows users will [need to install and use Windows Subsystem for Linux 2 (WSL 2)](./install-wsl-2.md).

[WSL 2](https://learn.microsoft.com/en-us/windows/wsl/about) lets you run a separate Linux (Ubuntu) environment on your Windows computer.
This includes a separate Linux file system that is mostly distinct from your regular Windows file system.
When contributing to OpenScPCA, all data analysis will happen on the WSL 2 side of your computer, which means:

- Your repository needs to be cloned into your WSL 2 file system
- All applications (e.g., GitKraken and RStudio) need to be installed as Linux applications on your WSL 2 system

_The OpenScPCA project does not support the use of Windows machines without WSL 2._
Note that all contributors have access to [Lightsail for Research (LSfR)](../software-platforms/aws/index.md#lsfr-virtual-computing-with-aws) instances which run Linux (specifically, Ubuntu), so you can use LSfR if you prefer not to install WSL 2.

!!! note "Learn more about WSL 2"
    For more information about WSL 2, please see this [Introduction to Windows Subsystem for Linux](https://learn.microsoft.com/en-us/training/modules/wsl-introduction/).

