# Environment setup

This section of the documentation walks you through additional setup steps that will allow you to contribute to OpenScPCA.

This section assumes that you have already taken these setup steps:

- Downloaded your [Git client of interest](../install-a-git-client.md)
    - If you are working on a Windows machine, you first need to install the [Windows Subsystem for Linux](STUB_LINK WSL instructions) before setting up Git.
    - TODO: Is this still true? If you are working on a [Lightsail for Research instance](STUB_LINK for LSfR), note that GitKraken will be pre-installed for you, but you will still need to sign in with your GitHub account.
- [Forked the `AlexsLemonade/OpenScPCA-analysis` repository](../fork-the-repo.md)
- [Cloned your fork](../clone-the-repo.md) to your computer

While contributing to OpenScPCA, you may end up working on your local computer and/or on a Lightsail for Research instance.
There are different instructions for getting set up to work on either platform.

## Setting up locally

If you are planning to work on your local computer, you will need to take these additional steps:

- [Download, install, and set up conda](./setup-conda.md)
    - OpenScPCA uses conda to manage your software environment and dependencies.
    - Installing conda also provides you with Python and other tools you need to write or contribute to Python-based modules
- [Set up `pre-commit`](./setup-precommit.md)
    - This will allow you to make contributions in Git
- [Optional] If you plan on writing or contributing to R-based modules, you should also [download and install R and RStudio, along with a few handy R packages](./install-r-rstudio.md)

## Setting up on Lightsail for Research

If you are planning to work on a Lightsail for Research instance, please refer to the [Lightsail for Research documentation](STUB_LINK for LSfR).
