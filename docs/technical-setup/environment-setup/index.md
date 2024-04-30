# Environment setup

This section of the documentation walks you through additional setup steps that will allow you to contribute to OpenScPCA.

This section assumes that you have already taken these setup steps:

- Download, install, and set up your [Git client](../install-a-git-client.md)
- [Fork the `AlexsLemonade/OpenScPCA-analysis` repository](../fork-the-repo.md)
- [Clone your fork](../clone-the-repo.md) to your computer

While contributing to OpenScPCA, you may end up working on your local computer and/or on a Lightsail for Research instance.
There are different instructions for getting set up to work on either platform.

## Setting up locally

If you are planning to work on your local computer, you will need to take these additional steps:

- [Download, install, and set up conda](./setup-conda.md)
    - OpenScPCA uses conda to manage your software environment and dependencies.
    - Installing conda also provides you with Python and other tools you need to write or contribute to Python-based modules
- [Set up `pre-commit`](./setup-precommit.md)
    - This will allow you to make contributions in Git
- [Configure the AWS CLI](./configure-aws-cli.md)
    - You can only complete this step if the Data Lab has [created an Amazon Web Services account for you](../../getting-started/accessing-resources/index.md).
- [Optional] [Install Docker](../../software-platforms/docker/index.md#how-to-install-docker)
    - This allows you to run your analysis modules inside a Docker container, ensuring reproducibility across machines.
- [Optional] If you plan on writing or contributing to R-based modules, you should also [download and install R and RStudio, along with a few handy R packages](./install-r-rstudio.md)

## Setting up on Lightsail for Research

If you are planning to work on a Lightsail for Research instance, please refer to the [Lightsail for Research documentation](../../software-platforms/aws/index.md#lightsail-for-research-virtual-computing-with-aws).
