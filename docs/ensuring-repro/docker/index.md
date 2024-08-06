# Working with Docker

## What is Docker?

Docker images are a helpful way to ensure that all contributors are able to run analysis in the same environment.
A Docker image compiles all packages, dependencies, and specified versions in a single portable image.
In short, Docker is a key tool we use to ensure _reproducibility_.

By using the same Docker image, contributors can run the same analysis on different machines and get the same results.

!!! note "More information on Docker"

    For more information on Docker, see:

    - [Why you should use Docker in your research](https://blog.zooniverse.org/2018/07/17/why-you-should-use-docker-in-your-research/)
    - [Digging into Data Science Tools: Docker](https://towardsdatascience.com/digging-into-data-science-tools-docker-bbb9f3579c87)
    - [A short guide to using Docker for your data science environment](https://towardsdatascience.com/a-short-guide-to-using-docker-for-your-data-science-environment-912617b3603e)
    - [Docker for Data Science](https://www.datacamp.com/tutorial/docker-for-data-science-introduction)
    - [Docker Desktop documentation](https://docs.docker.com/desktop/use-desktop/)

## Why do I need Docker?

Using Docker with OpenScPCA is optional, but highly recommended.
A docker image will be created and available for each analysis module.
This ensures reproducibility of that analysis module across any machine that is used.

For example, if you work on your analysis locally and on [Lightsail for Research](../../aws/index.md#lightsail-for-research-virtual-computing-with-aws), working in a Docker container ensures reproducible results.

We will also use Docker images to run the analysis modules in the `OpenScPCA-nf` workflow when generating final analysis results.

## How to install Docker

### macOS users

Follow the instructions to [Install Docker Desktop on Mac](https://docs.docker.com/desktop/install/mac-install/) on the Docker website.

#### Apple silicon (M-series) Mac users

After installing Docker, we recommend that Apple silicon Mac users follow these additional steps:

1. Install Rosetta 2 [as described in the "System Requirements" section in the "Mac with Apple silicon" tab](https://docs.docker.com/desktop/install/mac-install/#system-requirements) by running the following in [the terminal](../../getting-started/project-tools/using-the-terminal.md):

    ```sh
    # Install Rosetta 2, following instructions when prompted
    softwareupdate --install-rosetta
    ```

1. Enable Rosetta emulation in the Docker Desktop Settings.
In the ["General" settings](https://docs.docker.com/desktop/settings/mac/#general) tab, make sure the following are turned **on**:
    - "Use Virtualization framework"
    - "Use Rosetta for `x86_64/amd64` emulation on Apple Silicon"

### Windows users

To enable Docker on the WSL 2 side of your computer, you will need to turn on Docker's built-in [WSL 2 feature](https://docs.docker.com/desktop/wsl/).
Once enabled, you will be able to run Docker from the Ubuntu terminal.

Follow the instructions to [install Docker Desktop and enable the WSL 2 feature](https://docs.docker.com/desktop/wsl/#turn-on-docker-desktop-wsl-2) on the Docker website.

