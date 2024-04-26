# Install Docker

## What is Docker?

Docker images are a helpful way to ensure that all contributors are able to run analysis in the same environment.
A Docker image compiles all packages, dependencies, and specified versions in a single portable image.
In short, Docker is a key tool we use to ensure _reproducibility_.

By using the same Docker image, contributors can run the same analysis on different machines and get the same results.

!!! note "More information on Docker"

    For more information on Docker, see:

    - [Docker for Data Science](https://www.datacamp.com/tutorial/docker-for-data-science-introduction)
    - [Why you should use Docker in your research](https://blog.zooniverse.org/2018/07/17/why-you-should-use-docker-in-your-research/)
<!--
    - [Docker images and how we use them in OpenScPCA](STUB_LINK to what are docker images)
-->

## Why do I need Docker?

Using Docker with OpenScPCA is optional, but highly recommended.
A docker image will be created and available for each analysis module.
This ensures reproducibility of that analysis module across any machine that is used.

For example, if you work on your analysis locally and on [Lightsail for Research](../../software-platforms/aws/index.md#lsfr-virtual-computing-with-aws), working in a Docker container ensures reproducible results.

We will also use Docker images to run the analysis modules in the `OpenScPCA-nf` workflow when generating final analysis results.

## How to install Docker

Follow the platform-specific instructions [to download and install Docker Desktop](https://docs.docker.com/get-docker/) from the Docker website.
