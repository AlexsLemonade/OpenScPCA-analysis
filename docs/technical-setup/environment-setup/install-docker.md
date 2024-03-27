# Install Docker

## What is Docker?

Docker containers are a helpful way to ensure that all contributors are able to run analysis in the same environment.
A docker container compiles all packages, dependencies, and specified versions in a single portable image.

This means that contributors can use the same Docker image to run the same analysis on two different machines and ensure the results are reproducible.

Read more about [Docker images and how we use them in OpenScPCA](STUB-LINK to what are docker images).

## Why do I need Docker?

Using Docker with OpenScPCA is optional, but highly recommended.
A docker container will be created and available for each analysis module.
This ensures reproducibility of that analysis module across any machine that is used.

For example, if you work on your analysis locally and on [Lightsail for Research](STUB_LINK for LSfR), working in a Docker image ensures reproducible results.

## How to install Docker

Docker is available for Mac, Linux, and Windows.
Select your platform and follow the platform-specific instructions provided by Docker [here](https://docs.docker.com/get-docker/).
