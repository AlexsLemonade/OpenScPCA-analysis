# Creating Docker images

This page provides some guidance on creating Docker images for an OpenScPCA analysis module.

## Prerequisites

Before making a Docker image, you will likely want to have fully [defined the software environment](../../contributing-to-analyses/determining-requirements/determining-software-requirements.md) that you will be using for your analysis module.
Most likely, this will take the form of either:

- an [`renv.lock` file](../../contributing-to-analyses/determining-requirements/determining-software-requirements.md#determining-and-managing-software-dependencies-in-r) if you are working with R, or
- a [conda `environment.yml`](../../contributing-to-analyses/determining-requirements/determining-software-requirements.md#module-specific-conda-environments) and/or [`conda-lock.yml` file](../../contributing-to-analyses/determining-requirements/determining-software-requirements.md#freezing-dependencies-with-conda-lock) if you are working with Python or other software available through [Bioconda](https://bioconda.github.io).

The goal of creating a Docker image is to incorporate these software requirements into a container that can be used to run all parts of the analysis module without any installation steps other than installing Docker and pulling the image.

## Docker template files

When you create an analysis module, a template `Dockerfile` and a `.dockerignore` file will be created in the root of the analysis module.

The  `Dockerfile` includes a `FROM` statement that specifies the base image to use and a couple of `LABEL` statements to populate metadata about the image.

The `.dockerignore` file is used to specify files and directories that should *not* be included in the Docker build environment.
By default, this will be set to ignore all files with the exception of environment files (e.g., `renv.lock`, `environment.yml`, `conda-lock.yml`).
This will keep the build environment small and prevent unnecessary files from being included in the image.

## Base images

The base image you choose will depend on the software requirements of your analysis module.
You should use a publically available base image; ideally one that is in active development and well-maintained in order to ensure that your image is secure and up-to-date.

Below we present a few common base images that you might consider using for your analysis module.

### R images

If you are working with R, you will likely want to use one of the `rocker/r-ver` or `bioconductor/r-ver` images, which are a set of Docker images that provide complete R environments, making it easy to

