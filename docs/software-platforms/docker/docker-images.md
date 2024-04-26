# Creating Docker images

This page provides some guidance on creating Docker images for an OpenScPCA analysis module.

Note that creating docker images is not required for contributing to OpenScPCA, but it can be useful for testing and sharing your analysis module with others.
Each analysis module will ultimately be updated to include a Docker image to allow for running the analysis module as part of the `OpenScPCA-nf` workflow, but this does not need to be done by the module author.
Data Lab staff will create a Docker image as part of the process of adding the module to the `OpenScPCA-nf` workflow if it has not already been done.

## Prerequisites

Before making a Docker image, you will likely want to have fully [defined the software environment](../../contributing-to-analyses/determining-requirements/determining-software-requirements.md) that you will be using for your analysis module.
Most likely, this will take the form of either:

- an [`renv.lock` file](../../contributing-to-analyses/determining-requirements/determining-software-requirements.md#determining-and-managing-software-dependencies-in-r) if you are working with R, or
- a [conda `environment.yml`](../../contributing-to-analyses/determining-requirements/determining-software-requirements.md#module-specific-conda-environments) and/or [`conda-lock.yml` file](../../contributing-to-analyses/determining-requirements/determining-software-requirements.md#freezing-dependencies-with-conda-lock) if you are working with Python or other software available through [Bioconda](https://bioconda.github.io).

The goal of creating a Docker image is to incorporate these software requirements into a container that can be used to run all parts of the analysis module without any installation steps other than installing Docker and pulling the image.

You will also need to have [Docker installed on your computer](../../technical-setup/environment-setup/install-docker.md) to build and test the images.

## Docker template files

When you create an analysis module, a template `Dockerfile` and a `.dockerignore` file will be created in the root of the analysis module.

The  `Dockerfile` includes a `FROM` statement that specifies the base image to use and a couple of `LABEL` statements to populate metadata about the image.

The `.dockerignore` file is used to specify files and directories that should *not* be included in the Docker build environment.
By default, this will be set to ignore all files with the exception of environment files (e.g., `renv.lock`, `environment.yml`, `conda-lock.yml`).
This will keep the build environment small and prevent unnecessary files from being included in the image.

## Base images

The base image you choose will depend on the software requirements of your analysis module.
You should aim to use a publicly available base image that is in active development and well-maintained in order to ensure that your image is secure and up-to-date.
You should always use a *versioned* base image to ensure that your analysis module is reproducible, i.e., not the one tagged as `latest`.

## Additional software

If your analysis module requires additional software that is not available in the base image, you will need to install it in the Dockerfile.
The easiest way to do this is likely to be to use the same tool (`renv`, `conda`, etc.) that you used to define the software environment for your analysis module.
If necessary, you can also include statements to install other software using the system package manager (e.g., `apt-get` for Debian-based images).

## R-based images

If you are working with R, you will likely want to use one of the [`rocker/r-ver` images](https://hub.docker.com/r/rocker/r-ver/tags) or a[`bioconductor/r-ver` image](https://hub.docker.com/r/bioconductor/r-ver/tags) images, which are a set of Docker images that provide complete R environments, tied to a specific version of R and/or Bioconductor.
The Bioconductor images are especially handy as they include system software dependencies for Bioconductor packages that may not be installed in the `rocker` images.
They also support binary installation of Bioconductor packages, which can dramatically speed up build times.

You will want to include the following steps in your Dockerfile:

- Install `renv` (and `remotes` if you have any GitHub-based packages in your environment).
- Copy the `renv.lock` file from the host environment to the image.
- Use `renv::restore()` to reproduce the R environment from the `renv.lock` file.


A simple `Dockerfile` for an R-based analysis module that uses `renv` for its environment might look like this:

```Dockerfile
# Base image on the Bioconductor 3.18 image
FROM bioconductor/r-ver:3.18

# Install renv to enable later package installation
RUN Rscript -e "install.packages('renv')"

# Disable the renv cache to install packages directly into the R library
ENV RENV_CONFIG_CACHE_ENABLED FALSE

# Copy the renv.lock file from the host environment to the image
COPY renv.lock renv.lock

# restore from renv.lock file and clean up to reduce image size
RUN Rscript -e 'renv::restore()' && \
  rm -rf ~/.cache/R/renv && \
  rm -rf /tmp/downloaded_packages && \
  rm -rf /tmp/Rtmp*
```


## Conda-based images

If you are working with conda, you will likely want to use one of the [`continuumio/miniconda3` images](https://hub.docker.com/r/continuumio/miniconda3/tags).
Since the conda environment files should define the version of Python itself, you should not need to worry too much about the specific version of the base image, but you should still use a versioned image for reproducibility.

In the Dockerfile, you will want include the following steps:

- Install `conda-lock`.
- Copy the `conda-lock.yml` file from the host environment to the image.
- Use `conda-lock install` to create an environment from the `conda-lock.yml` file.
- Make sure the environment is activated when the container is launched by modifying the `.bashrc` file and setting the `ENTRYPOINT` to `bash -l -c`.

A simple `Dockerfile` for a conda-based analysis module that uses `conda-lock` for its environment might look like this:

```Dockerfile
# Dockerfile for the hello-python analysis
FROM continuumio/miniconda3:24.1.2-0

# set a name for the conda environment
ARG ENV_NAME=openscpca-analysis

# Install conda-lock to enable later package installation
RUN conda install --channel=conda-forge --name=base conda-lock

# Copy the conda-lock.yml file from the host environment to the image
COPY conda-lock.yml conda-lock.yml

# restore from conda-lock.yml file and clean up to reduce image size
RUN conda-lock install -n ${ENV_NAME} && \
  conda clean --all --yes

# Activate conda environment on bash launch
RUN echo "conda activate ${ENV_NAME}" >> ~/.bashrc

# Set entrypoint to bash to activate the environment for any commands
ENTRYPOINT ["bash", "-l", "-c"]
```
