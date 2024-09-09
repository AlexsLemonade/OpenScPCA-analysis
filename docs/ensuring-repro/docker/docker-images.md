# Creating Docker images

This page provides some guidance on creating Docker images for an OpenScPCA analysis module.

Having a Docker image for your analysis module is a key part of ensuring that your analysis is as reproducible as possible.

Note that creating docker images is not required for contributing to OpenScPCA, but it can be useful for testing and sharing your analysis module with others.
Each analysis module will ultimately be updated to include a Docker image before producing final results, but this does not need to be done by the module author.
Data Lab staff will create a Docker image as part of the process of adding the module to the `OpenScPCA-nf` workflow if it has not already been done.

## Prerequisites

Before making a Docker image, you will want to have [defined the software environment](../managing-software/index.md) that you will be using for your analysis module.
Having the requirements prepared will make creating a fully reproducible Docker image much more straightforward.

For the easiest setup, you should have one of the two following kind of environment definitions:

- an [`renv.lock` file](../managing-software/using-renv.md) if you are working with R, or
- a [conda `environment.yml`](../managing-software/using-conda.md) and/or [`conda-lock.yml` file](../managing-software/using-conda.md#conda-and-conda-lock) if you are working with Python or other software available through [Bioconda](https://bioconda.github.io).

You will also need to have [Docker installed on your computer](index.md#how-to-install-docker) to build and test the images.

!!! tip
    If you are new to Docker, you may want to check out the [Docker documentation](https://docs.docker.com/get-started/overview/) for a more in-depth introduction to Docker concepts and commands.

## Analysis module Dockerfiles

### Docker template files

When you [create an analysis module](../../contributing-to-analyses/analysis-modules/creating-a-module.md), a template `Dockerfile` and a `.dockerignore` file will be created in the root of the analysis module.

The `Dockerfile` includes a `FROM` statement that specifies the base image to use and a couple of `LABEL` statements to populate metadata about the image.

The `.dockerignore` file specifies files and directories that should *not* be included in the Docker build environment.
By default, we've set up this file to ignore all files except environment files (e.g., `renv.lock`, `environment.yml`, `conda-lock.yml`).
This will keep the build environment small and prevent unnecessary files from being included in the image.

### Base images

The base image you choose will depend on the software requirements of your analysis module.
You should aim to use a publicly available base image that is in active development and well-maintained in order to ensure that your image is secure and up-to-date.

!!! tip
    You should always use a *versioned* base image to ensure that your analysis module is reproducible, i.e., not the one tagged as `latest`.

We present some recommended base images for R and conda-based environments in the examples below, but you may need to choose a different base image depending on your specific requirements.

### Additional software

If your analysis module requires additional software that is not available in the base image, you will need to install it in the Dockerfile.
The easiest way to do this is likely to be to use the same tool (`renv`, `conda`, etc.) that you used to define the software environment for your analysis module.
If necessary, you can also include statements to install other software using the system package manager (e.g., `apt-get` for Debian-based images).

## Example Dockerfiles

Below are some example Dockerfiles for R and conda-based environments that use the `renv.lock` and `conda-lock.yml` files, respectively, to install required software in the image.
Note that these are relatively minimal examples, and do not include elements like `LABEL` statements that would normally be found in a final OpenScPCA Dockerfile.

### R-based images

If you are working with R, you will likely want to use one of the [`rocker/r-ver` images](https://hub.docker.com/r/rocker/r-ver/tags) or a [`bioconductor/r-ver` image](https://hub.docker.com/r/bioconductor/r-ver/tags).
These are Docker images that provide complete R environments tied to a specific version of R and/or Bioconductor.
The Bioconductor images are especially handy as they include system software dependencies for Bioconductor packages that may not be installed in the `rocker` images.
They also support binary installation of Bioconductor packages, which can dramatically speed up build times.

You will want to include the following steps in your Dockerfile:

- Install `renv` (and `remotes` if you have any GitHub-based packages in your environment).
- Copy the `renv.lock` file from the host environment to the image.
- Use `renv::restore()` to reproduce the R environment from the `renv.lock` file.


A simple `Dockerfile` for an R-based analysis module that uses `renv` for its environment might look like this:

```Dockerfile
# Base image on the Bioconductor 3.19 image
FROM bioconductor/r-ver:3.19

# Install renv to enable later package installation
RUN Rscript -e "install.packages('renv')"

# Disable the renv cache to install packages directly into the R library
ENV RENV_CONFIG_CACHE_ENABLED=FALSE

# Copy the renv.lock file from the host environment to the image
COPY renv.lock renv.lock

# restore from renv.lock file and clean up to reduce image size
RUN Rscript -e 'renv::restore()' \
  && rm -rf ~/.cache/R/renv \
  && rm -rf /tmp/downloaded_packages \
  && rm -rf /tmp/Rtmp*
```


### Conda-based images

If you are working with conda, you will likely want to use one of the [`continuumio/miniconda3` images](https://hub.docker.com/r/continuumio/miniconda3/tags).
Since the conda environment files should define the version of Python itself, you should not need to worry too much about the specific version of the base image, but you should still use a versioned image for reproducibility.

In the Dockerfile, you will want include the following steps:

- Install `conda-lock`.
- Copy the `conda-lock.yml` file from the host environment to the image.
- Use `conda-lock install` to create an environment from the `conda-lock.yml` file.
- Make sure the environment is activated when the container is launched by modifying the `.bashrc` file and setting the `CMD` to `/bin/bash`.

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
RUN conda-lock install -n ${ENV_NAME} \
  && conda clean --all --yes

# Activate conda environment on bash launch
RUN echo "conda activate ${ENV_NAME}" >> ~/.bashrc

# Set CMD to bash to activate the environment when launching
CMD ["/bin/bash"]
```


### Images with both R and conda environments

If your analysis module requires both R and conda environments, you may have to do just a bit more work to set up the Dockerfile, as there is not a single base image that includes both R and conda environments.
<!-- Should we make one, maybe? -->
We recommend starting with an R-based image and then installing conda manually.
In the example below, we have adapted the installation steps used in the [official Miniforge Dockerfile](https://github.com/conda-forge/miniforge-images/blob/master/ubuntu/Dockerfile).

```Dockerfile
# Dockerfile for an analysis module with both R and conda environments
# Base image on the Bioconductor 3.19 image
FROM bioconductor/r-ver:3.19

# set a name for the conda environment
ARG ENV_NAME=openscpca-analysis

# set environment variables to install conda
ENV PATH="/opt/conda/bin:${PATH}"

# Install conda via miniforge
# adapted from https://github.com/conda-forge/miniforge-images/blob/master/ubuntu/Dockerfile
RUN curl -L "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh" -o /tmp/miniforge.sh \
  && bash /tmp/miniforge.sh -b -p /opt/conda \
  && rm -f /tmp/miniforge.sh \
  && conda clean --tarballs --index-cache --packages --yes \
  && find /opt/conda -follow -type f -name '*.a' -delete \
  && find /opt/conda -follow -type f -name '*.pyc' -delete \
  && conda clean --force-pkgs-dirs --all --yes

# Activate conda environments in bash
RUN ln -s /opt/conda/etc/profile.d/conda.sh /etc/profile.d/conda.sh \
  && echo ". /opt/conda/etc/profile.d/conda.sh" >> /etc/skel/.bashrc \
  && echo ". /opt/conda/etc/profile.d/conda.sh" >> ~/.bashrc

# Install conda-lock
RUN conda install --channel=conda-forge --name=base conda-lock \
  && conda clean --all --yes

# Install renv
RUN Rscript -e "install.packages('renv')"

# Disable the renv cache to install packages directly into the R library
ENV RENV_CONFIG_CACHE_ENABLED=FALSE

# Copy conda lock file to image
COPY conda-lock.yml conda-lock.yml

# restore from conda-lock.yml file and clean up to reduce image size
RUN conda-lock install -n ${ENV_NAME} conda-lock.yml \
  && conda clean --all --yes

# Copy the renv.lock file from the host environment to the image
COPY renv.lock renv.lock

# restore from renv.lock file and clean up to reduce image size
RUN Rscript -e 'renv::restore()' \
  && rm -rf ~/.cache/R/renv \
  && rm -rf /tmp/downloaded_packages \
  && rm -rf /tmp/Rtmp*

# Activate conda environment on bash launch
RUN echo "conda activate ${ENV_NAME}" >> ~/.bashrc

# Set CMD to bash to activate the environment when launching
CMD ["/bin/bash"]
```

!!! tip
    Note that we install all of Miniconda, `conda-lock` and `renv` before installing any module-specific packages.
    This is because these steps are not likely to change often, so we can take advantage of Docker's caching to avoid re-installing them every time we build the image.

## Building Docker images

### Local testing

If you have a `Dockerfile` set up for your analysis module, you can build the Docker image locally to test it.
To match the infrastructure that the OpenScPCA project uses you should use the `docker buildx` system and the `--platform linux/amd64` flag when building the image.
The `Dockerfile` should be in the root of the analysis module directory.

The command below will build the Docker image using the `Dockerfile` in the module directory and tag it as `openscpca/your-module:latest`.
You should replace `your-module` with the name of your analysis module.

```bash
# move to the module directory
cd /path/to/your-module
# build the docker image
docker buildx build . -t openscpca/your-module:latest --platform linux/amd64
```

### Automated testing and deployment

Once a module has a complete Docker image, we can add it to the [OpenScPCA automated tests](../workflows/build-docker-gha.md) and the [OpenScPCA Docker registry](https://gallery.ecr.aws/openscpca).
These steps will usually be completed by Data Lab staff.

To add a module's Docker image to testing, we first activate the GitHub Action that builds the Docker image whenever the module's Dockerfile or environment files are updated.
In the `.github/workflows/` directory, there should be a `docker_{your-module}.yml` file that was created when the module was initialized.
To activate the GitHub Action, we will uncomment the `on:` block at the top of the file and update the file lists to include any files that should trigger a rebuild, such as files that define the software environment or that are otherwise copied into the Docker image.

The module name will also be added to the `modules:` list in the `.github/workflows/docker_all-modules.yml` workflow file.
This file is responsible for building all of the Docker images for the OpenScPCA project during periodic testing and when we release tagged versions of the repository.
