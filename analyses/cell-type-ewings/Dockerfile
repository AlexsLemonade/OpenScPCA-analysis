# Dockerfile for cell-type-ewings module
FROM bioconductor/r-ver:devel-R-4.4.0

# Labels following the Open Containers Initiative (OCI) recommendations
# For more information, see https://specs.opencontainers.org/image-spec/annotations/?v=v1.0.1
LABEL org.opencontainers.image.title="openscpca/cell-type-ewings"
LABEL org.opencontainers.image.description="Docker image for the OpenScPCA analysis module 'cell-type-ewings'"
LABEL org.opencontainers.image.authors="OpenScPCA scpca@ccdatalab.org"
LABEL org.opencontainers.image.source="https://github.com/AlexsLemonade/OpenScPCA-analysis/tree/main/analyses/cell-type-ewings"

# Disable the renv cache to install packages directly into the R library
ENV RENV_CONFIG_CACHE_ENABLED FALSE
# set a name for the conda environment
ARG ENV_NAME=openscpca-cell-type-ewings

# set environment variables to install miniconda
ENV PATH="/root/miniconda3/bin:${PATH}"
ARG PATH="/root/miniconda3/bin:${PATH}"

### Renv setup

# Copy the renv.lock file from the host environment to the image
COPY renv.lock renv.lock

# Install renv to enable later package installation
RUN Rscript -e "install.packages('renv')"

# Install dependencies for renv library
RUN apt-get -y update &&  \
    DEBIAN_FRONTEND=noninteractive \
    apt-get install --no-install-recommends -y \
    pandoc \
    wget \
    && rm -rf /var/lib/apt/lists/*

# restore from renv.lock file and clean up to reduce image size
RUN Rscript -e 'renv::restore()' && \
  rm -rf ~/.cache/R/renv && \
  rm -rf /tmp/downloaded_packages && \
  rm -rf /tmp/Rtmp*

### Conda setup
# Install miniconda
RUN wget \
    https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh \
    && mkdir /root/.conda \
    && bash Miniconda3-latest-Linux-x86_64.sh -b \
    && rm -f Miniconda3-latest-Linux-x86_64.sh

COPY environment.yml environment.yml

# restore from conda-lock.yml file and clean up to reduce image size
RUN conda env create -f environment.yml -n ${ENV_NAME} && \
  conda clean --all --yes

# Activate conda environment on bash launch
RUN conda init
RUN echo "conda activate ${ENV_NAME}" >> ~/.bashrc

# Set entrypoint to bash to activate the environment for any commands
ENTRYPOINT ["/bin/bash"]
