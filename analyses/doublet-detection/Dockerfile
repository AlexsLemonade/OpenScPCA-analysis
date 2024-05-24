# Base image on the Bioconductor 3.19 image
FROM bioconductor/r-ver:3.19

# Labels following the Open Containers Initiative (OCI) recommendations
# For more information, see https://specs.opencontainers.org/image-spec/annotations/?v=v1.0.1
LABEL org.opencontainers.image.authors="OpenScPCA scpca@ccdatalab.org"
LABEL org.opencontainers.image.source="https://github.com/AlexsLemonade/OpenScPCA-analysis/tree/main/templates/analysis-module"

# Install pandoc and wget
RUN apt-get -y update &&  \
    DEBIAN_FRONTEND=noninteractive \
    apt-get install --no-install-recommends -y \
    pandoc wget\
    && rm -rf /var/lib/apt/lists/*

# Set variables for conda usage, including environment name, and install miniconda3
ENV PATH="/root/miniconda3/bin:${PATH}"
ARG PATH="/root/miniconda3/bin:${PATH}"
ARG ENV_NAME=openscpca-doublet-detection

RUN wget \
    https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh \
    && mkdir /root/.conda \
    && bash Miniconda3-latest-Linux-x86_64.sh \
    && rm -f Miniconda3-latest-Linux-x86_64.sh

# Install renv to enable later package installation
RUN Rscript -e "install.packages('renv')"

# Disable the renv cache to install packages directly into the R library
ENV RENV_CONFIG_CACHE_ENABLED FALSE

# Copy the renv.lock file from the host environment to the image
COPY renv.lock renv.lock

# Restore from renv.lock file and clean up to reduce image size
RUN Rscript -e 'renv::restore()' && \
  rm -rf ~/.cache/R/renv && \
  rm -rf /tmp/downloaded_packages && \
  rm -rf /tmp/Rtmp*

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