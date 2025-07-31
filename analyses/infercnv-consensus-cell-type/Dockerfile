# Dockerfile for infercnv-consensus-cell-type module
# note: use 3.19 for compatibility with Rsamtools installation
FROM bioconductor/r-ver:3.19

# Labels following the Open Containers Initiative (OCI) recommendations
# For more information, see https://specs.opencontainers.org/image-spec/annotations/?v=v1.0.1
LABEL org.opencontainers.image.authors="OpenScPCA scpca@ccdatalab.org"
LABEL org.opencontainers.image.source="https://github.com/AlexsLemonade/OpenScPCA-analysis/tree/main/analyses/infercnv-consensus-cell-type"

# Set an environment variable to allow checking if we are in an OpenScPCA container
ENV OPENSCPCA_DOCKER=TRUE

# Disable the renv cache to install packages directly into the R library
ENV RENV_CONFIG_CACHE_ENABLED=FALSE

# Install dependencies for renv libraries
RUN apt-get -y update &&  \
  DEBIAN_FRONTEND=noninteractive \
  apt-get install --no-install-recommends -y \
  pandoc \
  libcurl4-openssl-dev \
  lzma-dev \
  && rm -rf /var/lib/apt/lists/*

# Install renv to enable later package installation
RUN Rscript -e "install.packages('renv')"

# Install Rsamtools directly to avoid build errors when installing from renv.lock
RUN Rscript -e 'BiocManager::install("Rsamtools")'

# Copy the renv.lock file from the host environment to the image
COPY renv.lock renv.lock

# restore from renv.lock file and clean up to reduce image size
RUN Rscript -e 'renv::restore()' && \
  rm -rf ~/.cache/R/renv && \
  rm -rf /tmp/downloaded_packages && \
  rm -rf /tmp/Rtmp*

# Set CMD to bash to activate the environment for any commands
CMD ["/bin/bash"]
