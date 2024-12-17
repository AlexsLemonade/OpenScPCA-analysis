# Dockerfile for cell-type-consensus module
FROM bioconductor/r-ver:3.20

# Labels following the Open Containers Initiative (OCI) recommendations
# For more information, see https://specs.opencontainers.org/image-spec/annotations/?v=v1.0.1
LABEL org.opencontainers.image.title="openscpca/cell-type-consensus"
LABEL org.opencontainers.image.description="Docker image for the OpenScPCA analysis module 'cell-type-consensus'"
LABEL org.opencontainers.image.authors="OpenScPCA scpca@ccdatalab.org"
LABEL org.opencontainers.image.source="https://github.com/AlexsLemonade/OpenScPCA-analysis/tree/main/analyses/cell-type-consensus"

# Set an environment variable to allow checking if we are in an OpenScPCA container
ENV OPENSCPCA_DOCKER=TRUE

# Disable the renv cache to install packages directly into the R library
ENV RENV_CONFIG_CACHE_ENABLED=FALSE

# Install dependencies for renv library
RUN apt-get -y update &&  \
  DEBIAN_FRONTEND=noninteractive \
  apt-get install --no-install-recommends -y \
  pandoc \
  && rm -rf /var/lib/apt/lists/*

# Install aws cli 
RUN curl "https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip" -o "awscliv2.zip"
RUN unzip awscliv2.zip
RUN ./aws/install

# Install renv to enable later package installation
RUN Rscript -e "install.packages('renv')"

# Copy the renv.lock file from the host environment to the image
COPY renv.lock renv.lock

# Temporarily install Rhtslib separately
RUN Rscript -e 'BiocManager::install("Rhtslib")'

# restore from renv.lock file and clean up to reduce image size
RUN Rscript -e 'renv::restore()' && \
  rm -rf ~/.cache/R/renv && \
  rm -rf /tmp/downloaded_packages && \
  rm -rf /tmp/Rtmp*

# Set CMD to bash to activate the environment for any commands
CMD ["/bin/bash"]
