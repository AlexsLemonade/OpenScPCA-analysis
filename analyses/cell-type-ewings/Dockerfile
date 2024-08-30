# Dockerfile for cell-type-ewings module
FROM bioconductor/r-ver:3.19

# Labels following the Open Containers Initiative (OCI) recommendations
# For more information, see https://specs.opencontainers.org/image-spec/annotations/?v=v1.0.1
LABEL org.opencontainers.image.title="openscpca/cell-type-ewings"
LABEL org.opencontainers.image.description="Docker image for the OpenScPCA analysis module 'cell-type-ewings'"
LABEL org.opencontainers.image.authors="OpenScPCA scpca@ccdatalab.org"
LABEL org.opencontainers.image.source="https://github.com/AlexsLemonade/OpenScPCA-analysis/tree/main/analyses/cell-type-ewings"

# Set an environment variable to allow checking if we are in an OpenScPCA container
ENV OPENSCPCA_DOCKER=TRUE

# Disable the renv cache to install packages directly into the R library
ENV RENV_CONFIG_CACHE_ENABLED=FALSE
# set a name for the conda environment
ARG ENV_NAME=openscpca-cell-type-ewings

# set environment variables to install miniconda
ENV PATH="/opt/conda/bin:${PATH}"

# Install dependencies for renv library
RUN apt-get -y update &&  \
  DEBIAN_FRONTEND=noninteractive \
  apt-get install --no-install-recommends -y \
  pandoc \
  wget \
  && rm -rf /var/lib/apt/lists/*

# Install conda via miniforge
# adapted from https://github.com/conda-forge/miniforge-images/blob/master/ubuntu/Dockerfile
RUN curl -L "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh" -o /tmp/miniforge.sh \
  && bash /tmp/miniforge.sh -b -p /opt/conda \
  && rm -f /tmp/miniforge.sh \
  && conda clean --tarballs --index-cache --packages --yes \
  && find /opt/conda -follow -type f -name '*.a' -delete \
  && find /opt/conda -follow -type f -name '*.pyc' -delete \
  && conda clean --force-pkgs-dirs --all --yes \
  && echo ". /opt/conda/etc/profile.d/conda.sh && conda activate base" >> /etc/skel/.bashrc \
  && echo ". /opt/conda/etc/profile.d/conda.sh && conda activate base" >> ~/.bashrc

# add conda-lock
RUN conda install conda-lock && conda clean --all --yes

# Install renv to enable later package installation
RUN Rscript -e "install.packages('renv')"

# Copy conda lock file to image
COPY conda-lock.yml conda-lock.yml

# restore from conda-lock.yml file and clean up to reduce image size
RUN conda-lock install -n ${ENV_NAME} conda-lock.yml \
  && conda clean --all --yes

# Activate conda environment on bash launch
RUN echo "conda activate ${ENV_NAME}" >> ~/.bashrc

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
