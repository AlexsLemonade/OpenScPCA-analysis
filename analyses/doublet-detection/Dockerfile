# Dockerfile for the doublet-detection module
FROM bioconductor/r-ver:3.19
# Labels following the Open Containers Initiative (OCI) recommendations
# For more information, see https://specs.opencontainers.org/image-spec/annotations/?v=v1.0.1
LABEL org.opencontainers.image.authors="OpenScPCA scpca@ccdatalab.org"
LABEL org.opencontainers.image.source="https://github.com/AlexsLemonade/OpenScPCA-analysis/tree/main/analyses/doublet-detection"

# set a name for the conda environment
ARG ENV_NAME=openscpca-doublet-detection

# Set an environment variable to allow checking if we are in an OpenScPCA container
ENV OPENSCPCA_DOCKER=TRUE

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
  && conda clean --force-pkgs-dirs --all --yes \
  && echo ". /opt/conda/etc/profile.d/conda.sh && conda activate base" >> /etc/skel/.bashrc \
  && echo ". /opt/conda/etc/profile.d/conda.sh && conda activate base" >> ~/.bashrc

# Install conda-lock
RUN conda install --channel=conda-forge --name=base conda-lock

# Install pandoc dependency
RUN apt-get -y update &&  \
  DEBIAN_FRONTEND=noninteractive \
  apt-get install --no-install-recommends -y \
  pandoc \
  && rm -rf /var/lib/apt/lists/*

# Install renv
RUN Rscript -e "install.packages('renv')"

# Disable the renv cache to install packages directly into the R library
ENV RENV_CONFIG_CACHE_ENABLED=FALSE

# Copy conda lock file to image
COPY conda-lock.yml conda-lock.yml

# restore from conda-lock.yml file and clean up to reduce image size
RUN conda-lock install -n ${ENV_NAME} conda-lock.yml && \
  conda clean --all --yes

# Copy the renv.lock file from the host environment to the image
COPY renv.lock renv.lock

# Install Rhtslib first to prevent package installation errors
RUN Rscript -e 'BiocManager::install("Rhtslib")'

# restore from renv.lock file and clean up to reduce image size
RUN Rscript -e 'renv::restore()' \
  && rm -rf ~/.cache/R/renv \
  && rm -rf /tmp/downloaded_packages \
  && rm -rf /tmp/Rtmp*

# Complete installation of zellkonverter conda env
ENV BASILISK_EXTERNAL_DIR /usr/local/renv/basilisk
RUN Rscript -e "proc <- basilisk::basiliskStart(env = zellkonverter::zellkonverterAnnDataEnv(), testload = 'anndata'); \
  basilisk::basiliskStop(proc); \
  basilisk.utils::cleanConda()"

# Activate conda environment on bash launch
RUN echo "conda activate ${ENV_NAME}" >> ~/.bashrc

# Set CMD to bash to activate the environment when launching
CMD ["/bin/bash"]
