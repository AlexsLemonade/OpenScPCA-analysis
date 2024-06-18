# Dockerfile for doublet-detection module
FROM bioconductor/r-ver:3.19
# Labels following the Open Containers Initiative (OCI) recommendations
# For more information, see https://specs.opencontainers.org/image-spec/annotations/?v=v1.0.1
LABEL org.opencontainers.image.authors="OpenScPCA scpca@ccdatalab.org"
LABEL org.opencontainers.image.source="https://github.com/AlexsLemonade/OpenScPCA-analysis/tree/main/analyses/doublet-detection"

# set a name for the conda environment
ARG ENV_NAME=openscpca-doublet-detection

# set environment variables to install miniconda
ENV PATH="/opt/conda/bin:${PATH}"

# Install Miniconda
# adapted from https://github.com/ContinuumIO/docker-images/blob/main/miniconda3/debian/Dockerfile
RUN curl -O https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh \
  && bash Miniconda3-latest-Linux-x86_64.sh -b -p /opt/conda \
  && ln -s /opt/conda/etc/profile.d/conda.sh /etc/profile.d/conda.sh  \
  && echo ". /opt/conda/etc/profile.d/conda.sh" >> ~/.bashrc  \
  && rm -f Miniconda3-latest-Linux-x86_64.sh \
  && find /opt/conda/ -follow -type f -name '*.a' -delete \
  && find /opt/conda/ -follow -type f -name '*.js.map' -delete

# Install conda-lock
RUN conda install --channel=conda-forge --name=base conda-lock

# Install renv
RUN Rscript -e "install.packages('renv')"

# Disable the renv cache to install packages directly into the R library
ENV RENV_CONFIG_CACHE_ENABLED FALSE

# Copy conda lock file to image
COPY conda-lock.yml conda-lock.yml

# restore from conda-lock.yml file and clean up to reduce image size
RUN conda-lock install -n ${ENV_NAME} conda-lock.yml && \
  conda clean --all --yes

# Copy the renv.lock file from the host environment to the image
COPY renv.lock renv.lock

# restore from renv.lock file and clean up to reduce image size
RUN Rscript -e 'renv::restore()' \
  && rm -rf ~/.cache/R/renv \
  && rm -rf /tmp/downloaded_packages \
  && rm -rf /tmp/Rtmp*

# Activate conda environment on bash launch
RUN echo "conda activate ${ENV_NAME}" >> ~/.bashrc

# Set entrypoint to bash to activate the environment for any commands
ENTRYPOINT ["bash", "-l", "-c"]
