# Dockerfile for metacell analysis
FROM bioconductor/r-ver:3.19

# Labels following the Open Containers Initiative (OCI) recommendations
# For more information, see https://specs.opencontainers.org/image-spec/annotations/?v=v1.0.1
LABEL org.opencontainers.image.authors="OpenScPCA scpca@ccdatalab.org"
LABEL org.opencontainers.image.source="https://github.com/AlexsLemonade/OpenScPCA-analysis/tree/main/analyses/metacells"

# Set an environment variable to allow checking if we are in an OpenScPCA container
ENV OPENSCPCA_DOCKER=TRUE

# set a name for the conda environment
ARG ENV_NAME=openscpca-metacells

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
RUN Rscript -e 'install.packages("renv")' \
  && rm -rf /tmp/downloaded_packages \
  && rm -rf /tmp/Rtmp*
ENV RENV_CONFIG_CACHE_ENABLED=FALSE

# restore from conda-lock.yml file and clean up to reduce image size
COPY conda-lock.yml conda-lock.yml
RUN conda-lock install -n ${ENV_NAME} conda-lock.yml \
  && conda clean --all --yes

# restore from renv.lock file and clean up to reduce image size
COPY renv.lock renv.lock
RUN Rscript -e 'renv::restore()' \
  && rm -rf ~/.cache/R/renv \
  && rm -rf /tmp/downloaded_packages \
  && rm -rf /tmp/Rtmp*

# Activate conda environment on bash launch
RUN echo "conda activate ${ENV_NAME}" >> ~/.bashrc

# Set CMD to bash to activate the environment when launching
CMD ["/bin/bash"]
