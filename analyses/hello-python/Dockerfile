# Dockerfile for the hello-python analysis
FROM condaforge/miniforge3:24.9.0-0

# Labels following the Open Containers Initiative (OCI) recommendations
# For more information, see https://specs.opencontainers.org/image-spec/annotations/?v=v1.0.1
LABEL org.opencontainers.image.authors="OpenScPCA scpca@ccdatalab.org"
LABEL org.opencontainers.image.source="https://github.com/AlexsLemonade/OpenScPCA-analysis/tree/main/analyses/hello-python"

# Set an environment variable to allow checking if we are in an OpenScPCA container
ENV OPENSCPCA_DOCKER=TRUE

# set a name for the conda environment
ARG ENV_NAME=openscpca-hello-python

# Install conda-lock to enable later package installation
RUN conda install --channel=conda-forge --name=base conda-lock

# Copy the conda-lock.yml file from the host environment to the image
COPY conda-lock.yml conda-lock.yml

# restore from conda-lock.yml file and clean up to reduce image size
RUN conda-lock install -n ${ENV_NAME} && \
  conda clean --all --yes

# Activate conda environment on bash launch
RUN echo "conda activate ${ENV_NAME}" >> ~/.bashrc

# Set CMD to bash to activate the environment when launching
CMD ["/bin/bash"]
