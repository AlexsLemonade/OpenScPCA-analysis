# Docker image for the simulate-sce analysis module
FROM bioconductor/bioconductor_docker:3.18

# Labels following the Open Containers Initiative (OCI) recommendations
# For more information, see https://specs.opencontainers.org/image-spec/annotations/?v=v1.0.1
LABEL org.opencontainers.image.title="ccdl/openscpca-simulate-sce"
LABEL org.opencontainers.image.description="Docker image for the OpenScPCA analysis module 'simulate-sce'"
LABEL org.opencontainers.image.authors="OpenScPCA scpca@ccdatalab.org"
LABEL org.opencontainers.image.source="https://github.com/AlexsLemonade/OpenScPCA-analysis/tree/main/analyses/simulate-sce"

ENV RENV_CONFIG_CACHE_ENABLED FALSE
RUN Rscript -e "install.packages(c('remotes', 'renv'))"

WORKDIR /usr/local/renv
COPY renv.lock renv.lock
# install packages with renv and remove caches
RUN Rscript -e "renv::restore()" && \
  rm -rf ~/.local/share/renv && \
  rm -rf /tmp/downloaded_packages && \
  rm -rf /tmp/Rtmp*

# Complete installation of zellkonverter conda env
ENV BASILISK_EXTERNAL_DIR /usr/local/renv/basilisk
RUN Rscript -e "proc <- basilisk::basiliskStart(env = zellkonverter::zellkonverterAnnDataEnv(), testload = 'anndata'); \
  basilisk::basiliskStop(proc)"

# set final workdir for commands
WORKDIR /home
