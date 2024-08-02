# pull base image
FROM bioconductor/tidyverse:3.19

# Set global R options
RUN echo "options(repos = 'https://cloud.r-project.org')" > $(R --no-echo --no-save -e "cat(Sys.getenv('R_HOME'))")/etc/Rprofile.site
ENV RETICULATE_MINICONDA_ENABLED=FALSE

RUN R --no-echo --no-restore --no-save -e "install.packages('remotes')"

RUN R -e "devtools::install_github('enblacar/SCpubr')"
RUN R -e "remotes::install_github('satijalab/seurat', 'seurat5', quiet = TRUE)"  # this also install patchwork (and others)
RUN R -e "remotes::install_github('satijalab/azimuth', quiet = TRUE)"            # this also install SingleCellExperiment, DT (and others)
RUN R  -e "remotes::install_github('cancerbits/DElegate')"
RUN R  -e "install.packages('viridis')"
RUN R  -e "install.packages('ggplotify')"
RUN R -e "BiocManager::install('edgeR')"

# make sure all R related binaries are in PATH in case we want to call them directly
ENV PATH ${R_HOME}/bin:$PATH
