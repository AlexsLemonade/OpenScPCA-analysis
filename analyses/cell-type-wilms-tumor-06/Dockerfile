# pull base image
FROM bioconductor/tidyverse:3.19


# Change the default CRAN mirror
RUN echo "options(repos = c(CRAN = 'https://mran.microsoft.com/snapshot/2022-02-01'), download.file.method = 'libcurl')" >> ${R_HOME}/etc/Rprofile.site


# Set global R options
RUN echo "options(repos = 'https://cloud.r-project.org')" > $(R --no-echo --no-save -e "cat(Sys.getenv('R_HOME'))")/etc/Rprofile.site
ENV RETICULATE_MINICONDA_ENABLED=FALSE


RUN R -e "devtools::install_github('enblacar/SCpubr')"
RUN R --no-echo --no-restore --no-save -e "install.packages('remotes')"
RUN R -e "remotes::install_github('satijalab/seurat', 'seurat5', quiet = TRUE)"  # this also install patchwork (and others)
RUN R -e "remotes::install_github('satijalab/azimuth', quiet = TRUE)"            # this also install SingleCellExperiment, DT (and others)


#RUN chmod -R a+rw ${R_HOME}/site-library # so that everyone can dynamically install more libraries within container
#RUN chmod -R a+rw ${R_HOME}/library

# add custom options for rstudio sessions
# make sure sessions stay alive forever
RUN echo "session-timeout-minutes=0" >> /etc/rstudio/rsession.conf
# make sure we get rstudio server logs in the container
# RUN echo $'[*]\nlog-level=warn\nlogger-type=file\n' > /etc/rstudio/logging.conf

# make sure all R related binaries are in PATH in case we want to call them directly
ENV PATH ${R_HOME}/bin:$PATH
