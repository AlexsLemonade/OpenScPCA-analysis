# pull base image
FROM rocker/tidyverse:4.3.0


RUN apt-get update \
  && apt-get install -y --no-install-recommends \
    lbzip2 \
    libfftw3-dev \
    libgdal-dev \
    libgeos-dev \
    libgsl0-dev \
    libgl1-mesa-dev \
    libglu1-mesa-dev \
    libhdf4-alt-dev \
    libhdf5-dev \
    libjq-dev \
   libgeos-dev \
    libpq-dev \
    libproj-dev \
    libprotobuf-dev \
    libnetcdf-dev \
    libsqlite3-dev \
    libssl-dev \
    libudunits2-dev \
    netcdf-bin \
    postgis \
    protobuf-compiler \
    sqlite3 \
    tk-dev \
    unixodbc-dev

# Change the default CRAN mirror
RUN echo "options(repos = c(CRAN = 'https://mran.microsoft.com/snapshot/2022-02-01'), download.file.method = 'libcurl')" >> ${R_HOME}/etc/Rprofile.site

RUN apt-get update -qq \
  && apt-get -y --no-install-recommends install \
  htop \
  nano \
  libigraph-dev \
  libcairo2-dev \
  libxt-dev \
  libcurl4-openssl-dev \
  libcurl4 \
  libxml2-dev \
  libxt-dev \
  openssl \
  libssl-dev \
  wget \
  curl \
  bzip2 \
  libbz2-dev \
  libpng-dev \
  libhdf5-dev \
  pigz \
  libudunits2-dev \
  libgdal-dev \
  libgeos-dev \
  libboost-dev \
  libboost-iostreams-dev \
  libboost-log-dev \
  libboost-system-dev \
  libboost-test-dev \
  libz-dev \
  libarmadillo-dev \
  libglpk-dev \
  jags \
  libgsl-dev \
  libharfbuzz-dev \
  libfribidi-dev \
  cmake \
  && apt-get clean \
  && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*






# Set global R options
RUN echo "options(repos = 'https://cloud.r-project.org')" > $(R --no-echo --no-save -e "cat(Sys.getenv('R_HOME'))")/etc/Rprofile.site
ENV RETICULATE_MINICONDA_ENABLED=FALSE

# Install Seurat's system dependencies
RUN apt-get update
RUN apt-get install -y \
    libhdf5-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libpng-dev \
    libboost-all-dev \
    libxml2-dev \
    openjdk-8-jdk \
    python3-dev \
    python3-pip \
    wget \
    git \
    libfftw3-dev \
    libgsl-dev \
    pkg-config



# Install UMAP
RUN LLVM_CONFIG=/usr/lib/llvm-10/bin/llvm-config pip3 install llvmlite
RUN pip3 install numpy
RUN pip3 install umap-learn

# Install FIt-SNE
RUN git clone --branch v1.2.1 https://github.com/KlugerLab/FIt-SNE.git
RUN g++ -std=c++11 -O3 FIt-SNE/src/sptree.cpp FIt-SNE/src/tsne.cpp FIt-SNE/src/nbodyfft.cpp  -o bin/fast_tsne -pthread -lfftw3 -lm



# Install latest Matrix
RUN R --no-echo --no-restore --no-save -e "install.packages('Matrix')"

# Install Seurat
RUN R --no-echo --no-restore --no-save -e "install.packages('remotes')"
RUN R -e "remotes::install_github('satijalab/seurat', 'seurat5', quiet = TRUE)"
RUN R -e "remotes::install_github(repo = 'satijalab/azimuth',  dependencies = TRUE)"

RUN R -e "devtools::install_github('enblacar/SCpubr')"
RUN R  -e "install.packages('scatterpie')"
RUN R  -e "install.packages('concaveman')"


RUN R -e "update.packages(ask = FALSE)"

RUN install2.r --error \
    DT \
    profvis \
    tictoc \
    markdown \
    plotly \
    bench \
    mclust 

RUN R -e "BiocManager::install(version = "3.17", ask = FALSE)"
RUN R -e "BiocManager::install('Biobase')"

RUN R -e "BiocManager::install('glmGamPoi')"

RUN R  -e "devtools::install_github('YingMa0107/CARD')"
RUN R -e "devtools::install_github('carter-allen/spruce')"
RUN R -e "devtools::install_github('carter-allen/maple')"
RUN R -e "BiocManager::install('TOAST')"
RUN R  -e "devtools::install_github('xuranw/MuSiC')"
RUN R  -e "remotes::install_github('cancerbits/DElegate')"
RUN R -e "BiocManager::install('edgeR')"


RUN R -e "BiocManager::install(c('clusterProfiler', 'msigdb'))"
RUN R  -e "devtools::install_github('wjawaid/enrichR')"

RUN R  -e "devtools::install_github('CBIIT-CGBB/scCorr')"
RUN R  -e "install.packages('corrplot')"
RUN R  -e "install.packages('DT')"
RUN R  -e "install.packages('pheatmap')"
#RUN R  -e "devtools::install_github('kassambara/ggpubr')"
#RUN R -e "BiocManager::install('scater')"
#RUN R  -e "remotes::install_github('XiaoZhangryy/iSC.MEB')"

RUN R  -e "install.packages('msigdbr')"

RUN apt-get install -y libmagick++-dev
RUN R  -e "remotes::install_github('jbergenstrahle/STUtility')"
RUN R  -e "devtools::install_github('davidgohel/ggiraph')"
RUN R  -e "remotes::install_github('ludvigla/semla')"
RUN R  -e "install.packages(c('shinyBS', 'beakr', 'colourpicker'))"
RUN R -e "devtools::install_github('kassambara/ggpubr')"
RUN R  -e "remotes::install_github('dgrun/RaceID3_StemID2_package')"
RUN R -e "BiocManager::install('TSCAN')"
RUN R -e "BiocManager::install('dittoSeq')"
RUN R -e "BiocManager::install('infercnv')"

RUN R -e "devtools::install_github('lme4/lme4',dependencies=TRUE)"
RUN R -e "devtools::install_github('hojsgaard/pbkrtest')"
RUN R -e "devtools::install_github('statistikat/VIM')"
RUN R -e "BiocManager::install('destiny')"
RUN R  -e "remotes::install_github('aet21/SCENT')"
RUN R -e "devtools::install_github('hwarden162/SCEnt')"
RUN R -e "BiocManager::install('scater')"

# Intsall SPATA2 and dependencies
RUN R -e "BiocManager::install(c('batchelor', 'Matrix.utils', 'EBImage'))"

RUN R -e "devtools::install_github('kueckelj/confuns')"
RUN R -e "devtools::install_github('theMILOlab/SPATAData')"
RUN R  -e "remotes::install_github('cvarrichio/Matrix.utils')"
RUN R  -e "install.packages('anndata')"
RUN R -e "devtools::install_github('theMILOlab/SPATA2')"

RUN R  -e "install.packages('glmnet')"
RUN R  -e "install.packages('simpleCache')"
RUN R  -e "install.packages('caret')"
RUN R  -e "install.packages('SoupX')"
RUN R  -e "remotes::install_github(repo = 'cancerbits/canceRbits')"
RUN R  -e "BiocManager::install('scRepertoire')"


RUN R  -e "remotes::install_github('satijalab/seurat-wrappers')"
RUN R -e "devtools::install_github('YuLab-SMU/ggtree')"
RUN R -e "BiocManager::install(c('clusterProfiler'))"

#RUN chmod -R a+rw ${R_HOME}/site-library # so that everyone can dynamically install more libraries within container
#RUN chmod -R a+rw ${R_HOME}/library

# add custom options for rstudio sessions
# make sure sessions stay alive forever
RUN echo "session-timeout-minutes=0" >> /etc/rstudio/rsession.conf
# make sure we get rstudio server logs in the container
# RUN echo $'[*]\nlog-level=warn\nlogger-type=file\n' > /etc/rstudio/logging.conf

# make sure all R related binaries are in PATH in case we want to call them directly
ENV PATH ${R_HOME}/bin:$PATH

