FROM clarkad/rdrd_cosia1.5:1.0.2

# location for mounting
#RUN mkdir /home/rstudio/data

##----------------------------------------------------------------------------##
##  R PACKAGES
##----------------------------------------------------------------------------##

#### General CRAN Packages ####
# #### Project Specific CRAN Packages ####

RUN R -e 'install.packages(c("dendextend","NMF", "SoupX", "fgsea", "gprofiler2", "circlize", "msigdbr", "enrichR", "lattice", "grr", "rafalib", "markdown", "reticulate", "fields", "KernSmooth", "ROCR", "paralell", "pheatmap", "ggplot2", "Matrix.utils", "readr", "dplyr", "igraph", "reshape2", "visNetwork", "RColorBrewer", "plotly", "data.table", "cowplot", "magrittr", "Matrix", "purrr", "tibble", "png", "hdf5r"))'

#### Project Specific GitHub Packages ####

RUN R -e 'remotes::install_github(repo = "netZoo/netZooR", ref = "master")'
RUN R -e 'remotes::install_github("llrs/topGO")'

#### Project Specific Bioconductor Packages ####

RUN R -e 'BiocManager::install(c("GOSemSim", "clusterProfiler" ,"DOSE", "enrichplot", "decoupleR", "GEOquery", "rrvgo", "AnnotationHub", "org.Mm.eg.db", "ComplexHeatmap", "scran", "limma", "Biobase", "scater", "edgeR", "S4vectors", "apeglm", "DESeq2"))'

##----------------------------------------------------------------------------##
## SYSTEM PACKAGES
##----------------------------------------------------------------------------##

RUN apt-get update --yes && \
  apt-get upgrade --yes

RUN apt-get install --yes \
  python-dev-is-python3 \
  python3-pip \
  python3-dev \
  python3-pip \
  python3-tk


RUN apt-get install --yes \
  libcurl4-openssl-dev \
  libssl-dev \
  libcairo2-dev libxt-dev \
  libxml2-dev \
  libudunits2-dev \
  libhdf5-dev \
  libv8-dev \
  libgdal-dev \
  xorg libx11-dev libglu1-mesa-dev \
  samtools libboost-all-dev \
  libgsl-dev

##----------------------------------------------------------------------------##
## PYTHON PACKAGES
##----------------------------------------------------------------------------##


