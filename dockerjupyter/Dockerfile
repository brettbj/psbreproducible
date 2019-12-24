FROM jupyter/r-notebook
MAINTAINER Brett Beaulieu-Jones

ARG DEBIAN_FRONTEND=noninteractive
USER root

RUN apt-get update
RUN apt-get -y install r-base 
RUN apt-get -y install r-recommended
RUN apt-get -y install r-base-dev
RUN apt-get -y install wget
RUN apt-get install build-essential libcurl4-gnutls-dev libxml2-dev libssl-dev

RUN Rscript -e 'install.packages("BiocManager", repos="http://cran.us.r-project.org");BiocManager::install(c("multtest", "affy", "AnnotationDbi"))'

RUN mkdir -p /data/
RUN wget -q http://mbni.org/customcdf/18.0.0/entrezg.download/hgu133plus2hsentrezgcdf_18.0.0.tar.gz -P /data
RUN R CMD INSTALL /data/hgu133plus2hsentrezgcdf_18.0.0.tar.gz 

RUN wget -q ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE47nnn/GSE47664/suppl/GSE47664_RAW.tar -P /data
RUN mkdir -p /data/CELS
RUN tar -xvf /data/GSE47664_RAW.tar -C /data/CELS
RUN mkdir -p /plots

RUN R -e "install.packages('testthat', repos='http://cran.us.r-project.org');"
RUN R -e "install.packages('knitr', repos='http://cran.us.r-project.org');"
RUN R -e "install.packages('devtools', repos='http://cran.us.r-project.org');"
RUN R -e "install.packages('roxygen2', repos='http://cran.us.r-project.org');"

