FROM ubuntu:18.04
MAINTAINER Brett Beaulieu-Jones

ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get update

# install R 
RUN apt-get -y install r-base 
RUN apt-get -y install r-recommended
RUN apt-get -y install r-base-dev
RUN apt-get -y install wget
RUN apt-get install build-essential libcurl4-gnutls-dev libxml2-dev libssl-dev

# install bioconductor 
RUN Rscript -e 'source("http://bioconductor.org/biocLite.R");biocLite(c("multtest", "affy", "AnnotationDbi"))'

# download version 18 of customcdf - packaged at time of building container
RUN mkdir -p /data/
RUN wget -q http://mbni.org/customcdf/18.0.0/entrezg.download/hgu133plus2hsentrezgcdf_18.0.0.tar.gz -P /data
RUN R CMD INSTALL /data/hgu133plus2hsentrezgcdf_18.0.0.tar.gz 

# download HeLa K/O gene expression
RUN wget -q ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE47nnn/GSE47664/suppl/GSE47664_RAW.tar -P /data
RUN mkdir -p /data/CELS
RUN tar -xvf /data/GSE47664_RAW.tar -C /data/CELS
RUN mkdir -p /plots

# install additional R packages
RUN R -e "install.packages('testthat');"
RUN R -e "install.packages('knitr');"
RUN R -e "install.packages('devtools');"
RUN R -e "install.packages('roxygen2');"

