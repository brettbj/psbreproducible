FROM ubuntu:18.04
MAINTAINER Brett Beaulieu-Jones

ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get update

RUN apt-get -y install r-base 
RUN apt-get -y install r-recommended
RUN apt-get -y install r-base-dev
RUN apt-get -y install wget

RUN Rscript -e 'source("http://bioconductor.org/biocLite.R"); biocLite(c("multtest", "affy", "AnnotationDbi"))'
RUN Rscript -e 'install.packages("devtools");' 
RUN Rscript -e 'install.packages("roxygen2");'
RUN Rscript -e 'install.packages("testthat");'
RUN Rscript -e  'install.packages("knitr");'

RUN mkdir -p /data/
RUN wget -q http://mbni.org/customcdf/18.0.0/entrezg.download/hgu133plus2hsentrezgcdf_18.0.0.tar.gz -P /data
RUN R CMD INSTALL /data/hgu133plus2hsentrezgcdf_18.0.0.tar.gz 

RUN wget -q ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE47nnn/GSE47664/suppl/GSE47664_RAW.tar -P /data