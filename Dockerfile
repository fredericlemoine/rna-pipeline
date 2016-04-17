FROM ubuntu

MAINTAINER Frederic Lemoine

#
# Install pre-requistes
#
RUN apt-get update --fix-missing && \
  apt-get install -q -y samtools python libcurl4-gnutls-dev libxml2 libxml2-dev libreadline6 libreadline6-dev wget gfortran g++ gcc make libpng-dev libjpeg-dev
  

RUN \
    wget -q https://cran.r-project.org/src/base/R-3/R-3.2.5.tar.gz -O- \
    | tar xz -C /opt/ && \
    cd /opt/R-3.2.5/ && \
    ./configure --with-x=no && \
    make && \
    make install && \
    rm -rf /opt/R-3.2.5

RUN \
  wget -q https://github.com/alexdobin/STAR/archive/2.5.1b.tar.gz -O- \
  | tar xz -C /opt/ && \
  cd /opt/STAR-2.5.1b/ && \
  make && \
  mkdir /opt/STAR && \
  cp /opt/STAR-2.5.1b/source/STAR /opt/STAR/ && \
  rm -rf /opt/STAR-2.5.1b/ 

RUN \
    wget -q http://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.5.7/sratoolkit.2.5.7-ubuntu64.tar.gz -O- \
    | tar xz -C /opt/ && \
    ln -s /opt/sratoolkit.2.5.7-ubuntu64/bin opt/sratoolkit

RUN \
    Rscript -e 'source("http://bioconductor.org/biocLite.R"); biocLite("DEXSeq")'

RUN \
    chmod +x /usr/local/lib/R/library/DEXSeq/python_scripts/* && \
    ln -s /usr/local/lib/R/library/DEXSeq/python_scripts/ /opt/dexseq

RUN echo 'alias dexseq_count="python /opt/dexseq/dexseq_count.py"' >> ~/.bashrc
RUN echo 'alias dexseq_prepare_annoation="python /opt/dexseq/dexseq_prepare_annotation.py"' >> ~/.bashrc

#
# Finalize environment
#
ENV PATH /usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/opt/STAR:/opt/sratoolkit