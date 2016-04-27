FROM ubuntu

MAINTAINER Frederic Lemoine

#
# Install pre-requistes
#
RUN apt-get update --fix-missing && \
  apt-get install -q -y samtools python libcurl4-gnutls-dev libxml2 libxml2-dev libreadline6 \
                        libreadline6-dev wget gfortran g++ gcc make libpng-dev libjpeg-dev \
			libcairo2-dev python-numpy python-matplotlib python-pip

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
  mkdir /opt/STAR && \
  cp /opt/STAR-2.5.1b/bin/Linux_x86_64_static/* /opt/STAR/ && \
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

RUN apt-get install patch
COPY dexseq.patch /tmp/

RUN patch -i /tmp/dexseq.patch /usr/local/lib/R/library/DEXSeq/python_scripts/dexseq_prepare_annotation.py 

RUN echo 'alias dexseq_count="python /opt/dexseq/dexseq_count.py"' >> ~/.bashrc
RUN echo 'alias dexseq_prepare_annotation="python /opt/dexseq/dexseq_prepare_annotation.py"' >> ~/.bashrc

RUN apt-get install -q -y python-dev


RUN pip install pysam && \
    pip install HTSeq

#
# Finalize environment
#
ENV PATH /bin/:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/opt/STAR:/opt/sratoolkit

RUN apt-get remove -q -y libcurl4-gnutls-dev libxml2-dev libreadline6-dev gfortran g++ gcc make libpng-dev libjpeg-dev libcairo2-dev python-pip patch python-dev && apt-get autoremove -y
