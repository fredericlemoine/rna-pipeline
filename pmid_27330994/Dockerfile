FROM ubuntu

MAINTAINER Frederic Lemoine

#
# Install pre-requistes
#
RUN apt-get update --fix-missing && \
  apt-get install -q -y libcurl4-gnutls-dev libxml2 libxml2-dev libreadline6 \
                        libreadline6-dev wget gfortran g++ gcc make libpng-dev libjpeg-dev \
			libcairo2-dev

RUN \
    wget -q https://cran.r-project.org/src/base/R-3/R-3.2.5.tar.gz -O- \
    | tar xz -C /opt/ && \
    cd /opt/R-3.2.5/ && \
    ./configure --with-x=no && \
    make && \
    make install && \
    rm -rf /opt/R-3.2.5 

RUN apt-get install unzip

RUN \
  wget "http://downloads.sourceforge.net/project/bowtie-bio/bowtie2/2.2.8/bowtie2-2.2.8-linux-x86_64.zip?r=https%3A%2F%2Fsourceforge.net%2Fprojects%2Fbowtie-bio%2Ffiles%2Fbowtie2%2F2.2.8%2F&ts=1461156134&use_mirror=ufpr" && \
  mv bowtie2-2.2.* bowtie2-2.2.8-linux-x86_64.zip && \
  unzip bowtie2-2.2.8-linux-x86_64.zip && \
  mkdir /opt/bowtie2/ && \
  cp bowtie2-2.2.8/bowtie2* /opt/bowtie2 && \
  rm -rf bowtie2-2.2.8/ bowtie2-2.2.8-linux-x86_64.zip

RUN \
  wget "http://downloads.sourceforge.net/project/subread/subread-1.5.0-p2/subread-1.5.0-p2-source.tar.gz?r=https%3A%2F%2Fsourceforge.net%2Fprojects%2Fsubread%2Ffiles%2Fsubread-1.5.0-p2%2F&ts=1461157751&use_mirror=iweb" && \
  mv subread* subread-1.5.0-p2-source.tar.gz && \
  tar -xzvf subread-1.5.0-p2-source.tar.gz && \
  cd subread-1.5.0-p2-source/src && \
  make -f Makefile.Linux && \
  mkdir /opt/subread/ && \
  cp ../bin/featureCounts /opt/subread && \
  cd ../ && \
  rm -rf subread-1.5.0-p2-source/ subread-1.5.0-p2-source.tar.gz

RUN apt-get install -y libncurses5-dev
RUN \
  wget http://github.com/samtools/samtools/releases/download/1.3.1/samtools-1.3.1.tar.bz2 && \
  tar -xjvf samtools-1.3.1.tar.bz2 && \
  cd samtools-1.3.1 && \
  ./configure && \
  make && \
  mkdir /opt/samtools/ && \
  mv samtools /opt/samtools/ && \
  cd .. && \
  rm -rf samtools-1.3.1 samtools-1.3.1.tar.bz2 

RUN \
    wget -q http://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.5.7/sratoolkit.2.5.7-ubuntu64.tar.gz -O- \
    | tar xz -C /opt/ && \
    ln -s /opt/sratoolkit.2.5.7-ubuntu64/bin opt/sratoolkit

RUN \
    Rscript -e 'source("http://bioconductor.org/biocLite.R"); biocLite("DESeq2")'

RUN \
    wget ftp://emboss.open-bio.org/pub/EMBOSS/EMBOSS-6.6.0.tar.gz && \
    tar -xzvf EMBOSS-6.6.0.tar.gz && \
    cd EMBOSS-6.6.0 && \
    ./configure --without-x --prefix=/opt/emboss && \
    make -j 10 && \
    mkdir /opt/emboss && \
    mkdir /opt/emboss/acd &&\
    mkdir /opt/emboss/data &&\
    cp ./ajax/acd/.libs/libacd.so.6 nucleus/.libs/libnucleus.so.6  ./ajax/ajaxdb/.libs/libajaxdb.so.6 ./ajax/ensembl/.libs/libensembl.so.6  ./ajax/graphics/.libs/libajaxg.so.6 ./ajax/zlib/.libs/libezlib.so.1 ./ajax/expat/.libs/libeexpat.so.1 ./ajax/pcre/.libs/libepcre.so.7 ./plplot/.libs/libeplplot.so.3 ./ajax/core/.libs/libajax.so.6 /opt/emboss/ && \
    cp ./emboss/acd/* /opt/emboss/acd && \
    cp ./emboss/data/Etags* /opt/emboss/data && \
    cp ./emboss/data/Efeatures* /opt/emboss/data && \
    cp emboss/.libs/seqret /opt/emboss && \
    cd .. && \
    rm -rf EMBOSS-6.6.0 EMBOSS-6.6.0.tar.gz

ENV EMBOSS_ACDROOT /opt/emboss/acd
ENV EMBOSS_ROOT /opt/emboss/

RUN apt-get install -y python
#
# Finalize environment
#
ENV PATH /bin/:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/opt/STAR:/opt/bowtie2:/opt/sratoolkit:/opt/samtools:/opt/subread:/opt/emboss
ENV LD_LIBRARY_PATH /opt/emboss

RUN apt-get remove -q -y unzip libcurl4-gnutls-dev libncurses5-dev libxml2-dev libreadline6-dev gfortran g++ gcc make libpng-dev libjpeg-dev libcairo2-dev  && apt-get autoremove -y
RUN apt-get install -y libgfortran3 libgomp1 libcairo2 libjpeg8
