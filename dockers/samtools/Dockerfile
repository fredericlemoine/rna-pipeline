FROM ubuntu

MAINTAINER Frederic Lemoine

#
# Install pre-requistes
#
RUN apt-get update --fix-missing
RUN apt-get install -q -y  wget g++ gcc make bzip2 libncurses5-dev zlib1g-dev

RUN \
  wget http://github.com/samtools/samtools/releases/download/1.3.1/samtools-1.3.1.tar.bz2 \
  && tar -xjvf samtools-1.3.1.tar.bz2 \
  && cd samtools-1.3.1 \
  && ./configure \
  && make \
  && mkdir /opt/samtools/ \
  && mv samtools /opt/samtools/ \
  && cd .. \
  && rm -rf samtools-1.3.1 samtools-1.3.1.tar.bz2
  
RUN apt-get remove -q -y zlib1g-dev libncurses5-dev wget g++ gcc make bzip2 \
    && apt-get autoremove -y

ENV PATH /bin/:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/opt/samtools/
