FROM pditommaso/dkrbase:1.2

MAINTAINER Paolo Di Tommaso <paolo.ditommaso@gmail.com>

#
# Install pre-requistes
#
RUN apt-get update --fix-missing && \
  apt-get install -q -y samtools python libstdc++6 r-base
  
#
# RNA-Seq tools 
# 
#RUN wget -q -O bowtie.zip https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.2.7/bowtie2-2.2.7-linux-x86_64.zip/download && \
#  unzip bowtie.zip -d /opt/ && \
#  ln -s /opt/bowtie2-2.2.7/ /opt/bowtie && \
#  rm bowtie.zip 
#  
#RUN \
#  wget -q http://cole-trapnell-lab.github.io/cufflinks/assets/downloads/cufflinks-2.2.1.Linux_x86_64.tar.gz -O- \
#  | tar xz -C /opt/ && \
#  ln -s /opt/cufflinks-2.2.1.Linux_x86_64/ /opt/cufflinks 
#  
#  
#RUN \
#  wget -q https://ccb.jhu.edu/software/tophat/downloads/tophat-2.1.0.Linux_x86_64.tar.gz -O- \
#  | tar xz -C /opt/ && \
#  ln -s /opt/tophat-2.1.0.Linux_x86_64/ /opt/tophat 
#
  
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


#
# Finalize environment
#
ENV PATH /usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/opt/STAR:/opt/sratoolkit
