# Base image for metavirs,
# uses Ubuntu Focal (LTS)
FROM ubuntu:20.04

# Depedencies of metvirs:
#  - CAT/5.2.3            # from src, installed: 5.2.3
#  - cutadapt/1.18        # apt-get, default: 2.8
#  - bowtie/2-2.3.4.1     # apt-get, deafult: 2.3.5.1
#  - fastq_screen/0.14.1  # from src, installed: 0.15.2
#  - fastqc/0.11.5        # apt-get, deafault: 0.11.9
#  - kraken/2.1.1         # from src, installed: 2.1.2
#  - kronatools/2.8       # from src, installed: 2.8.1
#  - megahit/1.2.9        # from src, installed: 1.2.9
#  - quast/5.0.2          # from src, installed: 5.0.2
#  - multiqc/1.12         # from src, installed: 1.12
#  - python/3.8           # apt-get, default: 3.8.2
#  - samtools/1.15        # apt-get, default: 1.10
#  - seqtk/1.3            # apt-get, deafult: 1.3
#  - spades/3.15.4        # from src, installed: 3.15.4
#  - ucsc/418             # personal mirror: snapshot Nov 2nd, 2020
LABEL maintainer=kuhnsa@nih.gov

############### INIT ################
# Create Container filesystem specific 
# working directory and opt directories
# to avoid collisions with the host's
# filesystem, i.e. /opt and /data
RUN mkdir -p /opt2 && mkdir -p /data2
WORKDIR /opt2 

# Set time zone to US east coast 
ENV TZ=America/New_York
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime \
        && echo $TZ > /etc/timezone

############### SETUP ################
# This section installs system packages 
# required for your project. If you need 
# extra system packages add them here.
RUN apt-get update \
    && apt-get -y upgrade \
    && DEBIAN_FRONTEND=noninteractive apt-get install -y \
        bc \
        build-essential \
        bzip2 \
        cmake \
        cpanminus \
        curl \
        gawk \
        git \
        g++ \
        gcc \
        gfortran \
        gzip \
        locales \
        make \
        openjdk-8-jdk \
        pandoc \
        # perl/5.30.0
        perl \
        libgd-gd2-perl \
        parallel \
        pigz \
        # python/2.7.17
        python2 \
        # python/3.8.2
        python3 \
        python3-pip \
        samtools \
        unzip \
        wget \
        zlib1g-dev \
        # metaquast dependencies below
        pkg-config \
        libfreetype6-dev \
        libpng-dev \
        libboost-all-dev \
        # spades dependencies below
        bamtools \
        bwa \
        libhat-trie0 \
        libnlopt0 \
        libssw0 \
    && apt-get clean && apt-get purge \
    && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# Make python3 the default interpreter
# and install Python Packages, metaquast 
# requires matplotlib python package, 
# install multiqc/1.12 (latest), spades
# requires joblib, pyyaml (yaml)
RUN ln -sf /usr/bin/python3.8 /usr/bin/python
RUN pip3 install --upgrade pip \
	&& pip3 install argparse \
	&& pip3 install numpy \
        && pip3 install pysam \
	&& pip3 install scipy \
        && pip3 install matplotlib \
        && pip3 install setuptools \
        && pip3 install multiqc==1.12 \
        && pip3 install pyyaml \
        && pip3 install joblib

# Installl Perl Libraries, GD::Graph::bars 
# is required for fastq_screen, metaquast
# requires Time::HiRes, NOTE: apt-get the 
# package libgd-gd2-perl prior to installing
# GD::Graph::bars library due to a known bug
RUN cpanm FindBin Term::ReadLine \
        && cpanm GD::Graph::bars \
        && cpanm Time::HiRes

# Set the locale
RUN localedef -i en_US -f UTF-8 en_US.UTF-8

############### INSTALL ################
# Install any bioinformatics tools
# available with apt-get on Ubuntu/20.04
RUN apt-get update \
    && apt-get -y upgrade \
    && DEBIAN_FRONTEND=noninteractive apt-get install -y \
        # cutadapt/2.8
        cutadapt \
        # bedtools/2.27.1
        bedtools \
        # bowtie/1.2.3
        bowtie \
        # bowtie2/2.3.5.1
        bowtie2 \
        # fastqc/0.11.9
        fastqc \
        # samtools/1.10
        samtools \
        # seqtk/1.3
        seqtk \
    && apt-get clean && apt-get purge \
    && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

############### MANUAL ################
# Install fastq_screen/0.15.2 manually, 
# not available to apt-get on Ubuntu/20.04.
# Dependencies already satisfied:
# https://github.com/StevenWingett/FastQ-Screen/tree/master/docs#installation
RUN wget https://github.com/StevenWingett/FastQ-Screen/archive/refs/tags/v0.15.2.tar.gz \
        && tar -xvzf /opt2/v0.15.2.tar.gz \
        && rm /opt2/v0.15.2.tar.gz
# Add fastq_screen to $PATH
ENV PATH="/opt2/FastQ-Screen-0.15.2:$PATH"
WORKDIR /opt2

# Install kraken/v2.1.2 (latest) manually,
# not available to apt-get on Ubuntu/20.04.
# Dependencies already satisfied:
# https://github.com/DerrickWood/kraken2/releases 
RUN wget https://github.com/DerrickWood/kraken2/archive/v2.1.2.tar.gz \
	    && tar -xvzf v2.1.2.tar.gz \
	    && rm /opt2/v2.1.2.tar.gz \
	    && cd /opt2/kraken2-2.1.2 \
 	    && ./install_kraken2.sh . \
	    && chmod a+rX /opt2/kraken2-2.1.2/*
# Add kraken2 to PATH
ENV PATH="/opt2/kraken2-2.1.2:$PATH"
WORKDIR /opt2

# Install kronatools/v2.8.1) manually, 
# not available to apt-get on Ubuntu/20.04.
# Dependencies already satisfied:
# https://github.com/marbl/Krona/releases
# Build Krona Taxonomy Databases
# ktImportTaxonomy relies on NCBI taxonomy: 
# https://github.com/marbl/Krona/wiki/Installing 
RUN wget https://github.com/marbl/Krona/releases/download/v2.8.1/KronaTools-2.8.1.tar \
	    && tar -xvf KronaTools-2.8.1.tar  \
	    && rm KronaTools-2.8.1.tar \ 
	    && cd KronaTools-2.8.1/ \
	    && ./install.pl  \
	    && mkdir -p /opt2/KronaTools-2.8.1/taxonomy \
	    && ./updateTaxonomy.sh
WORKDIR /opt2

# Install UCSC tools from personal mirror,
# Snapshot was taken on November 2nd, 2020
# https://github.com/skchronicles/ucsc_mirror/
RUN wget https://github.com/skchronicles/ucsc_mirror/raw/main/faSplit \
        && chmod a+x faSplit

# Install megahit/1.2.9 (latest) manually,
# not available to apt-get on Ubuntu/20.04.
# Dependencies already satisfied:
# https://github.com/voutcn/megahit
RUN wget https://github.com/voutcn/megahit/releases/download/v1.2.9/MEGAHIT-1.2.9-Linux-x86_64-static.tar.gz \
        && tar -xvzf MEGAHIT-1.2.9-Linux-x86_64-static.tar.gz \
        && rm /opt2/MEGAHIT-1.2.9-Linux-x86_64-static.tar.gz
# Add megahit to $PATH
ENV PATH="/opt2/MEGAHIT-1.2.9-Linux-x86_64-static/bin:$PATH"
WORKDIR /opt2

# Install metaquast/5.0.2 manually, 
# not available to apt-get on Ubuntu/20.04.
# Dependencies already satisfied:
# https://github.com/ablab/quast
RUN wget https://github.com/ablab/quast/releases/download/quast_5.0.2/quast-5.0.2.tar.gz \
        && tar -xvzf quast-5.0.2.tar.gz \
        && rm /opt2/quast-5.0.2.tar.gz \
        && cd /opt2/quast-5.0.2/ \
        && ./setup.py install
WORKDIR /opt2

# Install CAT/5.2.3 manually,
# not available to apt-get on Ubuntu/20.04.
# CAT has two dependencies that have 
# not been satisfied:
#  - prodigal
#  - diamond (bundled with CAT Prep ref)
# @Install prodigal/2.6.3
RUN wget https://github.com/hyattpd/Prodigal/releases/download/v2.6.3/prodigal.linux \
        && chmod +x /opt2/prodigal.linux \
        && mv /opt2/prodigal.linux /opt2/prodigal
# @Install CAT/5.2.3
RUN wget https://github.com/dutilh/CAT/archive/refs/tags/v5.2.3.tar.gz \
        && tar -xvzf v5.2.3.tar.gz \
        && rm /opt2/v5.2.3.tar.gz
# Add CAT to $PATH
ENV PATH="/opt2/CAT-5.2.3/CAT_pack:$PATH"
WORKDIR /opt2

# Install spades/3.15.4 manually, latest version
# not available to apt-get on Ubuntu/20.04.
# Dependencies already satisfied:
RUN wget http://cab.spbu.ru/files/release3.15.4/SPAdes-3.15.4-Linux.tar.gz \
        && tar -xvzf SPAdes-3.15.4-Linux.tar.gz  \
        && rm /opt2/SPAdes-3.15.4-Linux.tar.gz 
# Add metaspades to $PATH
ENV PATH="/opt2/SPAdes-3.15.4-Linux/bin:$PATH"
WORKDIR /opt2

# Add Dockerfile and export environment 
# variables and set java8 as default with links
# to alternative versions
ADD Dockerfile /opt2/metavirs.dockerfile
RUN chmod -R a+rX /opt2 \
        && ln -s /usr/lib/jvm/java-8-openjdk-amd64/bin/java /usr/bin/java8 \
        && ln -fs /usr/lib/jvm/java-8-openjdk-amd64/bin/java /usr/bin/java
ENV PATH="/opt2:$PATH"
WORKDIR /data2