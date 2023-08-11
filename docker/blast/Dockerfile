# Base image for blast/2.12.0,
# uses Ubuntu Jammy (LTS)
FROM ubuntu:22.04

# Depedencies of ncbi-blast+/2.12.0:
#  - libbz2-1.0: high-quality block-sorting file compressor library
#  - libc6 (>= 2.34): GNU C Library, shared libraries
#  - libgcc-s1 (>= 3.3.1): GCC support library
#  - libgcc-s1 (>= 3.4) [amd64]
#  - libgomp1 (>= 6): GCC OpenMP (GOMP) support library
#  - liblmdb0 (>= 0.9.7): Lightning memory-mapped database shared library
#  - libmbedcrypto7 (>= 2.28.0): lightweight crypto and SSL/TLS library
#  - libmbedtls14 (>= 2.28.0): lightweight crypto and SSL/TLS library
#  - libpcre3: Old Perl 5 Compatible Regular Expression Library
#  - libstdc++6 (>= 12): GNU Standard C++ Library v3
#  - ncbi-data: Platform-independent data for the NCBI toolkit
#  - perl: Larry Wall's Practical Extraction and Report Language
#  - python3: Interactive high-level object-oriented language
#  - zlib1g (>= 1:1.2.11): compression library 
#  - t-coffee: Multiple Sequence Alignment

LABEL maintainer="Skyler Kuhn <kuhnsa AT nih.gov>" \
    base_image="ubuntu:22.04" \
    version="v0.1.0"   \
    software="blast" \
    about.summary="Docker image for NCBI Blast+" \
    about.home="https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/" \
    about.documentation="https://blast.ncbi.nlm.nih.gov/Blast.cgi" \
    about.tags="Genomics"

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
        apt-transport-https \
        autotools-dev \
        automake \
        autoconf \
        build-essential \
        ca-certificates \
        cpanminus \
        curl \
        gawk \
        git \
        grep \
        gzip \
        locales \
        unzip \
        wget \
        # Installs ncbi-blast+/2.12.0
        ncbi-blast+=2.12.0+ds-3build1 \
        python3 \
        python3-pip \
        unzip \
        zlib1g-dev \
    && apt-get clean && apt-get purge \
    && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# Set the locale
RUN localedef -i en_US -f UTF-8 en_US.UTF-8

# Make python3 the default interpreter
# and install Python Packages 
# and install multiqc/1.14
RUN ln -sf /usr/bin/python3 /usr/bin/python
RUN pip3 install --upgrade pip \
	&& pip3 install argparse \
    && pip3 install matplotlib \
    && pip3 install multiqc==1.14 \
    && pip3 install pandas \
    && pip3 install scipy

# Perl fix issue
RUN cpanm FindBin Term::ReadLine LWP::Simple

############### INSTALL ################
# Install any bioinformatics tools
# available with apt-get on Ubuntu/22.04

############### MANUAL ################
# Install tools from src manually

################ POST #################
# Add Dockerfile and export environment 
# variables and update permissions
ADD Dockerfile /opt2/blast.dockerfile
RUN chmod -R a+rX /opt2
ENV PATH="/opt2:$PATH"
WORKDIR /data2