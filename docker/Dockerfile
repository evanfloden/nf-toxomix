FROM continuumio/miniconda

LABEL authors="evan.floden@gmail.com"

# Install container-wide requirements gcc, pip, zlib, libssl, make, libncurses, fortran77, g++, R
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        g++ \
        gcc \
        gfortran \
        libbz2-dev \
        libcurl4-openssl-dev \
        liblzma-dev \
        libncurses5-dev \
        libpcre3-dev \
        libreadline-dev \
        libssl-dev \
        make \
        python3 \
        python-dev \
        zlib1g-dev \
        procps \
    && rm -rf /var/lib/apt/lists/*

# Install pip
RUN curl -fsSL https://bootstrap.pypa.io/get-pip.py -o /opt/get-pip.py && \
    python /opt/get-pip.py && \
    rm /opt/get-pip.py

# Install Conda Packages
RUN conda config --add channels defaults \
 && conda config --add channels bioconda \
 && conda config --add channels conda-forge \
 && conda install -y kallisto \
                    r r-essentials \
                    r-irkernel \
                    r-sleuth \
                    fastqc=0.11.5 \
                    multiqc=1.1 \
                    pandas \
                    xlrd
