FROM openjdk:8

LABEL authors="evan.floden@gmail.com" \
    description="Docker image containing all requirements for NF-toxomix pipeline"

# Install container-wide requrements gcc, pip, zlib, libssl, make, libncurses, fortran77, g++, R
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        g++ \
        gcc \
        gfortran \
        libbz2-dev \
        libcurl4-openssl-dev \
        libgsl-dev \
        libgsl2 \
        liblzma-dev \
        libncurses5-dev \
        libpcre3-dev \
        libreadline-dev \
        libssl-dev \
        make \
        python-dev \
        zlib1g-dev \
        hdf5-tools \
        libhdf5-dev \
        hdf5-helpers \
        libhdf5-serial-dev \
    && rm -rf /var/lib/apt/lists/*

# Install pip
RUN curl -fsSL https://bootstrap.pypa.io/get-pip.py -o /opt/get-pip.py && \
    python /opt/get-pip.py && \
    rm /opt/get-pip.py

RUN wget http://repo.continuum.io/miniconda/Miniconda3-3.7.0-Linux-x86_64.sh -O ~/miniconda.sh && \
    bash ~/miniconda.sh -b -p /miniconda && \
    export PATH="/miniconda/bin:$PATH"

# Install Kallisto
RUN /miniconda/bin/conda config --add channels defaults \
 && /miniconda/bin/conda config --add channels conda-forge \
 && /miniconda/bin/conda config --add channels bioconda \
 && /miniconda/bin/conda install kallisto

# Install Sleuth
RUN R -e 'source("http://bioconductor.org/biocLite.R"); library(BiocInstaller); biocLite(c("XML","biomaRt")); biocLite("rhdf5"); install.packages("devtools", repos="http://cloud.r-project.org/"); devtools::install_github("pachterlab/sleuth")'


# Install FastQC
RUN curl -fsSL http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.5.zip -o /opt/fastqc_v0.11.5.zip && \
    unzip /opt/fastqc_v0.11.5.zip -d /opt/ && \
    chmod 755 /opt/FastQC/fastqc && \
    ln -s /opt/FastQC/fastqc /usr/local/bin/fastqc && \
    rm /opt/fastqc_v0.11.5.zip

# Install MultiQC
RUN pip install git+git://github.com/ewels/MultiQC.git
