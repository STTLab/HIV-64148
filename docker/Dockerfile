FROM ubuntu:22.04
ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update
RUN apt-get upgrade -y
RUN apt-get install -y apt-utils
RUN apt-get install -y build-essential
RUN apt-get install -y cmake
RUN apt-get install -y zlib1g-dev
RUN apt-get install -y libboost-all-dev
RUN apt-get install -y libcurl4-openssl-dev
RUN apt-get install -y wget
RUN apt-get install -y unzip
RUN apt-get install -y nano
RUN apt-get install -y python3-pip python-is-python3

# Install Git
RUN apt-get install -y git

# Install samtools
RUN apt-get install -y samtools

# Install daccord
WORKDIR /opt
RUN ["wget", "https://github.com/gt1/daccord/releases/download/0.0.10-release-20170526170720/daccord-0.0.10-release-20170526170720-x86_64-etch-linux-gnu.tar.gz"]
RUN ["tar", "-zvxf", "daccord-0.0.10-release-20170526170720-x86_64-etch-linux-gnu.tar.gz"]
RUN ["rm", "daccord-0.0.10-release-20170526170720-x86_64-etch-linux-gnu.tar.gz"]
RUN ["ln", "-fs", "/opt/daccord-0.0.10-release-20170526170720-x86_64-etch-linux-gnu/bin/daccord", "/usr/local/bin/daccord"]

# Install Minimap2
WORKDIR /opt
RUN ["git", "clone", "https://github.com/lh3/minimap2"]
WORKDIR /opt/minimap2
RUN make
RUN ["ln", "-fs", "/opt/minimap2/minimap2", "/usr/local/bin/minimap2"]

# Install spoa
WORKDIR /opt
RUN ["git", "clone", "https://github.com/rvaser/spoa"]
WORKDIR /opt/spoa/build
RUN ["cmake", "-DCMAKE_BUILD_TYPE=Release", ".."]
RUN make
RUN ["ln", "-fs", "/opt/spoa/build/bin/spoa", "/usr/local/bin/spoa"]

# Install MetaBAT
WORKDIR /opt
RUN ["wget", "https://bitbucket.org/berkeleylab/metabat/get/master.tar.gz"]
WORKDIR /opt/berkeleylab-metabat
WORKDIR /opt
RUN ["tar", "-zxf", "master.tar.gz", "-C", "berkeleylab-metabat", "--strip-components", "1"]
RUN ["rm", "master.tar.gz"]
RUN ["cmake", "/opt/berkeleylab-metabat"]
RUN make
RUN make test
RUN make install

# Install DAZZ_DB
WORKDIR /opt
RUN ["git", "clone", "https://github.com/thegenemyers/DAZZ_DB.git"]
RUN make -C /opt/DAZZ_DB
RUN make install -C /opt/DAZZ_DB DEST_DIR=/usr/local/bin

# Install DALIGNER
WORKDIR /opt
RUN ["git", "clone", "https://github.com/thegenemyers/DALIGNER.git"]
RUN make -C DALIGNER
RUN make install -C DALIGNER DEST_DIR=/usr/local/bin

# Install Strainline
WORKDIR /opt
RUN ["git", "clone", "https://github.com/HaploKit/Strainline.git"]
RUN ["ln", "-fs", "/opt/Strainline/src/*", "/usr/local/bin/"]


ENV MAMBA_ROOT_PREFIX="/micromamba/bin/"
ENV PATH="${PATH}:/usr/local/bin/:${MAMBA_ROOT_PREFIX}"

# Install mamba
WORKDIR /micromamba
RUN apt-get install -y curl
RUN curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/1.4.1 | tar -xvj bin/micromamba
# RUN ["./bin/micromamba", "shell", "init", "-s", "bash", "-p", "/micromamba"]
RUN ["micromamba", "--version"]

# Install Snippy
RUN ["micromamba", "create", "-n", "venv"]
RUN ["micromamba", "install", "-n", "venv", "-c", "bioconda", "-c", "conda-forge", "-y", "snippy==4.6.0"]
RUN ["micromamba", "install", "-n", "venv", "-c", "bioconda", "-y", "snpeff==4.5covid19"]

# Install BLAST+
RUN ["micromamba", "install", "-n", "venv", "-c", "bioconda", "-y","blast"]

# Install Canu and Medaka
# RUN ["micromamba ", "install", "-n", "venv", "-c", "bioconda", "-c", "conda-forge", "-y", "canu", "medaka"]
RUN ["micromamba", "clean", "--all", "-y"]

RUN ["pip","install", "biopython", "pandas", "numpy"]
RUN ["pip", "install", "--upgrade", "numpy"]

RUN apt-get remove -y build-essential

COPY ../scripts /hiv64148/scripts
COPY ../settings /hiv64148/settings

WORKDIR /workspace
CMD ["/bin/bash"]
