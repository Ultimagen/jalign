ARG TARGET="znver2"
ARG BUILD="all"

############################
FROM ubuntu:22.04 AS build-base

ENV TZ=US
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

RUN apt-get update && \
  apt-get -y install --no-install-recommends \
    libssl-dev \
    locales \
    build-essential \
    clang-11 \
    libboost-all-dev \
    git-lfs \
    pkg-config \
    bcftools \
    samtools \
	bedtools \
  && \
  rm -rf /var/lib/apt/lists/* && apt-get autoclean

# This is a hack to get around the apparent lack of /etc/alternatives support for clang-15
RUN ln -s clang-15 /usr/bin/clang && \
    ln -s clang++-15 /usr/bin/clang++

# locales
RUN locale-gen en_US.UTF-8
ENV LANG=en_US.UTF-8
ENV LANGUAGE=en_US:en
ENV LC_ALL=en_US.UTF-8


############################
FROM build-base AS dev

RUN apt-get update && \
  apt-get -y install --no-install-recommends \
    vim \
    ssh \
    zsh \
    wget \
    unzip \
    libparasail-dev

RUN apt-get update && \
  apt-get -y install --no-install-recommends \
    python3-pip \
	parallel
RUN ln -s /usr/bin/python3 /usr/bin/python
RUN pip install jupyterlab pysam networkx joblib


# download manta source code (no need to compule at this time)
RUN mkdir -p manta
RUN cd manta && wget --no-check-certificate https://github.com/Illumina/manta/archive/refs/tags/v1.6.0.zip && unzip v1.6.0.zip

# build executable
FROM dev AS build

# Create non root user
ARG USER_ID=1000
ARG GROUP_ID=$USER_ID
RUN groupadd -g $GROUP_ID jalign && \
    useradd -u $USER_ID -g $GROUP_ID -m -s /bin/bash jalign

# Create jalign working directory accessible by jalign user
RUN mkdir -p /jalign && \
    chown jalign:jalign /jalign && \
    chmod 755 /jalign

COPY . jalign/
RUN cd jalign/jump_align && ./mk.sh && cp jump_align /usr/local/bin && cp para_jalign /usr/local/bin

# prepare python stuff
RUN cd jalign && jupyter-nbconvert --to python cnv_realign.ipynb && chmod +x cnv_realign.py
RUN cd jalign && jupyter-nbconvert --to python sa_walker.ipynb && chmod +x sa_walker.py
RUN cd jalign && jupyter-nbconvert --to python pair_haps.ipynb && chmod +x pair_haps.py

ENV PATH="/jalign/:/jalign/jump_align:${PATH}"

USER jalign
WORKDIR /jalign

