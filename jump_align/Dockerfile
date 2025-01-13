ARG TARGET="znver2"
ARG BUILD="all"

############################
FROM ubuntu:22.04 AS build-base

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
    python2.7 \
    ssh \
    zsh \
    wget \
    unzip \
  ;

# download manta source code (no need to compule at this time)
RUN mkdir -p manta
RUN cd manta && wget --no-check-certificate https://github.com/Illumina/manta/archive/refs/tags/v1.6.0.zip && unzip v1.6.0.zip

# build executable
FROM dev AS build
RUN mkdir -p jump_align
COPY . jump_align/
RUN cd jump_align && ./mk.sh && cp jump_align /usr/local/bin

