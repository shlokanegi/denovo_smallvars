FROM ubuntu:20.04

# MAINTAINER Deprecated
LABEL maintainer="shnegi@ucsc.edu"

# Prevent dpkg from trying to ask any questions, ever
ENV DEBIAN_FRONTEND noninteractive
ENV DEBCONF_NONINTERACTIVE_SEEN true

RUN apt-get -y update \
    && apt-get install -y --no-install-recommends -y tabix wget git unzip make gcc g++ ca-certificates zlib1g-dev libcurl4-openssl-dev libbz2-dev libhdf5-dev liblzma-dev python3 \
    && rm -rf /var/lib/apt/lists/*

# get rtg tools
WORKDIR /apps

RUN wget https://github.com/RealTimeGenomics/rtg-tools/releases/download/3.12.1/rtg-tools-3.12.1-linux-x64.zip && \
    unzip rtg-tools-3.12.1-linux-x64.zip && \
    /bin/sh -c echo "yes" && \
    rm -rf rtg-tools-3.12.1-linux-x64.zip

RUN /bin/sh -c echo "RTG_TALKBACK=true \nRTG_USAGE=false \n" > /apps/rtg-tools-3.12.1/rtg.cfg && \
RUN ln -s /apps/rtg-tools-3.12.1/rtg /usr/local/bin/rtg

# get bcftools
ARG bcftoolsVer="1.3.1"
RUN wget https://github.com/samtools/bcftools/releases/download/${bcftoolsVer}/bcftools-${bcftoolsVer}.tar.bz2 && \
    tar -vxjf bcftools-${bcftoolsVer}.tar.bz2 && \
    rm bcftools-${bcftoolsVer}.tar.bz2 && \
    cd bcftools-${bcftoolsVer} && \
    make && \
    make install && \
    mkdir /data

WORKDIR /data
