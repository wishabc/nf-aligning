# Build base
FROM condaforge/mambaforge:latest AS build-base
RUN apt-get update && apt-get install -y \
      bash \
      rsync \
      build-essential

###########
# Kentutils
# Build from source to pin version to v438 
FROM ubuntu:18.04 AS build-kentutils
ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y \
      build-essential \
      rsync \
      libpng-dev \
      uuid-dev \
      libmysqlclient-dev \
      zlib1g-dev
RUN rsync -a rsync://hgdownload.soe.ucsc.edu/genome/admin/exe/userApps.archive/userApps.v438.src.tgz . \
      && tar xzf userApps.v438.src.tgz \
      && cd userApps \
      && make utils MACHTYPE=x86_64

##########
# Hotspot1
FROM build-base AS build-hotspot1
RUN apt-get install -y \
      build-essential \
      git \
      libgsl-dev \
      wget
RUN git clone https://github.com/StamLab/hotspot.git \
      && cd hotspot \
      && git checkout v4.1.1 \
      && cd hotspot-distr/hotspot-deploy \
      && make

##########
# Hotspot2
FROM build-base as build-hotspot2
RUN git clone https://github.com/Altius/hotspot2.git \
  && cd hotspot2 \
  && make \
  && cd / \
  && git clone https://github.com/StamLab/modwt.git \
  && cd modwt \
  && git checkout 28e9f479c737836ffc870199f2468e30659ab38d \
  && make

FROM build-base as build-conda
COPY ./environment.yml /environment.yml
RUN --mount=type=cache,target=/opt/conda/pkgs mamba env create -n babachi --file /environment.yml && echo 'conda activate babachi' >> ~/.bashrc
SHELL ["conda", "run", "--no-capture-output", "-n", "babachi", "/bin/bash", "-c"]

#######################
# Final image
FROM ubuntu:18.04 AS aligning-plus-hotspots
ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y \
      bash \
      bc \
      libcurl4 \
      libkrb5-dev \
      libgsl-dev \
      python-minimal

ENV HOTSPOT_DIR=/hotspot1
COPY --from=build-hotspot1 /hotspot/hotspot-distr/ $HOTSPOT_DIR
ENV PATH=$HOTSPOT_DIR/hotspot-deploy/bin:$PATH
COPY --from=build-hotspot2 /hotspot2/bin /usr/local/bin/
COPY --from=build-hotspot2 /hotspot2/scripts /usr/local/bin/
COPY --from=build-hotspot2 /modwt/bin /usr/local/bin/
COPY --from=build-kentutils /userApps/bin/ /usr/local/bin/
COPY --from=build-conda /opt/conda /opt/conda
ENV PATH=/opt/conda/envs/babachi/bin:$PATH
ARG PATH=/opt/conda/envs/babachi/bin:$PATH
