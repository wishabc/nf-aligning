# Build base
FROM condaforge/mambaforge:latest AS build-base
RUN apt-get update && apt-get install -y \
      bash \
      rsync \
      build-essential

###########
# Kentutils
FROM build-base as build-kentutils
RUN mkdir /scripts && rsync -azvP rsync://hgdownload.soe.ucsc.edu/genome/admin/exe/linux.x86_64/ /scripts


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


#######################
# Final image
FROM condaforge/mambaforge:latest as aligning-plus-hotspots
ENV DEBIAN_FRONTEND=noninteractive
SHELL ["/bin/bash", "-c"]
COPY ./environment.yml /environment.yml
RUN mamba env create -n babachi --file /environment.yml && echo 'conda activate babachi' >> ~/.bashrc

COPY --from=build-hotspot2 /hotspot2/bin /usr/local/bin/
COPY --from=build-hotspot2 /hotspot2/scripts /usr/local/bin/
COPY --from=build-hotspot2 /modwt/bin /usr/local/bin/
COPY --from=build-kentutils /scripts/ /usr/local/bin/