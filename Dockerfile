# Build base
FROM ubuntu:18.04 as build-base
RUN apt-get update
RUN apt-get install -y \
      build-essential \
      git \
      curl \
      wget \
      zlib1g-dev

FROM condaforge/mambaforge:latest as conda
RUN echo 'Building conda'
###########
# miniconda
SHELL ["/bin/bash", "-c"]
COPY ./environment.yml /scripts/environment.yml
RUN mamba env create -n babachi --file /scripts/environment.yml
RUN echo 'conda activate babachi' >> ~/.bashrc

###########
# Kentutils
# FROM build-base as build-kentutils
# RUN apt-get install -y \
#       build-essential \
#       git \
#       libmysqlclient-dev \
#       libpng-dev \
#       libssh-dev \
#       wget \
#       zlib1g-dev
# RUN wget --quiet https://github.com/ENCODE-DCC/kentUtils/archive/v302.0.0.tar.gz \
#       && tar xf v302.0.0.tar.gz \
#       && cd kentUtils-302.0.0 \
#       && make

##########
# Hotspot2
# FROM build-base as build-hotspot2
# RUN git clone https://github.com/Altius/hotspot2.git \
#   && cd hotspot2 \
#   && make \
#   && cd / \
#   && git clone https://github.com/StamLab/modwt.git \
#   && cd modwt \
#   && git checkout 28e9f479c737836ffc870199f2468e30659ab38d \
#   && make


#######################
# Final image for DNase
# FROM ubuntu:18.04 as aligning-dnase

# ENV DEBIAN_FRONTEND=noninteractive
# RUN apt-get update
# RUN apt-get install -y \
#       bash \
#       build-essential \
#       coreutils \
#       gawk \
#       libboost-dev \
#       libgsl-dev \
#       littler \
#       openjdk-8-jre \
#       zlib1g-dev




# COPY --from=build-hotspot2 /hotspot2/bin /usr/local/bin/
# COPY --from=build-hotspot2 /hotspot2/scripts /usr/local/bin/
# COPY --from=build-hotspot2 /modwt/bin /usr/local/bin/
#COPY --from=build-kentutils /kentUtils-302.0.0/bin/ /usr/local/bin/