# Making Docker Conatiner
# sudo docker build --no-cache -f seurat.Dockerfile -t seurat .
# docker run -v $PWD:/mnt --name seurat -d -i -t seurat /bin/bash
# docker exec -it seurat bash
#
# Pushing to Docker Hub
# docker commit af35c94f4420 seurat:5.2.0
# docker tag 7de01c30b1f8 swingett/seurat:5.2.0
# docker push swingett/seurat:5.2.0


FROM ubuntu:plucky-20241213

LABEL image.author.name "Steven Wingett"

SHELL ["/bin/bash", "-c"]

ENV DEBIAN_FRONTEND=noninteractive 

RUN apt update -y

RUN apt install r-base  -y
RUN apt install r-base-dev -y
RUN apt install curl -y
RUN apt install apt-show-versions -y
RUN apt install libssl-dev -y
RUN apt install libcurl4-openssl-dev -y
RUN apt install libxml2-dev -y
RUN apt install libfontconfig1-dev -y 
RUN apt install libharfbuzz-dev -y
RUN apt install libfribidi-dev -y
RUN apt install libfreetype6-dev -y
RUN apt install libpng-dev -y
RUN apt install libtiff5-dev -y
RUN apt install libjpeg-dev -y
RUN apt install libcairo2-dev -y
RUN apt install libmagick++-dev -y


# library(devtools)
RUN Rscript -e 'install.packages("devtools")'

# library(tidyverse)
RUN Rscript -e 'install.packages("xml2")'
RUN Rscript -e 'install.packages("systemfonts")'
RUN Rscript -e 'install.packages("textshaping")'
RUN Rscript -e 'install.packages("ragg")'
RUN Rscript -e 'install.packages("rvest")'
RUN Rscript -e 'install.packages("tidyverse")'


# library(Matrix)
RUN Rscript -e 'install.packages("Matrix")'

# library(Seurat)
RUN Rscript -e 'install.packages("Seurat")'

