#sudo docker build --no-cache -f pbp.Dockerfile -t pbp .
#docker run -v $PWD:/mnt --name pbp -d -i -t pbp /bin/bash
#docker exec -it pbp bash


FROM ubuntu:noble-20241011

#FROM conda/miniconda3

LABEL image.author.name "Steven Wingett"
#LABEL image.author.email "your@email.here"

SHELL ["/bin/bash", "-c"]

#RUN apt update

RUN apt update
RUN apt install -y lftp
RUN apt install -y zip
RUN apt install -y wget

RUN mkdir -p /opt/myfiles
WORKDIR "/opt/myfiles"

RUN lftp -e 'cd /pub/swingett/software_download/; mget PBP.1.4.0.zip; bye' -u anonymous,'' ftp.mrc-lmb.cam.ac.uk
RUN unzip PBP.1.4.0.zip
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-py311_23.5.2-0-Linux-x86_64.sh
RUN bash Miniconda3-py311_23.5.2-0-Linux-x86_64.sh -b -p /opt/miniconda3/



RUN lftp -e 'cd /pub/swingett/software_download/; mget conda_bashrc; bye' -u anonymous,'' ftp.mrc-lmb.cam.ac.uk






ENV CONDA_DIR /opt/miniconda3
# Put conda in path so we can use conda activate
ENV PATH=$CONDA_DIR/bin:$PATH


#RUN echo HI
RUN conda create -n spipe python=3.10 -y
RUN conda init bash

RUN source ./conda_bashrc && conda activate spipe
#RUN conda init bash && source ~/.bashrc && conda activate spipe


#
#RUN conda activate spipe
#RUN echo "conda activate myenv" >> ~/.bashrc

WORKDIR "/opt/myfiles/ParseBiosciences-Pipeline.1.4.0"



RUN bash ./install_dependencies_conda.sh -i -y
RUN pip install --no-cache-dir ./
RUN split-pipe --version


CMD echo hi there!