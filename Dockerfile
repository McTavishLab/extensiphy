# syntax=docker/dockerfile:1
FROM ubuntu:20.04
WORKDIR /usr/src/app

#needed for samtools install
RUN apt-get update && DEBIAN_FRONTEND="noninteractive" TZ="America/New_York" apt-get install -y tzdata

#install python 3 and alias to python
RUN apt-get install -y python3
RUN apt-get install -y python3-pip
RUN ln -s /usr/bin/python3 /usr/bin/python

#intall an bunch o deps
RUN apt-get install raxml
RUN apt-get install seqtk
RUN apt-get install bcftools -y
RUN apt-get install samtools -y

#fastx toolkit and bwa need wget
RUN apt-get install gcc g++ pkg-config wget -y
RUN wget http://hannonlab.cshl.edu/fastx_toolkit/fastx_toolkit_0.0.13_binaries_Linux_2.6_amd64.tar.bz2
RUN tar -xjf fastx_toolkit_0.0.13_binaries_Linux_2.6_amd64.tar.bz2
RUN cp ./bin/* /usr/local/bin


RUN wget https://github.com/bwa-mem2/bwa-mem2/releases/download/v2.0pre2/bwa-mem2-2.0pre2_x64-linux.tar.bz2
RUN tar -xjf bwa-mem2-2.0pre2_x64-linux.tar.bz2
##Either OR!
#RUN cp bwa-mem2-2.0pre2_x64-linux/bwa-mem* /usr/local/bin/

ENV PATH="/usr/src/app/bwa-mem2-2.0pre2_x64-linux:${PATH}"

#use git to grab EP
RUN apt-get install git -y
RUN pip install dendropy
RUN git clone https://github.com/McTavishLab/extensiphy.git

WORKDIR /usr/src/app/extensiphy
RUN git checkout dev



## Speciifcs example
# docker build --tag eptest .
# docker run --name myspecificeptest -it -v /home/ejmctavish/projects/docker_demo/imaginarydata:/usr/src/app/tmpdata eptest bash
# docker start -i  myspecificeptest 



put the data you want to use in a folder named ...
##Generalized example
# docker build --tag {general container name} .
# docker run --name {name for this specific instance} -it -v {full path to dir you want shared with docker}:/usr/src/app/linked_host_dir {general container name} bash
# docker start -i  {name for this specific instance}


