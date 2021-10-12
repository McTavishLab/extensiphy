# syntax=docker/dockerfile:1
FROM ubuntu:latest
WORKDIR /usr/src/app


RUN apt-get update && DEBIAN_FRONTEND="noninteractive" TZ="America/New_York" apt-get install -y tzdata

#install python and symlink python3 to python
RUN apt-get install -y python3
RUN apt-get install -y python3-pip
RUN ln -s /usr/bin/python3 /usr/bin/python
RUN pip3 install dendropy

#install dependencies available from apt-get
RUN apt-get install raxml
RUN apt-get install git -y
RUN apt-get install seaview -y
RUN apt-get install wget -y
RUN apt install default-jre -y

#install parsnp
RUN wget https://github.com/marbl/parsnp/releases/download/v1.2/parsnp-Linux64-v1.2.tar.gz
RUN tar -xzf parsnp-Linux64-v1.2.tar.gz
RUN cp ./Parsnp-Linux64-v1.2/parsnp /usr/local/bin/

#install spades
RUN wget http://cab.spbu.ru/files/release3.15.3/SPAdes-3.15.3-Linux.tar.gz
RUN tar -xzf SPAdes-3.15.3-Linux.tar.gz
RUN cp ./SPAdes-3.15.3-Linux/bin/spades* /usr/local/bin/

#install BBDUK and repair
RUN wget https://sourceforge.net/projects/bbmap/files/BBMap_38.94.tar.gz
RUN tar -xzf BBMap_38.94.tar.gz
ENV PATH="/usr/src/app/bbmap:/usr/src/app/SPAdes-3.15.3-Linux/bin/:/usr/src/app/Parsnp-Linux64-v1.2/:$PATH"

# install Extensiphy
RUN git clone https://github.com/McTavishLab/extensiphy.git

WORKDIR /usr/src/app/extensiphy
