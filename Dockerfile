FROM nvidia/cuda:11.7.1-cudnn8-devel-ubuntu20.04

ENV DEBIAN_FRONTEND noninteractive

RUN apt-get update
RUN apt-get install -y zstd wget git libtool automake build-essential

RUN mkdir /usr/work
WORKDIR /usr/work

RUN cd /usr/work && git clone https://github.com/smhafiz/private_queries_it_pir.git .
RUN wget https://gmplib.org/download/gmp/gmp-6.3.0.tar.zst && tar -xvf gmp-6.3.0.tar.zst
RUN cd gmp-6.3.0 && ./configure && make && make install && make check
RUN wget https://libntl.org/ntl-11.5.1.tar.gz && tar -xvf ntl-11.5.1.tar.gz
RUN ldconfig && cd ntl-11.5.1/src && ./configure && make && make install

RUN git clone https://github.com/hstraub/socketxx
RUN apt-get install -y texinfo
RUN cd socketxx && ./autogen && ./configure && make && make install

RUN apt-get install -y libgcrypt-dev
COPY percy++-1.0.0 /usr/work/percy 
RUN ldconfig && cd percy && make

RUN apt-get install -y zip tmux
RUN cp percy/itserver.cc itserver.cc

#CMD tail -f /dev/null
