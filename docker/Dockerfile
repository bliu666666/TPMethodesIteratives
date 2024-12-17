FROM ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive

RUN apt update && apt install -y \
    build-essential \
    make \
    gcc \
    gfortran \
    libblas-dev \
    liblapack-dev \
    liblapacke-dev \
 && rm -rf /var/lib/apt/lists/*

WORKDIR /Poisson1D

COPY . /Poisson1D

RUN mkdir -p bin

ENV HOSTNAME=lbw-virtual-machine

RUN make all

CMD ["make", "run_tpPoisson1D_direct"]
