FROM ubuntu:22.04

RUN apt-get update
RUN apt-get install -y \
    gcc \
    g++ \
    gdb \
    valgrind \
    libboost-all-dev \
    openmpi-bin \
    libopenmpi-dev \
    git \
    python3 \
    python3-pip \
    libfmt-dev \
    libtbb-dev \
    mc
RUN apt-get install -y catch2  cmake
RUN apt-get install -y vim curl wget sudo 
RUN apt-get install -y clang-format clang-tidy clangd

RUN groupadd -g 1000 user && \
    useradd -u 1000 -g 1000 -m -s /bin/bash user && \
    echo "user ALL=(ALL) NOPASSWD:ALL" >> /etc/sudoers && \
    echo "Set disable_coredump false" >> /etc/sudo.conf

USER user
WORKDIR /app
CMD tail -f /dev/null
