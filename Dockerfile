FROM ubuntu:24.04

RUN apt-get update
RUN apt-get install -y gcc-14  g++-14 cmake gdb valgrind
RUN apt-get install -y libboost-dev
RUN apt-get install -y openmpi-bin libopenmpi-dev
RUN apt-get install -y git
RUN apt-get install -y python3 python3-pip
RUN apt-get install -y libfmt-dev libtbb-dev catch2
RUN apt-get install -y mc vim curl wget sudo 
RUN apt-get install -y clang-format clang-tidy clangd
RUN apt-get install -y iputils-ping

RUN pip3 install --break-system-packages --no-cache-dir black
RUN pip3 install --break-system-packages --no-cache-dir cmake-format

RUN update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-14 100 && \
    update-alternatives --install /usr/bin/g++ g++ /usr/bin/g++-14 100 && \
    update-alternatives --install /usr/bin/gcov gcov /usr/bin/gcov-14 100

RUN usermod -l user ubuntu && \
    groupmod -n user ubuntu && \
    usermod -d /home/user -m user && \
    echo "user ALL=(ALL) NOPASSWD:ALL" >> /etc/sudoers

USER user

RUN echo 'export CC=/usr/bin/gcc-14' >> /home/user/.bashrc
RUN echo 'export CXX=/usr/bin/g++-14' >> /home/user/.bashrc 

RUN mkdir -p /home/user/.vscode-server

WORKDIR /app
CMD tail -f /dev/null
