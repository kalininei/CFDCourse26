FROM ubuntu:24.04

# ================ INSTALL DEPS
RUN apt-get update  && \
    apt-get install -y gcc-14  g++-14 cmake gdb valgrind  && \
    update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-14 100 && \
    update-alternatives --install /usr/bin/g++ g++ /usr/bin/g++-14 100 && \
    update-alternatives --install /usr/bin/gcov gcov /usr/bin/gcov-14 100 && \
    apt-get install -y libboost-dev  && \
    apt-get install -y openmpi-bin libopenmpi-dev  && \
    apt-get install -y git  && \
    apt-get install -y python3 python3-pip  && \
    apt-get install -y libfmt-dev libtbb-dev catch2  && \
    apt-get install -y mc vim curl wget sudo   && \
    apt-get install -y clang-format  && \
    apt-get install -y iputils-ping  && \
    apt-get install -y rsync  && \
    apt-get install -y libsuitesparse-dev libxml2-dev && \
    pip3 install --break-system-packages --no-cache-dir black && \
    pip3 install --break-system-packages --no-cache-dir cmake-format  && \
    pip3 install --break-system-packages --no-cache-dir decorator

# ================ SET USER
RUN usermod -l user ubuntu && \
    groupmod -n user ubuntu && \
    usermod -d /home/user -m user && \
    echo "user ALL=(ALL) NOPASSWD:ALL" >> /etc/sudoers

ARG HOST_UID=1000
ARG HOST_GID=1000

RUN if [ "$HOST_UID" != "1000" ]; then usermod -u $HOST_UID user; fi && \
    if [ "$HOST_GID" != "1000" ]; then groupmod -g $HOST_GID user; fi

USER user

# ================ misc
RUN echo 'export CC=/usr/bin/gcc-14' >> /home/user/.bashrc  && \
    echo 'export CXX=/usr/bin/g++-14' >> /home/user/.bashrc && \
    mkdir -p /home/user/.vscode-server && \
    printf '\n\
# run mc remembering last mc directory on exit\n\
function mcc() {\n\
    MC_TMP_FILE=~/.mclast.txt\n\
    /usr/bin/mc --printwd="$MC_TMP_FILE"\n\
    if [ $? -eq 0 ] && [ -f "$MC_TMP_FILE" ]; then\n\
        NEW_DIR=$(cat "$MC_TMP_FILE" 2>/dev/null)\n\
        if [ -d "$NEW_DIR" ] && [ "$NEW_DIR" != "$ORIG_DIR" ]; then\n\
            echo $NEW_DIR\n\
            cd "$NEW_DIR" || return\n\
        fi\n\
    fi\n\
    rm -f "$MC_TMP_FILE"\n\
}\n' >> ~/.bashrc

CMD /app/scripts/init/init_container.sh
