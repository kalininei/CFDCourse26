FROM texlive/texlive:latest

# Устанавливаем утилиты
RUN apt-get update && apt-get install -y \
    make \
    python3 \
    python3-pip \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /workspace

CMD ["/bin/bash"]
