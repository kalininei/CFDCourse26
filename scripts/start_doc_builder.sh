#!/bin/bash
if ! docker image inspect cfd26-doc-image &> /dev/null; then
    echo "Building cfd26-doc-image"
    docker build doc -t cfd26-doc-image
fi
docker run -it --rm --name cfd26-doc -v "$(pwd):/workspace" cfd26-doc-image /bin/bash
