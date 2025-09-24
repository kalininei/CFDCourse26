#!/usr/bin/bash

SRC_DIR=/app/extern/~hybmesh_build_dir
HM_ARC=/app/extern/HybMesh.tar.gz

# remember current directory
ORIG_DIR=$(pwd)

# Check if archive exists
if [ ! -f "$HM_ARC" ]; then
    echo "Error: Archive $HM_ARC not found!"
    exit 1
fi

# if $SRC_DIR does not exist -> unpack HM_ARC
if [ ! -d "$SRC_DIR" ]; then
    echo "Extracting $HM_ARC to $SRC_DIR"
    mkdir -p "$SRC_DIR"
    tar -xzf "$HM_ARC" -C "$SRC_DIR" --strip-components=1
fi

# go to hybmesh_build_dir
cd "$SRC_DIR" || { echo "Failed to enter $SRC_DIR"; exit 1; }

# build procedure
mkdir -p build
cd build || { echo "Failed to enter build directory"; exit 1; }
cmake .. -DCMAKE_BUILD_TYPE=Release || { echo "CMake failed"; exit 1; }
make -j8
sudo make install || { echo "Make install failed"; exit 1; }
sudo chown -R user:user "$SRC_DIR"

# go back to the original directory
cd "$ORIG_DIR"
echo "Build completed successfully. Back in: $(pwd)"
