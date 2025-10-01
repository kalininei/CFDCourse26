#!/bin/bash
sudo chown -R user:user /app

mkdir -p /app/.vscode
mkdir -p /app/.vscode-server

rsync -a --ignore-existing /app/scripts/init/.vscode/*        /app/.vscode
rsync -a --ignore-existing /app/scripts/init/.vscode-server/* /app/.vscode-server

sudo chmod -R 777 /app/.vscode
sudo chmod -R 777 /app/.vscode-server

# HybMesh
HYBMESH_BIN="/usr/local/bin/hybmesh"
HYBMESH_STAB_BIN="/app/scripts/init/hybmesh.stab"
if [ ! -f "$HYBMESH_BIN" ]; then
    echo "HybMesh not found at $HYBMESH_BIN"
    sudo cp "$HYBMESH_STAB_BIN" "$HYBMESH_BIN"
    echo "Created HybMesh stab script"
fi

echo "init_container.sh - DONE"
exec tail -f /dev/null
