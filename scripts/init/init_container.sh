#!/bin/bash
mkdir -p /app/.vscode
mkdir -p /app/.vscode-server

sudo chown user:user /app/.vscode
sudo chown user:user /app/.vscode-server

rsync -a --ignore-existing /app/scripts/init/.vscode/*        /app/.vscode
rsync -a --ignore-existing /app/scripts/init/.vscode-server/* /app/.vscode-server

exec tail -f /dev/null
