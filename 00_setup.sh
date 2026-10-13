#!/usr/bin/env bash
# 00_setup.sh - build the conda env. RUN ON A LOGIN NODE, IN THE FOREGROUND, INSIDE tmux.
# Compute nodes have no outbound network, and backgrounded jobs get reaped on login nodes.
# Small batches, one at a time, as on this account.
set -euo pipefail
cd "$(dirname "$0")"
source ./env.sh

echo "== creating env at $ENVPATH =="
if [ ! -x "$PY" ]; then
  conda create -y -p "$ENVPATH" -c conda-forge python=3.11
fi
conda activate "$ENVPATH"

echo "== batch 1: numerical core =="
conda install -y -p "$ENVPATH" -c conda-forge numpy scipy pandas

echo "== batch 2: stats and plotting =="
conda install -y -p "$ENVPATH" -c conda-forge statsmodels scikit-learn matplotlib

echo "== batch 3: pip extras =="
# diptest is small and has no conda-forge build on every platform; pip is fine here.
"$ENVPATH/bin/pip" install --no-input diptest || \
  echo "WARNING: diptest unavailable; 04_discreteness.py will fall back to GMM/BIC only"

echo "== verify =="
"$PY" -c "import numpy,scipy,pandas,statsmodels,sklearn,matplotlib;print('core ok')"
"$PY" -c "import diptest;print('diptest ok')" || true
echo "Setup complete. ENVPATH=$ENVPATH"
