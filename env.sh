# env.sh - source this at the top of every script and in interactive shells.
# Usage:  source /ocean/projects/bio260082p/kanigant/vmh-leptin-project/env.sh

export ACCOUNT=bio260082p
export PROJECT=/ocean/projects/bio260082p/kanigant/vmh-leptin-project

# conda base is the shared system install; the account has no `conda init`,
# so source the hook and activate by FULL PREFIX PATH, never by name.
export CONDA_HOOK=/opt/packages/anaconda3-2024.10-1/etc/profile.d/conda.sh
export ENVPATH=/ocean/projects/bio260082p/kanigant/conda/envs/vmhlepr
export CONDA_PKGS_DIRS=/ocean/projects/bio260082p/kanigant/conda/pkgs

# Always call python by full path. Do not rely on `python` or `python3` on PATH.
export PY=$ENVPATH/bin/python

export DATA=$PROJECT/data
export RESULTS=$PROJECT/results
export FIGURES=$PROJECT/figures
export LOGS=$PROJECT/logs

if [ -f "$CONDA_HOOK" ]; then
  # shellcheck disable=SC1090
  . "$CONDA_HOOK"
  conda activate "$ENVPATH" 2>/dev/null || true
fi

mkdir -p "$DATA" "$RESULTS" "$FIGURES" "$LOGS"
