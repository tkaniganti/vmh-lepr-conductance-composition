# env.sh - source this at the top of every script and in interactive shells.
# Usage:  source /ocean/projects/bio260082p/kanigant/vmh-leptin-project/env.sh

export ACCOUNT=bio260082p
# PROJECT defaults to the directory this file lives in, so a fresh clone runs
# unmodified. Set PROJECT in the environment to override.
_envsh="${BASH_SOURCE[0]:-$0}"
_guess="$(cd "$(dirname "$_envsh")" 2>/dev/null && pwd)"
[ -f "$_guess/vmhlib.py" ] || _guess="${SLURM_SUBMIT_DIR:-$PWD}"
[ -f "$_guess/vmhlib.py" ] || _guess="$PWD"
export PROJECT="${PROJECT:-$_guess}"
unset _guess
unset _envsh
if [ ! -f "$PROJECT/vmhlib.py" ]; then
    echo "env.sh: PROJECT=$PROJECT does not contain vmhlib.py." >&2
    echo "Set PROJECT explicitly to the checkout directory." >&2
    return 1 2>/dev/null || exit 1
fi

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
