#!/usr/bin/env bash
#SBATCH --job-name=fe2s2_step
#SBATCH --cpus-per-task=16
#SBATCH --mem=160G
#SBATCH --time=12:00:00
#SBATCH --output=fe2s2-%j.out
#SBATCH --error=fe2s2-%j.err
set -euo pipefail
# Submit from the repository root; Slurm copies this script to a spool directory.
cd "${SLURM_SUBMIT_DIR:?Submit this script using sbatch from the repository root}"
source environment/runtime.sh
step="${1:?Usage: sbatch environment/submit_step.sh 02_parent_uks --workdir /scratch/run}"
shift
[[ "$step" =~ ^[0-9][0-9]_[a-z0-9_]+$ && -f "$step/run.py" ]] || { echo "Unknown step: $step" >&2; exit 2; }
exec python "$step/run.py" "$@"
