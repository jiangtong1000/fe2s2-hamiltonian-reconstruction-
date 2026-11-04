#!/usr/bin/env bash
# Source your site's chemistry environment only when explicitly configured.
if [[ -n "${FE2S2_RUNTIME:-}" ]]; then
    set +u
    source "$FE2S2_RUNTIME"
    set -u
fi
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-${OMP_NUM_THREADS:-1}}"
export MKL_NUM_THREADS="$OMP_NUM_THREADS"
export OPENBLAS_NUM_THREADS=1
export PYTHONUNBUFFERED=1
export MPLBACKEND=Agg
