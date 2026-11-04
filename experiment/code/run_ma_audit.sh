#!/bin/bash
# Background driver for the P0 audit (one process).
set -u
C=$(cd "$(dirname "$0")" && pwd); PY=${PY:-$C/../../review/repro/venv/bin/python}
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2 PYTHONDONTWRITEBYTECODE=1 MPLBACKEND=Agg
cd "$C"; $PY -B ma_audit.py --reps 200 --thin-outer 100 --thin-inner 20 --grid-reps 100 --sim-reps 200 > ../logs/rev2_audit.stdout 2>&1
echo "rc=$?" >> ../logs/rev2_audit.stdout
