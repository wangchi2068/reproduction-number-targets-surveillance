#!/bin/bash
# Authors' pipeline (publication budgets as in run_revision.py) rerun on the copied
# code/data, one sequential lane. The design stage is not run: Table S7 is withdrawn
# (see findings.md). Calibration is not rerun; its reproduced outputs are copied.
set -u
C=$(cd "$(dirname "$0")" && pwd); PY=${PY:-$C/../../review/repro/venv/bin/python}
O=$C/../analysis/authors; L=$C/../logs; mkdir -p "$O" "$L"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 PYTHONIOENCODING=utf-8 MPLBACKEND=Agg
cd "$C"
run(){ name=$1; shift; echo "$(date -u '+%F %T') start $name" >> $L/authors_timeline.log; $PY -B "$@" > $L/authors_$name.log 2>&1; echo "$(date -u '+%F %T') end $name rc=$?" >> $L/authors_timeline.log; }
$PY -B -c "from revision_data import audit_inputs; from pathlib import Path; audit_inputs(Path('$O'))" > $L/authors_audit.log 2>&1
run empirical revision_empirical.py --output-dir $O --bootstrap 2000 --secondary-bootstrap 1000
run fixed_effects s03c_state_fe.py --output-dir $O --bootstrap 1000
run fixed_effects_block2 s03c_state_fe.py --output-dir $O/block2 --bootstrap 1000 --block-length 2
run fixed_effects_block3 s03c_state_fe.py --output-dir $O/block3 --bootstrap 1000 --block-length 3
run theory revision_theory.py --output-dir $O --replicates 200000
run secondary revision_secondary.py --output-dir $O --bootstrap 1000
run thinning s03b_thinning.py --output-dir $O --bootstrap 1000 --point-reps 400 --inner-reps 50 --max-inner-reps 400
mkdir -p $O/calibration_reproduced
cp ../../review/repro/out/calibration*.json $O/calibration_reproduced/
echo "copied from review/repro/out (reproduce-check rerun, bit-identical to authors' committed calibration)" > $O/calibration_reproduced/PROVENANCE.txt
echo "$(date -u '+%F %T') all done" >> $L/authors_timeline.log
