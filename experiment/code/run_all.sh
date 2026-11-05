#!/bin/bash
# Full pipeline: every number, table and figure of the manuscript from the data in ../data.
#   bash code/run_all.sh            (from the experiment directory, or from anywhere)
#   PY=/path/to/python bash code/run_all.sh
#   SKIP_RENDER=1 bash code/run_all.sh   (analyses only; then python code/ma_render.py)
# Four lanes run in parallel on at most 7 cores, then the renderer, the test table and the
# gap-closing analyses (interval calibration, share-grid tests, held-out reference). Any failing
# step stops the run with a non-zero exit code; logs are in ../logs/run_*.log.
set -euo pipefail
C=$(cd "$(dirname "$0")" && pwd); E=$(dirname "$C"); L=$E/logs; O=$E/analysis
PY=${PY:-$E/venv/bin/python}; [ -x "$PY" ] || PY=$(command -v python3)
export PYTHONDONTWRITEBYTECODE=1 PYTHONIOENCODING=utf-8 MPLBACKEND=Agg
export MA_EXT_LOG=run_extensions_progress.log MA_FE_LOG=run_fe_decomp_progress.log
mkdir -p "$L" "$O/authors"; cd "$C"
T=$L/run_timeline.log
step(){ local name=$1 thr=$2; shift 2
  echo "$(date -u '+%F %T') start $name" >> "$T"
  if OMP_NUM_THREADS=$thr OPENBLAS_NUM_THREADS=$thr MKL_NUM_THREADS=$thr "$PY" -B "$@" > "$L/run_$name.log" 2>&1; then
    echo "$(date -u '+%F %T') end $name rc=0" >> "$T"
  else rc=$?; echo "$(date -u '+%F %T') FAIL $name rc=$rc" >> "$T"; return $rc; fi; }

lane_authors(){  # authors' estimators at publication budgets (design stage and calibration withdrawn)
  step authors_audit 1 -c "from revision_data import audit_inputs; from pathlib import Path; audit_inputs(Path('$O/authors'))"
  step authors_empirical 1 revision_empirical.py --output-dir "$O/authors" --bootstrap 2000 --secondary-bootstrap 1000
  step authors_fixed_effects 1 s03c_state_fe.py --output-dir "$O/authors" --bootstrap 1000
  step authors_fixed_effects_block2 1 s03c_state_fe.py --output-dir "$O/authors/block2" --bootstrap 1000 --block-length 2
  step authors_fixed_effects_block3 1 s03c_state_fe.py --output-dir "$O/authors/block3" --bootstrap 1000 --block-length 3
  step authors_theory 1 revision_theory.py --output-dir "$O/authors" --replicates 200000
  step authors_secondary 1 revision_secondary.py --output-dir "$O/authors" --bootstrap 1000
  step authors_thinning 1 s03b_thinning.py --output-dir "$O/authors" --bootstrap 1000 --point-reps 400 --inner-reps 50 --max-inner-reps 400
  step thinning_exact 1 ma_thinning_exact.py; }
lane_audit(){ step audit 2 ma_audit.py --reps 200 --thin-outer 100 --thin-inner 20 --grid-reps 100 --sim-reps 200 --log-name run_audit_progress.log
  step china 2 china_covid_empirical.py; }
lane_empirical(){ rm -f "$O/ma_extensions.json"
  step extensions 1 ma_extensions.py table1 ratios within sensitivity seasons validation
  step fe_decomp 1 ma_fe_decomp.py; }
lane_extra(){ step audit_extra 1 ma_audit_extra.py --reps 200 --workers 3 --log-name run_audit_extra_progress.log
  step coverage 1 ma_coverage.py --reps 400 --boot 200 --workers 3 --log-name run_coverage_progress.log
  step audit_national 1 ma_audit_national.py --workers 3 --log-name run_audit_national_progress.log; }

echo "$(date -u '+%F %T') run_all start (python $PY)" >> "$T"
lane_authors & P1=$!; lane_audit & P2=$!; lane_empirical & P3=$!; lane_extra & P4=$!
rc=0; for p in $P1 $P2 $P3 $P4; do wait "$p" || rc=1; done
[ $rc -eq 0 ] || { echo "$(date -u '+%F %T') run_all FAILED in a lane (see run_timeline.log)" >> "$T"; exit 1; }
[ -n "${SKIP_RENDER:-}" ] || { step render 1 ma_render.py; step tests 1 ma_tests.py; }
# Gap-closing analyses (need the coverage panels, ma_extensions.json and the tests/ entries of results.json); 7 workers
step gap_calibration 1 ma_calibration.py --workers 7
[ -n "${SKIP_RENDER:-}" ] || step gap_share_tests 1 ma_share_tests.py --workers 7
step gap_heldout_sim 1 ma_heldout_sim.py --workers 7
[ -n "${SKIP_RENDER:-}" ] || step gaps_merge 1 ma_gaps.py
echo "$(date -u '+%F %T') run_all done" >> "$T"
